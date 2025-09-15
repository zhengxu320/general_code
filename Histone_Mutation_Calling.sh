#!/usr/bin/env bash
set -euo pipefail

#!/usr/bin/env bash
set -euo pipefail

# wes_histone_pipeline_min.sh
# Minimal, modular WES pipeline (150bp PE):
#   FastQC -> BWA MEM -> sort -> MarkDuplicates -> GATK HaplotypeCaller -> histone-region VCF
#
# Requirements in PATH: fastqc, bwa, samtools, gatk (v4+), bcftools, bedtools
#
# HISTONE BED (you provide; required):
#   - Standard BED (0-based, half-open)
#   - Columns: chrom  start  end  [optional: gene_name]
#   - Reference naming must match your FASTA (“chr1” vs “1”)
#   - Exon-only or merged gene spans are both fine; we just intersect VCF with this BED
#
# Example:
#   chr1    226062282   226063281   HIST1H3A
#   chr1    226073281   226074050   HIST1H4A
#
# Example usage:
#   bash wes_histone_pipeline_min.sh \
#     --fq1 S1_R1.fastq.gz --fq2 S1_R2.fastq.gz \
#     --sample S1 --ref GRCh38.fa \
#     --exome-bed targets/exome.bed \
#     --histone-bed annotations/histones.hg38.bed \
#     --threads 16 --mem-gb 32 --outdir results

############################################
# Tunables (override via CLI)
############################################
THREADS=8
MEM_GB=24                # Java heap for GATK steps
OUTDIR="results"

SAMPLE=""
FQ1=""
FQ2=""
REF_FA=""
EXOME_BED=""            # optional capture targets for HaplotypeCaller (-L)
HISTONE_BED=""          # required
PLATFORM="ILLUMINA"
LIB="lib1"
UNIT="unit1"
RGID=""                 # defaults to SAMPLE if empty
FORCE=0                 # overwrite outputs if 1

############################################
# Tiny helpers
############################################
log(){ echo "[$(date '+%F %T')] $*" >&2; }
die(){ echo "ERROR: $*" >&2; exit 1; }
ensure_dir(){ mkdir -p "$@"; }
java_mem(){ echo "-Xmx${MEM_GB}g"; }
check_file(){ [[ -r "$1" ]] || die "Cannot read $1"; }

run_if_needed(){ # target, then command string
  local target="$1"; shift
  if [[ -s "$target" && "$FORCE" -eq 0 ]]; then
    log "Exists (skip): $target"
  else
    log "Running: $*"
    bash -c "$*"
    [[ -s "$target" ]] || die "Expected output missing: $target"
  fi
}

check_ref_indexes(){
  local fa="$1"; check_file "$fa"
  [[ -s "${fa}.fai" ]] || { log "samtools faidx $fa"; samtools faidx "$fa"; }
  local dict="${fa%.fa}.dict"; dict="${dict%.fasta}.dict"
  [[ -s "$dict" ]] || { log "gatk CreateSequenceDictionary"; gatk $(java_mem) CreateSequenceDictionary -R "$fa" -O "$dict"; }
}

############################################
# Parse CLI
############################################
while [[ $# -gt 0 ]]; do
  case "$1" in
    --fq1) FQ1="$2"; shift 2;;
    --fq2) FQ2="$2"; shift 2;;
    --sample) SAMPLE="$2"; shift 2;;
    --ref) REF_FA="$2"; shift 2;;
    --exome-bed) EXOME_BED="$2"; shift 2;;
    --histone-bed) HISTONE_BED="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --mem-gb) MEM_GB="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --platform) PLATFORM="$2"; shift 2;;
    --lib) LIB="$2"; shift 2;;
    --unit) UNIT="$2"; shift 2;;
    --rgid) RGID="$2"; shift 2;;
    --force) FORCE=1; shift;;
    -h|--help) grep -n "^# " "$0" | sed 's/^# //'; exit 0;;
    *) die "Unknown arg: $1";;
  esac
done

# Required inputs
[[ -n "$FQ1" && -n "$FQ2" && -n "$REF_FA" && -n "$HISTONE_BED" ]] || die "Required: --fq1 --fq2 --ref --histone-bed (and usually --sample)"
check_file "$FQ1"; check_file "$FQ2"; check_file "$HISTONE_BED"
check_ref_indexes "$REF_FA"
[[ -z "$EXOME_BED" ]] || check_file "$EXOME_BED"

if [[ -z "$SAMPLE" ]]; then
  base=$(basename "$FQ1")
  SAMPLE="${base%%_*}"
  log "Inferred SAMPLE=$SAMPLE"
fi
[[ -n "$RGID" ]] || RGID="$SAMPLE"

############################################
# Layout
############################################
QC_DIR="$OUTDIR/qc"
ALIGN_DIR="$OUTDIR/align"
VCF_DIR="$OUTDIR/vcf"
ensure_dir "$OUTDIR" "$QC_DIR" "$ALIGN_DIR" "$VCF_DIR"

############################################
# Step 1: FastQC
############################################
step_fastqc(){
  # Use a .done flag because FastQC produces *_fastqc.html per input
  local flag="$QC_DIR/.fastqc_${SAMPLE}.done"
  if [[ -s "$flag" && "$FORCE" -eq 0 ]]; then
    log "FastQC already done."
    return 0
  fi
  log "Running FastQC..."
  fastqc -t "$THREADS" -o "$QC_DIR" "$FQ1" "$FQ2"
  touch "$flag"
}

############################################
# Step 2: Align + sort (hard-pinned to bwa mem)
############################################
step_align_sort(){
  local rg="@RG\tID:${RGID}\tSM:${SAMPLE}\tPL:${PLATFORM}\tLB:${LIB}\tPU:${UNIT}"
  local bam="$ALIGN_DIR/${SAMPLE}.sorted.bam"

  # Ensure BWA index exists (safe to re-run; does nothing if present)
  bwa index "$REF_FA" >/dev/null 2>&1 || true

  run_if_needed "$bam" "
    bwa mem -t $THREADS -R '$rg' '$REF_FA' '$FQ1' '$FQ2' \
    | samtools sort -@ $THREADS -m 2G -o '$bam' -
  "
  run_if_needed "${bam}.bai" "samtools index -@ $THREADS '$bam'"
  echo "$bam"
}

############################################
# Step 3: Mark duplicates
############################################
step_markdups(){
  local inbam="$1"
  local outbam="$ALIGN_DIR/${SAMPLE}.md.bam"
  local metrics="$ALIGN_DIR/${SAMPLE}.md.metrics.txt"
  run_if_needed "$outbam" "
    gatk $(java_mem) MarkDuplicates \
      -I '$inbam' -O '$outbam' -M '$metrics' \
      --CREATE_INDEX true --VALIDATION_STRINGENCY LENIENT
  "
  echo "$outbam"
}

############################################
# Step 4: HaplotypeCaller (single-sample)
############################################
step_haplotypecaller(){
  local inbam="$1"
  local vcf="$VCF_DIR/${SAMPLE}.hc.vcf.gz"
  local restrict=""
  [[ -n "$EXOME_BED" ]] && restrict="-L '$EXOME_BED'"

  run_if_needed "$vcf" "
    gatk $(java_mem) HaplotypeCaller \
      -R '$REF_FA' -I '$inbam' $restrict \
      -O '$vcf' --native-pair-hmm-threads $THREADS
  "
  run_if_needed "${vcf}.tbi" "bcftools index -f '$vcf'"
  echo "$vcf"
}

############################################
# Step 5: Extract histone-region variants
############################################
step_extract_histone(){
  local vcf="$1"
  local hv="$VCF_DIR/${SAMPLE}.histone.vcf.gz"
  local tsv="$VCF_DIR/${SAMPLE}.histone.tsv"

  run_if_needed "$hv" "bcftools view -R '$HISTONE_BED' -Oz -o '$hv' '$vcf'"
  run_if_needed "${hv}.tbi" "bcftools index -f '$hv'"

  # Simple per-variant summary; fields may be '.' if absent
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%SAMPLE\t%GT]\t%INFO/AC\t%INFO/AN\n' "$hv" > "$tsv" || true
  echo "$hv"
}

############################################
# MAIN
############################################
log "Starting pipeline: SAMPLE=$SAMPLE"
step_fastqc
sorted_bam=$(step_align_sort)
md_bam=$(step_markdups "$sorted_bam")
vcf=$(step_haplotypecaller "$md_bam")
step_extract_histone "$vcf"
log "Done. Outputs in: $OUTDIR"