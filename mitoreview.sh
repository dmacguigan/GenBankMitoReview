#!/bin/bash
# End-to-end mitogenome codon review:
#   1) Download mitogenomes from GenBank for an organism group into <prefix>.gb
#   2) Summarize start/stop codon usage by gene, flagging non-canonical codons
#
# Automatically activates the 'mitoreview' mamba environment if not already active.
#
# Usage:
#   bash mitoreview.sh [OPTIONS] '<query>' <output_prefix>
#
# Example:
#   bash mitoreview.sh '"Percidae"[Organism]' percidae
#   bash mitoreview.sh --outdir results --api-key MYKEY '"Percidae"[Organism]' percidae
#
# The 'mitoreview' environment must exist before running this script.
# Create it once with:
#   mamba env create -f environment.yml

set -euo pipefail
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# If not already running inside the mitoreview environment, re-exec via
# 'mamba run' so the user does not need to activate manually.
if [[ "${CONDA_DEFAULT_ENV:-}" != "mitoreview" ]]; then
    if ! command -v mamba &>/dev/null; then
        echo "ERROR: mamba not found in PATH." >&2
        echo "Install mamba or activate the mitoreview environment manually:" >&2
        echo "  mamba activate mitoreview" >&2
        exit 1
    fi
    exec mamba run -n mitoreview bash "${BASH_SOURCE[0]}" "$@"
fi

print_help() {
  cat <<'EOF'
Usage:
  bash mitoreview.sh [OPTIONS] '<query>' <output_prefix>

Description:
  End-to-end mitogenome codon review pipeline:
    1) Downloads mitochondrial genomic records from NCBI GenBank matching the
       given query, retaining only multi-gene mitogenomes (via
       download_mitogenomes.sh + filter_gb.py).
    2) Summarizes start/stop codon usage by gene across all downloaded
       mitogenomes, flagging non-canonical codons according to each record's
       /transl_table genetic code (via summarize_codons.py).

Arguments:
  '<query>'          Advanced GenBank search term (single-quoted), e.g.
                       '"Percidae"[Organism]'
                       '"Percidae"[Organism] AND "PRJNA720393"[BioProject]'
  <output_prefix>    Base name for all output files (no path component;
                     the directory is set with --outdir)

Options:
  --outdir <dir>     Directory for all output files (default: current directory).
                     Created automatically if it does not exist. Accepts full or
                     relative paths.
  --api-key <key>    NCBI API key. Raises the Entrez request rate limit from
                     3 to 10 requests/second. Register at
                     https://www.ncbi.nlm.nih.gov/account/
  --complete         Restrict download to records titled "complete genome"
  --refseq           Restrict download to RefSeq records (NC_* accessions)
  -h, --help         Show this help message and exit

Output:
  <outdir>/<prefix>.raw.gb      Unfiltered downloaded records
  <outdir>/<prefix>.gb          Multi-gene mitogenomes after filtering
  <outdir>/<prefix>.codons.tsv  Per-CDS table: gene, codon, canonical status, flags
  <outdir>/<prefix>.summary.txt Per-gene codon summary (same as stdout)
  <outdir>/<prefix>.report.md   Markdown report with figures
  <outdir>/<prefix>.report.pdf  PDF rendered from the Markdown report
  stdout                        Per-gene codon distribution summary with [FLAG] lines

Examples:
  bash mitoreview.sh '"Percidae"[Organism]' percidae
  bash mitoreview.sh --outdir results '"Percidae"[Organism]' percidae
  bash mitoreview.sh --outdir /data/mito --refseq '"Percidae"[Organism]' percidae

See also:
  download_mitogenomes.sh -h    Download step help
  python summarize_codons.py -h Summarize step help
  python filter_gb.py -h        Filter step help
EOF
}

# parse options
OUTDIR="."
API_KEY=""
DOWNLOAD_OPTS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)  print_help; exit 0 ;;
    --complete) DOWNLOAD_OPTS+=("--complete"); shift ;;
    --refseq)   DOWNLOAD_OPTS+=("--refseq"); shift ;;
    --outdir)   [[ $# -ge 2 ]] || { echo "ERROR: --outdir requires a value" >&2; exit 1; }
                OUTDIR="$2"; shift 2 ;;
    --api-key)  [[ $# -ge 2 ]] || { echo "ERROR: --api-key requires a value" >&2; exit 1; }
                API_KEY="$2"
                DOWNLOAD_OPTS+=("--api-key" "$2"); shift 2 ;;
    --) shift; break ;;
    -*) echo "ERROR: unknown option '$1'" >&2; echo "Use -h for help." >&2; exit 1 ;;
    *)  break ;;
  esac
done

if [[ $# -ne 2 ]]; then
  echo "Usage: bash mitoreview.sh [--outdir DIR] [--api-key KEY] [--complete] [--refseq] '<query>' <output_prefix>" >&2
  echo "Use -h for full help." >&2
  exit 1
fi

QUERY="$1"
PREFIX="$2"
OUTDIR="${OUTDIR%/}"
FULL_PREFIX="${OUTDIR}/${PREFIX}"

bash "${SCRIPT_DIR}/download_mitogenomes.sh" "${DOWNLOAD_OPTS[@]}" --outdir "${OUTDIR}" "${QUERY}" "${PREFIX}"

echo ""
echo "Summarizing codon usage..."
echo ""
SUMMARIZE_OPTS=()
[[ -n "${API_KEY}" ]] && SUMMARIZE_OPTS+=("--api-key" "${API_KEY}")
python "${SCRIPT_DIR}/summarize_codons.py" "${FULL_PREFIX}.gb" -o "${FULL_PREFIX}" "${SUMMARIZE_OPTS[@]}"
