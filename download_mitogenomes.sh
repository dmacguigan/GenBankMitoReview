#!/bin/bash
# Download mitogenomes from GenBank for a group of organisms into a .gb file.
#
# author: Dan MacGuigan
#
# Run inside the 'mitoreview' conda/mamba environment (see environment.yml):
#   mamba activate mitoreview
#
# Dependencies (provided by that environment):
#   Entrez Direct tools (esearch, efilter, efetch)
#   python + biopython (for filter_gb.py)
#
# The query is an advanced GenBank search term, e.g. a taxon ("Percidae"[Organism])
# or a BioProject ("PRJNA720393"[BioProject]); see
# https://www.ncbi.nlm.nih.gov/nuccore/advanced
# Multiple terms can be combined, e.g. '"Percidae"[Organism] AND "PRJNA720393"[BioProject]'
#
# By default this downloads all mitochondrial genomic records matching the query,
# then drops records with one or fewer genes, leaving multi-gene mitogenomes.
#
# Usage:
#   bash download_mitogenomes.sh [--complete] [--refseq] '<query>' <output_prefix>
#
#   --complete   restrict to records titled "complete genome"
#   --refseq     restrict to RefSeq records (NC_* accessions)
#
# Example:
#   bash download_mitogenomes.sh '"Percidae"[Organism]' percidae
#   bash download_mitogenomes.sh --complete '"Percidae"[Organism]' percidae
#   bash download_mitogenomes.sh --refseq '"Percidae"[Organism]' percidae
#
# (NOTE: the query must be wrapped in single quotes)

set -euo pipefail

# locate this script's directory so filter_gb.py is found regardless of cwd
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Prepend the bundled curl shim so Entrez Direct uses HTTP/1.1 instead of its
# HTTP/1.0 default, which fails with "curl (56) unexpected eof" on OpenSSL 3.x.
PATH="${SCRIPT_DIR}/bin:${PATH}"

print_help() {
  cat <<'EOF'
Usage:
  bash download_mitogenomes.sh [OPTIONS] '<query>' <output_prefix>

Description:
  Downloads mitochondrial genomic records from NCBI GenBank matching the given
  query and writes them to <output_prefix>.gb, retaining only multi-gene records
  (i.e. true mitogenomes). A raw unfiltered copy is saved as <output_prefix>.raw.gb.

  The query is an advanced GenBank search term. It must be wrapped in single
  quotes. See https://www.ncbi.nlm.nih.gov/nuccore/advanced for query syntax.

Arguments:
  '<query>'          Advanced GenBank search term (single-quoted), e.g.
                       '"Percidae"[Organism]'
                       '"Percidae"[Organism] AND "PRJNA720393"[BioProject]'
  <output_prefix>    Prefix for all output files (e.g. 'percidae' produces
                     percidae.raw.gb and percidae.gb)

Options:
  --outdir <dir>     Directory for all output files (default: current directory).
                     Created automatically if it does not exist. Accepts full or
                     relative paths.
  --api-key <key>    NCBI API key. Raises the Entrez request rate limit from
                     3 to 10 requests/second. Register at
                     https://www.ncbi.nlm.nih.gov/account/
  --complete         Restrict to records whose title contains "complete genome"
  --refseq           Restrict to RefSeq records (NC_* accessions)
  -h, --help         Show this help message and exit

Output:
  <outdir>/<prefix>.raw.gb    All matching records before gene-count filtering
  <outdir>/<prefix>.gb        Multi-gene mitogenomes only (records with >1 gene feature)

Examples:
  bash download_mitogenomes.sh '"Percidae"[Organism]' percidae
  bash download_mitogenomes.sh --outdir results '"Percidae"[Organism]' percidae
  bash download_mitogenomes.sh --outdir /data/mito --complete '"Percidae"[Organism]' percidae

Dependencies (provided by the 'mitoreview' mamba environment):
  Entrez Direct (esearch, efilter, efetch), python, biopython
EOF
}

# INPUT / OPTIONS
COMPLETE=0
REFSEQ=0
OUTDIR="."
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)  print_help; exit 0 ;;
    --complete) COMPLETE=1; shift ;;
    --refseq)   REFSEQ=1; shift ;;
    --outdir)   [[ $# -ge 2 ]] || { echo "ERROR: --outdir requires a value" >&2; exit 1; }
                OUTDIR="$2"; shift 2 ;;
    --api-key)  [[ $# -ge 2 ]] || { echo "ERROR: --api-key requires a value" >&2; exit 1; }
                export NCBI_API_KEY="$2"; shift 2 ;;
    --) shift; break ;;
    -*) echo "ERROR: unknown option '$1'" >&2; echo "Use -h for help." >&2; exit 1 ;;
    *)  break ;;
  esac
done

if [[ $# -ne 2 ]]; then
  echo "Usage: bash download_mitogenomes.sh [--outdir DIR] [--api-key KEY] [--complete] [--refseq] '<query>' <output_prefix>" >&2
  echo "Use -h for full help." >&2
  exit 1
fi

USER_QUERY="$1"
PREFIX="$2"

# create output directory if needed and build the full file prefix
mkdir -p "${OUTDIR}"
OUTDIR="${OUTDIR%/}"   # strip any trailing slash
FULL_PREFIX="${OUTDIR}/${PREFIX}"

# build the GenBank query
QUERY="(${USER_QUERY})"
if [[ "${COMPLETE}" -eq 1 ]]; then
  QUERY="${QUERY} AND \"complete genome\"[Title]"
fi
if [[ "${REFSEQ}" -eq 1 ]]; then
  QUERY="${QUERY} AND refseq[filter]"
fi

echo "Downloading mitochondrial genomic records from NCBI GenBank using the query:"
echo "  ${QUERY}"
if [[ -n "${NCBI_API_KEY:-}" ]]; then
  echo "  (NCBI API key set — 10 requests/second limit)"
fi
echo ""

# count first so the user knows the scale
N_REC=$( esearch -db nuccore -query "${QUERY}" \
  | efilter -location mitochondrion -molecule genomic \
  | grep "Count" | head -1 | cut -f2 -d">" | cut -f1 -d"<" )

echo "Found ${N_REC} matching mitochondrial record(s). Downloading now..."
echo "This may take a few minutes."
echo ""

esearch -db nuccore -query "${QUERY}" \
  | efilter -location mitochondrion -molecule genomic \
  | efetch -format gb > "${FULL_PREFIX}.raw.gb"

echo "Raw records written to ${FULL_PREFIX}.raw.gb"
echo "Filtering to multi-gene mitogenomes..."
echo ""

# drop records with <=1 gene; final mitogenome set is <full_prefix>.gb
python "${SCRIPT_DIR}/filter_gb.py" "${FULL_PREFIX}.raw.gb" "${FULL_PREFIX}.gb"

echo ""
echo "Done. Mitogenome records written to ${FULL_PREFIX}.gb"
echo "Summarize codon usage with:"
echo "  python ${SCRIPT_DIR}/summarize_codons.py ${FULL_PREFIX}.gb -o ${FULL_PREFIX}"
