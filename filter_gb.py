#!/usr/bin/env python3
"""
Filter a GenBank file down to true multi-gene mitochondrial records that
contain at least one complete (non-partial) gene.

A record is kept only if BOTH conditions hold:
  1. More than one CDS/rRNA/tRNA feature (bare 'gene' wrapper features are not
     counted to avoid double-counting the single-gene case, where GenBank
     typically annotates both a 'gene' and a 'CDS' for the same locus).
  2. At least one CDS/rRNA/tRNA feature has exact (non-fuzzy) positions at
     both ends — i.e. no '<' or '>' in its location. This excludes records
     like partial single-gene submissions that slip through condition 1 and
     records composed entirely of truncated genes.

Example usage:
  python filter_gb.py input.raw.gb output.gb
"""

import argparse
from Bio import SeqIO
from Bio.SeqFeature import AfterPosition, BeforePosition

# Only count true product features; bare 'gene' features are wrapper
# annotations that always co-occur with CDS/rRNA/tRNA for the same locus.
CODING_FEATURE_TYPES = {"CDS", "rRNA", "tRNA"}


def count_coding(record):
    """Count CDS/rRNA/tRNA features in a record."""
    return sum(1 for f in record.features if f.type in CODING_FEATURE_TYPES)


def _position_is_fuzzy(pos):
    return isinstance(pos, (BeforePosition, AfterPosition))


def _location_is_complete(loc):
    """True if every part of a (possibly compound) location has exact ends."""
    parts = loc.parts if hasattr(loc, "parts") else [loc]
    return not any(
        _position_is_fuzzy(p.start) or _position_is_fuzzy(p.end)
        for p in parts
    )


def has_complete_gene(record):
    """True if at least one CDS/rRNA/tRNA feature has fully exact coordinates."""
    return any(
        _location_is_complete(f.location)
        for f in record.features
        if f.type in CODING_FEATURE_TYPES
    )


def filter_genbank(input_file, output_file):
    kept = []
    n_total = 0
    n_single = 0
    n_no_complete = 0

    with open(input_file, "r") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            n_total += 1
            if count_coding(record) <= 1:
                n_single += 1
            elif not has_complete_gene(record):
                n_no_complete += 1
            else:
                kept.append(record)

    with open(output_file, "w") as handle:
        SeqIO.write(kept, handle, "genbank")

    print(f"Filtered {n_total} record(s): kept {len(kept)} mitogenome(s), "
          f"dropped {n_single} with <=1 gene, "
          f"{n_no_complete} with no complete gene")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Filter a GenBank file to multi-gene mitochondrial records "
                    "with at least one complete gene.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Details:
  A record is kept if it passes BOTH checks:
    1. More than one CDS/rRNA/tRNA feature. Bare 'gene' features are excluded
       from the count because GenBank typically annotates both a 'gene' and a
       'CDS' for the same locus, which would cause single-gene records to
       appear multi-gene.
    2. At least one CDS/rRNA/tRNA feature has exact (non-fuzzy) coordinates
       at both ends (no '<' or '>' partial indicators). This excludes records
       composed entirely of truncated/partial genes.

  This script is called automatically by download_mitogenomes.sh but can also
  be run standalone on any GenBank file.

Examples:
  python filter_gb.py downloaded.raw.gb mitogenomes.gb
  python filter_gb.py input.gb output.gb
""")
    ap.add_argument("input",  help="input GenBank file to filter")
    ap.add_argument("output", help="output GenBank file (filtered records only)")
    args = ap.parse_args()
    filter_genbank(args.input, args.output)
