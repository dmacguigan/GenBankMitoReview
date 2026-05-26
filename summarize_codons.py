#!/usr/bin/env python3
"""
Summarize start and stop codon usage by gene across all mitogenomes in a
GenBank file, flagging codons that are non-canonical for the appropriate
genetic code.

Unlike a naive "first/last 3 bp of a FASTA sequence" approach, this reads the
CDS features directly from the GenBank records, so it respects:
  - strand and spliced (join) locations  -> Bio.SeqFeature.extract()
  - reading frame                         -> the /codon_start qualifier
  - the appropriate genetic code          -> the /transl_table qualifier
    (https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes)

For each CDS the start codon (first 3 bp of the coding frame) and stop codon
(last 3 bp) are compared against the start/stop codon sets defined by that
record's genetic code, via Bio.Data.CodonTable. Codons outside those sets are
flagged. Codons that cannot be judged (partial/truncated ends, common for mito
CDS whose stop is completed by polyadenylation) are flagged separately rather
than reported as non-canonical.

Example usage:
  python summarize_codons.py mitogenomes.gb -o mitogenomes
"""

import argparse
import base64
import datetime
import io
import os
import re
import sys
from collections import Counter, defaultdict

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.SeqFeature import AfterPosition, BeforePosition

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import numpy as np
    _HAVE_MPL = True
except ImportError:
    _HAVE_MPL = False

# --- protein-coding-gene name normalization ---------------------------------
# Mitochondrial PCGs appear under many synonyms in GenBank /gene and /product
# qualifiers. This maps them to a single canonical label so codon usage is
# grouped consistently across any metazoan.
#
# Ported from the MitoPilot package's normalize_pcg / pcg_from_product:
# https://github.com/Smithsonian/MitoPilot/blob/main/R/blast_ref_utils.R
# Canonical labels follow MitoPilot: nad1-6, nad4l, cox1-3, atp6, atp8, cob.

# Direct alias lookup. Keys are lower-cased; separators are also collapsed
# before lookup, so "nd-1" / "ATP 6" / "CO_I" all resolve.
PCG_LOOKUP = {
    # NADH dehydrogenase subunits
    "nad1": "nad1", "nd1": "nad1", "ndh1": "nad1", "nadh1": "nad1",
    "nad2": "nad2", "nd2": "nad2", "ndh2": "nad2", "nadh2": "nad2",
    "nad3": "nad3", "nd3": "nad3", "ndh3": "nad3", "nadh3": "nad3",
    "nad4": "nad4", "nd4": "nad4", "ndh4": "nad4", "nadh4": "nad4",
    "nad4l": "nad4l", "nd4l": "nad4l", "ndh4l": "nad4l", "nadh4l": "nad4l",
    "nad5": "nad5", "nd5": "nad5", "ndh5": "nad5", "nadh5": "nad5",
    "nad6": "nad6", "nd6": "nad6", "ndh6": "nad6", "nadh6": "nad6",
    # cytochrome c oxidase
    "cox1": "cox1", "co1": "cox1", "coi": "cox1", "coxi": "cox1",
    "cox2": "cox2", "co2": "cox2", "coii": "cox2", "coxii": "cox2",
    "cox3": "cox3", "co3": "cox3", "coiii": "cox3", "coxiii": "cox3",
    # ATP synthase
    "atp6": "atp6", "atpase6": "atp6", "atpsynthase6": "atp6",
    "atp8": "atp8", "atpase8": "atp8", "atpsynthase8": "atp8",
    # cytochrome b
    "cob": "cob", "cytb": "cob", "cyb": "cob", "cytob": "cob",
    "cytochromeb": "cob",
    # mtMutS / msh1 -- DNA mismatch-repair gene found in octocoral (and some
    # other non-bilaterian) mitogenomes, beyond the standard 13 PCGs.
    "mtmuts": "mtMutS", "muts": "mtMutS", "mutsl": "mtMutS",
    "msh1": "mtMutS", "msh": "mtMutS", "mtmsh": "mtMutS",
    "mutshomolog": "mtMutS", "mutsprotein": "mtMutS",
    # ATP synthase subunit 9 -- retained in mitogenomes of some sponges and
    # cnidarians; transferred to nucleus in most bilaterians.
    "atp9": "atp9", "atpase9": "atp9",
    # Twin-arginine translocase subunit C (tatC) -- present in some sponge
    # mitogenomes (e.g. Oscarella, Amphimedon).
    "tatc": "tatC",
    # Reverse transcriptase (rvt) -- group II intron RT domain; found in bryozoan,
    # annelid, and sponge mitogenomes as part of an intron-encoded protein (IEP).
    "rvt": "rvt", "rt": "rvt", "revtrans": "rvt",
    # Intron maturase (im) -- group II intron X/maturase domain; promotes splicing
    # of the host intron; co-occurs with rvt in cheilostome bryozoan mitogenomes.
    "im": "im", "mat": "im", "maturase": "im",
}


def _collapse(s):
    """Collapse whitespace, underscores and dashes."""
    return re.sub(r"[\s_-]+", "", s)


def pcg_from_product(text):
    """Parse a free-text product description into a canonical PCG symbol.

    Returns the symbol (e.g. 'nad4l', 'cox1', 'atp6', 'cob') or None.
    """
    p = text.lower().strip()
    if not p:
        return None

    # NADH dehydrogenase subunit X (X may be "4L")
    m = re.search(
        r"nadh[^a-z0-9]*(?:dehydrogenase)?[^a-z0-9]*(?:subunit)?[^a-z0-9]*([0-9]+l?)",
        p)
    if m:
        return "nad" + m.group(1)

    # cytochrome (c) oxidase subunit X (Roman or Arabic)
    m = re.search(
        r"cytochrome[^a-z0-9]*(?:c)?[^a-z0-9]*oxidase[^a-z0-9]*(?:subunit)?[^a-z0-9]*(i{1,3}|[0-9]+)",
        p)
    if m:
        sub_id = m.group(1)
        arabic = {"i": "1", "ii": "2", "iii": "3"}.get(sub_id, sub_id)
        return "cox" + arabic

    # (apo)cytochrome b
    if re.search(r"(?:apo)?cytochrome[^a-z0-9]*b\b", p):
        return "cob"

    # ATP synthase F0 subunit X / ATPase X. Take the LAST number after the
    # anchor so the "0" in "F0" is skipped (e.g. "ATP synthase F0 subunit 6").
    # Only subunits 6 and 8 are universally mitochondrially-encoded in metazoans.
    m = re.search(r"(?:atp[^a-z0-9]*synthase|atpase)(.*)$", p)
    if m:
        nums = re.findall(r"[0-9]+", m.group(1))
        if nums and nums[-1] in ("6", "8", "9"):
            return "atp" + nums[-1]

    # mtMutS / MutS / DNA mismatch-repair protein (octocoral mitogenomes).
    # \bmt?muts also catches the "mtmuts" form (mt- prefix not separated).
    if re.search(r"\bmt?muts\b", p) or "mismatch repair" in p:
        return "mtMutS"

    # Twin-arginine translocase subunit C (sponge mitogenomes).
    if re.search(r"twin[^a-z0-9]*arginine[^a-z0-9]*translocase", p) or \
            re.search(r"\btatc?\b", p):
        return "tatC"

    # Reverse transcriptase / group II intron RT domain (bryozoan/annelid/sponge IEP).
    if re.search(
            r"(?:group[^a-z0-9]*ii[^a-z0-9]*(?:intron[^a-z0-9]*)?)?reverse[^a-z0-9]*transcriptase",
            p):
        return "rvt"

    # Intron maturase / RNA maturase / intron-encoded protein.
    if re.search(r"(?:intron[^a-z0-9]*|rna[^a-z0-9]*)?maturase", p) or \
            re.search(r"intron[^a-z0-9]*encoded[^a-z0-9]*protein", p):
        return "im"

    return None


def normalize_pcg(name, product=None):
    """Map a raw /gene name (with /product fallback) to a canonical PCG label.

    Returns the canonical symbol or None if unrecognized. Mirrors MitoPilot's
    normalize_pcg: alias table (with separator collapse) -> product-text regex
    -> retry against the product field.
    """
    if name:
        # strip mt-/MT-/mtDNA-/m- HGNC-style prefixes
        name = re.sub(r"^(?:mt[-_ ]?dna[-_ ]|mt[-_ ]|m[-_ ])", "",
                      name.strip(), flags=re.IGNORECASE)
        key = name.lower().strip()
        if key in PCG_LOOKUP:
            return PCG_LOOKUP[key]
        if _collapse(key) in PCG_LOOKUP:
            return PCG_LOOKUP[_collapse(key)]
        result = pcg_from_product(name)
        if result:
            return result

    if product and product.strip() and product != name:
        pkey = product.lower().strip()
        if pkey in PCG_LOOKUP:
            return PCG_LOOKUP[pkey]
        if _collapse(pkey) in PCG_LOOKUP:
            return PCG_LOOKUP[_collapse(pkey)]
        result = pcg_from_product(product)
        if result:
            return result

    return None


def get_gene_label(feature):
    """Canonical PCG label for a CDS feature, plus its raw name.

    Uses /gene as the primary name and /product as fallback. Unrecognized
    genes are returned with a '?' prefix so anomalies stay visible (and grouped
    by their own name) rather than silently merged or lost.
    """
    q = feature.qualifiers
    gene = q["gene"][0] if "gene" in q else None
    product = q["product"][0] if "product" in q else None
    raw = gene or product or ""

    norm = normalize_pcg(gene or product, product)
    if norm:
        return norm, raw
    return ("?" + raw.strip()) if raw.strip() else "?unknown", raw


def get_transl_table(feature, record):
    """Genetic-code table id for a CDS: feature /transl_table, else the source
    feature's, else None (unknown)."""
    if "transl_table" in feature.qualifiers:
        return int(feature.qualifiers["transl_table"][0])
    for f in record.features:
        if f.type == "source" and "transl_table" in f.qualifiers:
            return int(f.qualifiers["transl_table"][0])
    return None


def end_is_partial(location, five_prime):
    """True if the requested end of the location is fuzzy (< or >), accounting
    for strand. five_prime=True checks the biological start of the CDS.
    On the minus strand the biological 5' end is location.end, 3' is start."""
    minus = location.strand == -1
    if five_prime:
        pos = location.end if minus else location.start
    else:
        pos = location.start if minus else location.end
    return isinstance(pos, (BeforePosition, AfterPosition))


def analyze_cds(feature, record):
    """Return a dict of analysis fields for one CDS feature, or None to skip."""
    gene, gene_raw = get_gene_label(feature)
    table_id = get_transl_table(feature, record)

    coding = feature.extract(record.seq)
    codon_start = int(feature.qualifiers.get("codon_start", ["1"])[0])
    coding = coding[codon_start - 1:]
    coding = str(coding).upper()
    length = len(coding)

    start_codon = coding[:3] if length >= 3 else coding

    start_partial = end_is_partial(feature.location, five_prime=True)
    stop_partial = end_is_partial(feature.location, five_prime=False)

    remainder = length % 3
    if remainder == 0:
        stop_codon = coding[-3:] if length >= 3 else coding
    else:
        stop_codon = coding[-remainder:]
        stop_partial = True

    # T or TA trailing bases are the start of TAA, completed by polyadenylation.
    # Treat these as complete, canonical stops rather than partial CDS.
    polya_stop = stop_codon in ("T", "TA") and stop_partial
    if polya_stop:
        stop_partial = False

    start_canonical = None
    stop_canonical = None
    table = None
    if table_id is not None:
        try:
            table = CodonTable.unambiguous_dna_by_id[table_id]
        except KeyError:
            table = None
    if table is not None:
        if not start_partial:
            start_canonical = start_codon in table.start_codons
        if not stop_partial:
            stop_canonical = stop_codon in table.stop_codons
    if polya_stop:
        stop_canonical = True  # TAA is canonical in all metazoan mito genetic codes

    notes = []
    if table_id is None:
        notes.append("no transl_table")
    elif table is None:
        notes.append(f"unknown transl_table {table_id}")
    if polya_stop:
        notes.append("stop codon completed by poly-A tail (TAA)")
    elif length % 3 != 0:
        notes.append(f"length {length} not multiple of 3")

    return {
        "record_id": record.id,
        "species": record.annotations.get("organism", ""),
        "gene": gene,
        "gene_raw": gene_raw,
        "transl_table": table_id if table_id is not None else "",
        "cds_length": length,
        "codon_start": codon_start,
        "start_codon": start_codon,
        "start_canonical": start_canonical,
        "start_flag": flag_for(start_canonical, start_partial),
        "stop_codon": stop_codon,
        "stop_canonical": stop_canonical,
        "stop_flag": "polyA" if polya_stop else flag_for(stop_canonical, stop_partial),
        "notes": "; ".join(notes),
    }


def flag_for(canonical, partial):
    """Human-readable flag for a codon."""
    if partial:
        return "partial"
    if canonical is None:
        return "unknown"
    return "ok" if canonical else "NON-CANONICAL"


def _get_taxon_id(record):
    """Return the NCBI taxon ID string from a record's source /db_xref, or None."""
    for f in record.features:
        if f.type == "source":
            for xref in f.qualifiers.get("db_xref", []):
                if xref.startswith("taxon:"):
                    return xref.split(":", 1)[1]
    return None


def _fetch_taxonomy(taxon_ids, api_key=None):
    """Fetch phylum/superclass/order/family from NCBI Taxonomy for a set of taxon IDs.

    Returns {taxon_id_str: {"phylum": ..., "superclass": ..., "order": ..., "family": ...}}.
    Missing ranks are returned as empty strings. Failures are warned and return {}.
    """
    from Bio import Entrez
    Entrez.email = "mitoreview@example.com"
    if api_key:
        Entrez.api_key = api_key

    RANKS = {"phylum", "order", "family"}
    result = {}
    ids = [str(t) for t in taxon_ids if t]
    if not ids:
        return result

    try:
        for i in range(0, len(ids), 200):
            chunk = ids[i:i + 200]
            handle = Entrez.efetch(db="taxonomy", id=",".join(chunk), retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            for rec in records:
                tid  = str(rec["TaxId"])
                info = {r: "" for r in RANKS}
                for item in rec.get("LineageEx", []):
                    rank = item["Rank"].lower()
                    if rank in RANKS:
                        info[rank] = item["ScientificName"]
                result[tid] = info
    except Exception as e:
        sys.stderr.write(f"WARNING: NCBI taxonomy fetch failed: {e}\n")

    return result


def write_tsv(rows, path):
    cols = ["record_id", "phylum", "order", "family", "species",
            "gene", "gene_raw", "transl_table",
            "cds_length", "codon_start", "start_codon", "start_canonical",
            "start_flag", "stop_codon", "stop_canonical", "stop_flag", "notes"]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) if r[c] is not None else "" for c in cols) + "\n")

def _build_summary_text(rows):
    """Return the full per-gene codon summary as a string (for stdout and files)."""
    lines = []
    w = lines.append

    table_records = defaultdict(set)
    for r in rows:
        table_records[r["transl_table"]].add(r["record_id"])

    w("=" * 80)
    w("GENETIC CODE TABLE SUMMARY")
    w("=" * 80)
    for tt in sorted(table_records,
                     key=lambda x: (x == "", int(x) if x != "" else 0)):
        recs = table_records[tt]
        if tt == "":
            w(f"\n  Table: unknown/missing   ({len(recs)} record(s))")
            continue
        try:
            tbl = CodonTable.unambiguous_dna_by_id[int(tt)]
            starts = ", ".join(sorted(tbl.start_codons))
            stops  = ", ".join(sorted(tbl.stop_codons))
            w(f"\n  Table {tt}: {tbl.names[0]}   ({len(recs)} record(s))")
            w(f"    Canonical start codons: {starts}")
            w(f"    Canonical stop codons:  {stops}")
        except (KeyError, ValueError):
            w(f"\n  Table {tt}: unrecognized   ({len(recs)} record(s))")
    w("")

    by_gene = defaultdict(list)
    for r in rows:
        by_gene[r["gene"]].append(r)

    w("=" * 80)
    w("MITOGENOME START/STOP CODON SUMMARY BY GENE")
    w("=" * 80)
    w(f"\nTotal CDS analyzed: {len(rows)}  across {len(by_gene)} gene group(s)")
    w("Canonical status uses each record's /transl_table genetic code.")
    w("'polyA' = T/TA stop codon completed to canonical TAA by polyadenylation.")
    w("'partial' = genuinely truncated/fuzzy end; not judged.\n")

    for gene in sorted(by_gene):
        cds      = by_gene[gene]
        complete = [r for r in cds
                    if r["start_flag"] != "partial" and r["stop_flag"] != "partial"]
        partial  = [r for r in cds
                    if r["start_flag"] == "partial" or r["stop_flag"] == "partial"]
        n = len(cds)

        lengths  = [r["cds_length"] for r in cds]
        mean_l   = sum(lengths) / n
        std_l    = (sum((x - mean_l) ** 2 for x in lengths) / n) ** 0.5

        w("=" * 80)
        w(f"Gene: {gene}   (n = {n} CDS, {len(complete)} complete, {len(partial)} partial)")
        w(f"  Length (bp):  min={min(lengths)}   max={max(lengths)}"
          f"   mean={mean_l:.1f}   sd={std_l:.1f}")
        w("-" * 80)

        if complete:
            nc = len(complete)
            starts   = Counter(r["start_codon"] for r in complete)
            stops    = Counter(r["stop_codon"]   for r in complete)
            nc_start = [r for r in complete if r["start_flag"] == "NON-CANONICAL"]
            nc_stop  = [r for r in complete if r["stop_flag"]  == "NON-CANONICAL"]

            w(f"  Complete CDS ({nc}):")
            w("    Start codons:")
            for codon, c in sorted(starts.items(), key=lambda x: (-x[1], x[0])):
                w(f"      {codon:<5s} {c:4d} ({100*c/nc:5.1f}%)")
            w(f"    Non-canonical start: {len(nc_start)}/{nc} ({100*len(nc_start)/nc:.1f}%)")

            w("    Stop codons:")
            for codon, c in sorted(stops.items(), key=lambda x: (-x[1], x[0])):
                tag = "  [polyA-completed]" if codon in ("T", "TA") else ""
                w(f"      {codon:<5s} {c:4d} ({100*c/nc:5.1f}%){tag}")
            n_polya = sum(1 for r in complete if r["stop_flag"] == "polyA")
            if n_polya:
                w(f"    polyA-completed stop (T/TA→TAA): {n_polya}/{nc} ({100*n_polya/nc:.1f}%)")
            w(f"    Non-canonical stop: {len(nc_stop)}/{nc} ({100*len(nc_stop)/nc:.1f}%)")

        if partial:
            w(f"  Partial CDS ({len(partial)}):")
            for r in partial:
                s = f"start {r['start_codon']} ({r['start_flag']})"
                t = f"stop  {r['stop_codon']} ({r['stop_flag']})"
                w(f"    {r['record_id']} ({r['species']}): {s}, {t}")

        w("")

    return "\n".join(lines)


def print_summary(rows, path=None):
    text = _build_summary_text(rows)
    sys.stdout.write(text + "\n")
    if path is not None:
        with open(path, "w") as fh:
            fh.write(text + "\n")


# ── figure helpers (require matplotlib) ─────────────────────────────────────

def _fig_to_file(fig, path):
    fig.savefig(path, bbox_inches="tight", dpi=150)
    plt.close(fig)


def _heatmap_fig(complete, codon_key, flag_key, title):
    """Frequency heatmap genes × codons; red border on non-canonical cells."""
    by_gene = defaultdict(list)
    for r in complete:
        by_gene[r["gene"]].append(r)
    genes      = sorted(by_gene)
    all_codons = sorted(set(r[codon_key] for r in complete))
    if not genes or not all_codons:
        return None

    mat      = np.zeros((len(genes), len(all_codons)))
    nc_cells = set()
    for i, g in enumerate(genes):
        gc    = by_gene[g]
        n     = len(gc)
        counts = Counter(r[codon_key] for r in gc)
        nc_set = {r[codon_key] for r in gc if r[flag_key] == "NON-CANONICAL"}
        for j, cod in enumerate(all_codons):
            mat[i, j] = 100 * counts.get(cod, 0) / n
            if cod in nc_set:
                nc_cells.add((i, j))

    fig, ax = plt.subplots(figsize=(max(6.0, len(all_codons) * 0.60),
                                    max(3.0, len(genes) * 0.45)))
    im = ax.imshow(mat, aspect="auto", cmap="Blues", vmin=0, vmax=100)
    plt.colorbar(im, ax=ax, label="% of complete CDS", shrink=0.8)
    ax.set_xticks(range(len(all_codons)))
    ax.set_xticklabels(all_codons, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=9)
    for i in range(len(genes)):
        for j in range(len(all_codons)):
            if mat[i, j] > 0:
                color = "white" if mat[i, j] > 60 else "black"
                ax.text(j, i, f"{mat[i, j]:.0f}", ha="center", va="center",
                        fontsize=7, color=color)
    for (i, j) in nc_cells:
        rect = mpatches.FancyBboxPatch(
            (j - 0.5, i - 0.5), 1, 1,
            boxstyle="square,pad=0", fill=False,
            edgecolor="red", linewidth=2.5, zorder=3)
        ax.add_patch(rect)
    ax.set_title(f"{title} frequency by gene  (complete CDS; red border = non-canonical)")
    fig.tight_layout()
    return fig


def _stacked_bar_fig(complete):
    """Two-panel stacked bar chart: start (top) and stop (bottom) codons by gene."""
    by_gene = defaultdict(list)
    for r in complete:
        by_gene[r["gene"]].append(r)
    genes = sorted(by_gene)
    if not genes:
        return None

    start_codons = sorted(set(r["start_codon"] for r in complete))
    stop_codons  = sorted(set(r["stop_codon"]  for r in complete))
    nc_starts    = {r["start_codon"] for r in complete if r["start_flag"] == "NON-CANONICAL"}
    nc_stops     = {r["stop_codon"]  for r in complete if r["stop_flag"]  == "NON-CANONICAL"}

    cmap   = plt.get_cmap("tab20")
    n_c    = max(len(start_codons), len(stop_codons), 1)
    sc     = {c: cmap(i / n_c) for i, c in enumerate(start_codons)}
    tc     = {c: cmap(i / n_c) for i, c in enumerate(stop_codons)}

    fig, axes = plt.subplots(2, 1,
                              figsize=(max(8.0, len(genes) * 0.7 + 2), 9),
                              sharex=True)
    x = np.arange(len(genes))

    for ax, ck, fk, codons, nc_set, colors, stitle in [
        (axes[0], "start_codon", "start_flag", start_codons, nc_starts, sc, "Start codons"),
        (axes[1], "stop_codon",  "stop_flag",  stop_codons,  nc_stops,  tc, "Stop codons"),
    ]:
        bottoms = np.zeros(len(genes))
        for cod in codons:
            freqs = np.array([
                100 * Counter(r[ck] for r in by_gene[g]).get(cod, 0) / len(by_gene[g])
                for g in genes
            ])
            label = cod + (" ★" if cod in nc_set else "")
            ax.bar(x, freqs, 0.65, bottom=bottoms,
                   color=colors[cod], hatch="//" if cod in nc_set else None,
                   edgecolor="white", linewidth=0.5, label=label)
            bottoms += freqs
        ax.set_ylim(0, 115)
        ax.set_ylabel("% complete CDS")
        ax.set_title(f"{stitle}  (// hatch = non-canonical)")
        ax.legend(loc="upper right", fontsize=7, ncol=2,
                  title="codon  (★ = non-canonical)", framealpha=0.8)

    axes[1].set_xticks(x)
    axes[1].set_xticklabels(genes, rotation=45, ha="right", fontsize=9)
    fig.tight_layout()
    return fig


def _gaussian_kde(data, x_pts, bandwidth=None):
    """Evaluate a Gaussian KDE at x_pts without scipy."""
    data = np.asarray(data, dtype=float)
    n = len(data)
    if bandwidth is None:
        std = np.std(data)
        bandwidth = max(1.06 * std * n ** -0.2, 10.0)  # Silverman, floor 10 bp
    kde = np.sum(
        np.exp(-0.5 * ((x_pts[:, None] - data[None, :]) / bandwidth) ** 2),
        axis=1,
    )
    kde /= n * bandwidth * np.sqrt(2 * np.pi)
    return kde


def _length_ridgeline_fig(rows):
    """Ridgeline (joy) plot of CDS length distributions, one ridge per gene.

    Genes are sorted by median CDS length so the plot reads bottom-to-top
    from shortest to longest gene.
    """
    by_gene = defaultdict(list)
    for r in rows:
        by_gene[r["gene"]].append(r["cds_length"])
    genes = sorted(by_gene)
    if not genes:
        return None

    all_lengths = [r["cds_length"] for r in rows]
    x_min = max(0, min(all_lengths) - 150)
    x_max = max(all_lengths) + 150
    x_pts = np.linspace(x_min, x_max, 600)

    n_genes    = len(genes)
    row_height = 1.0
    overlap    = 0.85   # KDE peak scaled to this fraction of row_height

    cmap   = plt.get_cmap("tab20")
    colors = [cmap(i / max(n_genes, 1)) for i in range(n_genes)]

    fig, ax = plt.subplots(figsize=(10, max(4, n_genes * 0.55 + 1)))

    for i, gene in enumerate(genes):
        lengths = np.array(by_gene[gene], dtype=float)
        y_base  = i * row_height
        kde     = _gaussian_kde(lengths, x_pts)
        if kde.max() > 0:
            y_vals = y_base + kde * (overlap * row_height / kde.max())
        else:
            y_vals = np.full_like(x_pts, y_base)
        color = colors[i]
        ax.fill_between(x_pts, y_base, y_vals, alpha=0.75, color=color, linewidth=0)
        ax.plot(x_pts, y_vals, color=color, linewidth=0.9)

    ax.set_yticks([i * row_height for i in range(n_genes)])
    ax.set_yticklabels(genes, fontsize=9)
    ax.set_ylim(-0.5 * row_height, n_genes * row_height)
    ax.set_xlabel("CDS length (bp)")
    ax.set_title("CDS length distribution by gene")
    for spine in ("left", "right", "top"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="y", length=0)
    fig.tight_layout()
    return fig


def _build_summary_md(rows):
    """Per-gene codon summary formatted as Markdown for embedding in the report."""
    lines = []
    w = lines.append

    table_records = defaultdict(set)
    for r in rows:
        table_records[r["transl_table"]].add(r["record_id"])

    w("## Genetic Code Summary")
    w("")
    for tt in sorted(table_records,
                     key=lambda x: (x == "", int(x) if x != "" else 0)):
        recs = table_records[tt]
        if tt == "":
            w(f"**Table: unknown/missing** — {len(recs)} record(s)")
            w("")
            continue
        try:
            tbl    = CodonTable.unambiguous_dna_by_id[int(tt)]
            starts = " · ".join(f"`{c}`" for c in sorted(tbl.start_codons))
            stops  = " · ".join(f"`{c}`" for c in sorted(tbl.stop_codons))
            w(f"**Table {tt}: {tbl.names[0]}** — {len(recs)} record(s)")
            w(f"- Canonical start codons: {starts}")
            w(f"- Canonical stop codons: {stops}")
        except (KeyError, ValueError):
            w(f"**Table {tt}: unrecognized** — {len(recs)} record(s)")
        w("")

    w("---")
    w("")
    by_gene = defaultdict(list)
    for r in rows:
        by_gene[r["gene"]].append(r)

    w("## Codon Usage by Gene")
    w("")
    w(f"*{len(rows)} CDS analyzed across {len(by_gene)} gene group(s). "
      "Canonical status uses each record's `/transl_table` genetic code.*")
    w("")

    for gene in sorted(by_gene):
        cds      = by_gene[gene]
        complete = [r for r in cds
                    if r["start_flag"] != "partial" and r["stop_flag"] != "partial"]
        partial  = [r for r in cds
                    if r["start_flag"] == "partial" or r["stop_flag"] == "partial"]
        n = len(cds)

        lengths = [r["cds_length"] for r in cds]
        mean_l  = sum(lengths) / n
        std_l   = (sum((x - mean_l) ** 2 for x in lengths) / n) ** 0.5

        w(f"### `{gene}`")
        w("")
        w(f"**{n} CDS** · {len(complete)} complete · {len(partial)} partial")
        w("")
        w(f"**Length (bp):** min={min(lengths)} · max={max(lengths)} · "
          f"mean={mean_l:.1f} · SD={std_l:.1f}")
        w("")

        if complete:
            nc       = len(complete)
            starts   = Counter(r["start_codon"] for r in complete)
            stops    = Counter(r["stop_codon"]   for r in complete)
            nc_start = [r for r in complete if r["start_flag"] == "NON-CANONICAL"]
            nc_stop  = [r for r in complete if r["stop_flag"]  == "NON-CANONICAL"]

            w("**Start codons**")
            w("")
            w("| Codon | n | % |")
            w("|:------|--:|--:|")
            for codon, c in sorted(starts.items(), key=lambda x: (-x[1], x[0])):
                w(f"| `{codon}` | {c} | {100*c/nc:.1f}% |")
            w("")
            if nc_start:
                w(f"**Non-canonical start: {len(nc_start)}/{nc} ({100*len(nc_start)/nc:.1f}%)**")
            else:
                w(f"Non-canonical start: {len(nc_start)}/{nc} ({100*len(nc_start)/nc:.1f}%)")
            w("")

            w("**Stop codons**")
            w("")
            w("| Codon | n | % |")
            w("|:------|--:|--:|")
            n_polya = sum(1 for r in complete if r["stop_flag"] == "polyA")
            for codon, c in sorted(stops.items(), key=lambda x: (-x[1], x[0])):
                label = f"`{codon}` *(polyA)*" if codon in ("T", "TA") else f"`{codon}`"
                w(f"| {label} | {c} | {100*c/nc:.1f}% |")
            w("")
            if n_polya:
                w(f"*polyA-completed (T/TA → TAA): {n_polya}/{nc} ({100*n_polya/nc:.1f}%)*")
                w("")
            if nc_stop:
                w(f"**Non-canonical stop: {len(nc_stop)}/{nc} ({100*len(nc_stop)/nc:.1f}%)**")
            else:
                w(f"Non-canonical stop: {len(nc_stop)}/{nc} ({100*len(nc_stop)/nc:.1f}%)")
            w("")

        if partial:
            w(f"*{len(partial)} partial CDS (truncated/fuzzy ends — not assessed)*")
            w("")

        w("---")
        w("")

    return "\n".join(lines)


# ── Markdown + PDF report ────────────────────────────────────────────────────

def _md_table(rows):
    """Markdown table of all CDS rows, NON-CANONICAL rows sorted first."""
    PRIORITY = {"NON-CANONICAL": 3, "partial": 2, "unknown": 1, "ok": 0, "polyA": 0}
    cols = ["record_id", "phylum", "order", "family", "species", "gene", "transl_table",
            "start_codon", "start_flag", "stop_codon", "stop_flag", "notes"]
    header = "| " + " | ".join(cols) + " |"
    sep    = "| " + " | ".join("---" for _ in cols) + " |"
    lines  = [header, sep]
    for r in sorted(rows, key=lambda r: -max(PRIORITY.get(r["start_flag"], 0),
                                              PRIORITY.get(r["stop_flag"], 0))):
        cells = []
        for c in cols:
            val = str(r[c]) if r[c] is not None else ""
            if c in ("start_flag", "stop_flag") and val == "NON-CANONICAL":
                val = f"**{val}**"
            cells.append(val)
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def write_markdown_report(rows, out_prefix, genbank_path, cmd_str=None):
    """Write <prefix>.report.md (with PNG figures) and <prefix>.report.pdf.

    Requires matplotlib (figures) and pandoc + weasyprint (PDF conversion).
    Markdown is always written; PDF is skipped gracefully if pandoc is absent.
    """
    import subprocess

    if not _HAVE_MPL:
        sys.stderr.write("WARNING: matplotlib not available; skipping report.\n")
        return None, None

    complete = [r for r in rows
                if r["start_flag"] != "partial" and r["stop_flag"] != "partial"]

    # save figure PNGs alongside the other output files
    figs = {}
    for tag, build in [
        ("fig_start_heatmap", lambda: _heatmap_fig(complete, "start_codon", "start_flag", "Start codon") if complete else None),
        ("fig_stop_heatmap",  lambda: _heatmap_fig(complete, "stop_codon",  "stop_flag",  "Stop codon")  if complete else None),
        ("fig_stacked_bars",  lambda: _stacked_bar_fig(complete) if complete else None),
        ("fig_geneLen_ridgeline", lambda: _length_ridgeline_fig(rows) if rows else None),
    ]:
        fig = build()
        if fig is not None:
            path = f"{out_prefix}.{tag}.png"
            _fig_to_file(fig, path)
            figs[tag] = os.path.abspath(path)

    n_records = len(set(r["record_id"] for r in rows))
    date_str  = datetime.date.today().isoformat()
    tsv_path  = f"{out_prefix}.codons.tsv"

    header = [
        "# Mitogenome Codon Review",
        "",
        f"**Generated:** {date_str}  |  "
        f"**Records:** {n_records}  |  "
        f"**CDS analyzed:** {len(rows)}",
        "",
    ]
    if cmd_str:
        header += [f"**Command:** `{cmd_str}`", ""]
    header += [
        f"**Input:** `{genbank_path}`",
        "",
        f"**Full per-CDS data table:** `{tsv_path}`",
        "",
        "---",
        "",
    ]

    md = header + [_build_summary_md(rows)]

    section_map = [
        ("fig_geneLen_ridgeline", "## 1. CDS Length Distribution",
         "CDS length distribution"),
        ("fig_start_heatmap", "## 2. Start Codon Frequency Heatmap",
         "Start codon heatmap"),
        ("fig_stop_heatmap",  "## 3. Stop Codon Frequency Heatmap",
         "Stop codon heatmap"),
        ("fig_stacked_bars",  "## 4. Start and Stop Codon Distributions",
         "Codon distributions"),
    ]
    for tag, heading, alt in section_map:
        if tag in figs:
            md += [heading, "", f"![{alt}]({figs[tag]})", ""]

    md_path  = f"{out_prefix}.report.md"
    pdf_path = f"{out_prefix}.report.pdf"
    with open(md_path, "w") as fh:
        fh.write("\n".join(md))

    cmd = ["pandoc", md_path, "-o", pdf_path, "--pdf-engine=weasyprint"]
    try:
        subprocess.run(cmd, check=True, capture_output=True)
    except FileNotFoundError:
        sys.stderr.write(
            f"WARNING: pandoc not found; PDF not generated. "
            f"Markdown written to {md_path}\n")
        pdf_path = None
    except subprocess.CalledProcessError as e:
        sys.stderr.write(
            f"WARNING: pandoc failed (exit {e.returncode}); PDF not generated.\n"
            f"  stderr: {e.stderr.decode().strip()}\n"
            f"  Markdown written to {md_path}\n")
        pdf_path = None

    return md_path, pdf_path


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Summarize start and stop codon usage by gene across all mitogenomes\n"
            "in a GenBank file, flagging codons that are non-canonical for the\n"
            "appropriate genetic code."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Details:
  Codons are read directly from CDS features, not from raw FASTA sequence, so
  the analysis correctly handles strand, spliced (join) locations, reading
  frame (/codon_start), and the per-record genetic code (/transl_table).
  Canonical codon sets come from Bio.Data.CodonTable, matching the NCBI table
  at https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes

Flags in output:
  ok             codon is in the canonical set for this record's genetic code
  NON-CANONICAL  codon is NOT in the canonical set  [FLAG]
  polyA          T or TA stop codon, completed to canonical TAA by polyadenylation
  partial        genuinely truncated/fuzzy end; not judged
  unknown        /transl_table absent or unrecognized; cannot judge

Output files:
  <prefix>.codons.tsv         one row per CDS (record_id, phylum, order,
                                family, species, gene, gene_raw,
                                transl_table, cds_length, codon_start,
                                start_codon, start_canonical, start_flag,
                                stop_codon, stop_canonical, stop_flag, notes)
  <prefix>.summary.txt        plain-text summary (same as stdout)
  <prefix>.report.md          Markdown report with embedded figure references
  <prefix>.report.pdf         PDF rendered from the Markdown (needs pandoc +
                                weasyprint; skipped gracefully if absent)
  <prefix>.fig_*.png          individual figure PNGs referenced by the report

Examples:
  python summarize_codons.py mitogenomes.gb
  python summarize_codons.py mitogenomes.gb -o results/percidae
""")
    ap.add_argument("genbank", help="input GenBank (.gb) file of mitogenomes")
    ap.add_argument("-o", "--out-prefix", default=None,
                    help="output file prefix (default: input filename without extension)")
    ap.add_argument("--api-key", default=None, metavar="KEY",
                    help="NCBI API key (raises taxonomy fetch rate from 3 to 10 req/s)")
    args = ap.parse_args()

    taxon_by_record = {}
    rows = []
    with open(args.genbank, "r") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            taxon_by_record[record.id] = _get_taxon_id(record)
            for feature in record.features:
                if feature.type != "CDS":
                    continue
                try:
                    rows.append(analyze_cds(feature, record))
                except Exception as e:
                    sys.stderr.write(
                        f"WARNING: skipping a CDS in {record.id}: {e}\n")

    if not rows:
        print("No CDS features found in the input. Nothing to summarize.")
        return

    # fetch higher taxonomy from NCBI for all unique taxon IDs
    unique_tids = {t for t in taxon_by_record.values() if t}
    print(f"Fetching taxonomy for {len(unique_tids)} taxon ID(s) from NCBI...")
    taxonomy = _fetch_taxonomy(unique_tids, api_key=args.api_key)
    for row in rows:
        info = taxonomy.get(str(taxon_by_record.get(row["record_id"]) or ""), {})
        row["phylum"]     = info.get("phylum",     "")
        row["order"]      = info.get("order",      "")
        row["family"]     = info.get("family",     "")

    out_prefix   = args.out_prefix or args.genbank.rsplit(".", 1)[0]
    tsv_path     = f"{out_prefix}.codons.tsv"
    summary_path = f"{out_prefix}.summary.txt"

    write_tsv(rows, tsv_path)
    print_summary(rows, summary_path)

    print("=" * 80)
    print(f"Per-CDS table:  {tsv_path}")
    print(f"Summary:        {summary_path}")

    argv = sys.argv[:]
    for i, tok in enumerate(argv):
        if tok == "--api-key" and i + 1 < len(argv):
            argv[i + 1] = "***"
    cmd_str = " ".join(argv)
    md_path, pdf_path = write_markdown_report(rows, out_prefix, args.genbank, cmd_str)
    if md_path:
        print(f"Markdown report: {md_path}")
    if pdf_path:
        print(f"PDF report:      {pdf_path}")
    print("=" * 80)


if __name__ == "__main__":
    main()
