# GenBankMitoReview

Download mitogenomes from GenBank for a group of organisms and review their
start/stop codon usage by gene, flagging codons that are non-canonical for the
appropriate genetic code.

## What it does

1. **Download** all mitochondrial genomic records matching a GenBank query into
   a `.gb` file, dropping records with one or fewer genes so only multi-gene
   mitogenomes remain.
2. **Summarize** start and stop codon usage per gene. Each
   record's `/transl_table` qualifier selects the
   [NCBI genetic code](https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes)
   used to judge whether a codon is canonical. Non-canonical start/stop codons
   are flagged.
3. **Calculate** and summarize protein-coding gene lengths.
4. **Summarize tRNA, rRNA, and control region** counts and lengths per genome,
   with dedicated plots — automatically enabled when `--complete` or `--refseq`
   is used.

Codons are read from the CDS features themselves, respecting strand, spliced
(`join`) locations, reading frame (`/codon_start`), and the genetic code.

## Getting started

```bash
git clone https://github.com/dmacguigan/GenBankMitoReview.git
cd GenBankMitoReview
```

## Setup

Create the `mitoreview` environment with mamba or conda (mamba is faster):

```bash
mamba env create -f environment.yml   # preferred
# or
conda env create -f environment.yml
```

This installs Entrez Direct (`esearch`/`efilter`/`efetch`), Biopython, matplotlib,
and the report-generation tools at pinned versions. No manual activation needed —
`mitoreview.sh` detects mamba or conda automatically and re-invokes itself inside
the environment.

## Usage

```bash
bash mitoreview.sh [--complete] [--refseq] '<query>' <output_prefix>
```

- `<query>`: an advanced GenBank search term in **single quotes**, e.g.
  `'"Percidae"[Organism]'` or `'"Percidae"[Organism] AND "PRJNA720393"[BioProject]'`
  (see <https://www.ncbi.nlm.nih.gov/nuccore/advanced>).
- `--complete`: restrict to records titled "complete genome"; also enables
  tRNA/rRNA/control region summary stats and plots.
- `--refseq`: restrict to RefSeq records (`NC_*` accessions); also enables
  tRNA/rRNA/control region summary stats and plots.
- `--outdir <dir>`: write all output to this directory (created if absent).
- `--api-key <key>`: NCBI API key for higher request rate.

Example:

```bash
bash mitoreview.sh '"Percidae"[Organism]' percidae
```

The two steps can also be run separately:

```bash
bash download_mitogenomes.sh '"Percidae"[Organism]' percidae   # -> percidae.gb
python summarize_codons.py percidae.gb -o percidae             # -> per-CDS table + report
```

## Output

- `<prefix>.gb` — the downloaded multi-gene mitogenomes.
- `<prefix>.raw.gb` — unfiltered download (kept for reference).
- `<prefix>.codons.tsv` — one row per CDS: gene (normalized), genetic-code
  table, codon-start, start/stop codon, canonical status, flag, and notes.
- `<prefix>.summary.txt` — per-gene codon summary (same text as stdout).
- `<prefix>.report.md` + `<prefix>.report.pdf` — Markdown/PDF report with
  figures and per-gene length statistics.
- `<prefix>.fig_start_heatmap.png` — start codon frequency heatmap (genes × codons).
- `<prefix>.fig_stop_heatmap.png` — stop codon frequency heatmap (genes × codons).
- `<prefix>.fig_stacked_bars.png` — stacked bar chart of start and stop codon distributions.
- `<prefix>.fig_geneLen_ridgeline.png` — ridgeline plot of CDS length distributions by gene.
- A per-gene text summary printed to screen: codon distributions, length stats,
  and the count and identity of non-canonical calls.

With `--complete` or `--refseq`, the pipeline also produces:

- `<prefix>.noncoding.tsv` — one row per tRNA, rRNA, or control region feature:
  `record_id`, `feature_type`, `gene_name`, `length`.
- `<prefix>.fig_noncoding_counts.png` — box + strip plot of tRNA, rRNA, and
  control region counts per genome.
- `<prefix>.fig_trna_lengths.png` — ridgeline plot of tRNA length distributions
  by gene (sorted by median length).
- `<prefix>.fig_rrna_cr_lengths.png` — box plot of rRNA (12S/16S) lengths and
  histogram of control region lengths.

A `polyA` flag marks T or TA stop codons completed to canonical TAA by
polyadenylation; these are counted as complete and canonical in all summaries and
figures. A `partial` flag marks genuinely truncated or fuzzy CDS ends.

Gene names are normalized to canonical metazoan PCG labels from the `/gene` and
`/product` qualifiers, so synonyms group together. The normalizer is ported from
the MitoPilot package's
[`normalize_pcg`](https://github.com/Smithsonian/MitoPilot/blob/main/R/blast_ref_utils.R).

Recognized genes (bilaterian standard 13 + non-bilaterian extras):

| Symbol | Gene | Taxon notes | Reference |
|--------|------|-------------|-----------|
| `nad1`-`nad6`, `nad4l` | NADH dehydrogenase subunits | universal metazoan mito | — |
| `cox1`-`cox3` | cytochrome c oxidase subunits | universal metazoan mito | — |
| `atp6`, `atp8` | ATP synthase subunits 6/8 | universal metazoan mito | — |
| `cob` | cytochrome b | universal metazoan mito | — |
| `mtMutS` | MutS mismatch-repair | octocorals; MutS/msh1/mtmsh synonyms | [Pont-Kingdon et al. 1995](https://doi.org/10.1038/375109b0) |
| `atp9` | ATP synthase subunit 9 | some sponges, cnidarians | [Cardona et al. 2014](https://doi.org/10.1016/j.ympev.2013.07.016) |
| `tatC` | twin-arginine translocase subunit C | some sponges (Oscarellidae) | [Gazave et al. 2010](https://doi.org/10.1371/journal.pone.0014290) |
| `rvt` | reverse transcriptase | bryozoans, annelids, sponges; group II intron-encoded | [Jenkins et al. 2022](https://doi.org/10.1038/s41598-022-14477-3) |
| `im` | intron maturase | bryozoans (cheilostomes); group II intron-encoded | [Jenkins et al. 2022](https://doi.org/10.1038/s41598-022-14477-3) |

Every CDS is analyzed regardless of name; any protein-coding gene beyond the
recognized set is still summarized, keeping its raw name prefixed with `?` so
non-standard or novel genes stay visible rather than being dropped.

## Files

- `environment.yml` — pinned mamba/conda environment.
- `mitoreview.sh` — main entry point; handles env activation + runs both steps.
- `download_mitogenomes.sh` — GenBank download + gene-count filter.
- `filter_gb.py` — drops records with one or fewer genes.
- `summarize_codons.py` — per-gene start/stop codon analysis, flagging, and reports.
- `bin/curl` — curl shim that forces HTTP/1.1 (fixes an OpenSSL 3.x / edirect issue).
