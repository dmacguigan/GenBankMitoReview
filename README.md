# GenBankMitoReview

Download mitogenomes from GenBank for a group of organisms and review their
start/stop codon usage by gene, flagging codons that are non-canonical for the
appropriate genetic code.

## What it does

1. **Download** all mitochondrial genomic records matching a GenBank query into
   a `.gb` file, dropping records with one or fewer genes so only multi-gene
   mitogenomes remain.
2. **Summarize** start and stop codons per gene across those mitogenomes. Each
   record's `/transl_table` qualifier selects the
   [NCBI genetic code](https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes)
   used to judge whether a codon is canonical. Non-canonical start/stop codons
   are flagged.

Codons are read from the CDS features themselves, respecting strand, spliced
(`join`) locations, reading frame (`/codon_start`), and the genetic code, so the
calls are accurate rather than a literal first/last 3 bp of a sequence.

## Setup

```bash
mamba env create -f environment.yml
```

This installs Entrez Direct (`esearch`/`efilter`/`efetch`), Biopython, matplotlib,
and the report-generation tools at pinned versions. No manual activation needed —
`mitoreview.sh` handles it automatically.

## Usage

End to end:

```bash
bash mitoreview.sh [--complete] [--refseq] '<query>' <output_prefix>
```

- `<query>`: an advanced GenBank search term in **single quotes**, e.g.
  `'"Percidae"[Organism]'` or `'"Percidae"[Organism] AND "PRJNA720393"[BioProject]'`
  (see <https://www.ncbi.nlm.nih.gov/nuccore/advanced>).
- `--complete`: restrict to records titled "complete genome".
- `--refseq`: restrict to RefSeq records (`NC_*` accessions).
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
  heatmaps, stacked bar charts, and the full CDS data table.
- A per-gene text summary printed to screen: codon distributions and the count
  and identity of non-canonical calls.

A `partial` flag (rather than non-canonical) marks truncated/fuzzy CDS ends,
such as mitochondrial stop codons completed by polyadenylation.

Gene names are normalized to canonical metazoan PCG labels from the `/gene` and
`/product` qualifiers, so synonyms group together. The normalizer is ported from
the MitoPilot package's
[`normalize_pcg`](https://github.com/Smithsonian/MitoPilot/blob/main/R/blast_ref_utils.R).

Recognized genes (bilaterian standard 13 + non-bilaterian extras):

| Symbol | Gene | Taxon notes |
|--------|------|-------------|
| `nad1`-`nad6`, `nad4l` | NADH dehydrogenase subunits | universal metazoan mito |
| `cox1`-`cox3` | cytochrome c oxidase subunits | universal metazoan mito |
| `atp6`, `atp8` | ATP synthase subunits 6/8 | universal metazoan mito |
| `cob` | cytochrome b | universal metazoan mito |
| `mtMutS` | MutS mismatch-repair | octocorals; MutS/msh1/mtmsh synonyms |
| `nad7` | NADH dehydrogenase subunit 7 | some cnidarians (e.g., Hydra) |
| `atp9` | ATP synthase subunit 9 | Trichoplax, some sponges/cnidarians |
| `sdh2`, `sdh3`, `sdh4` | succinate dehydrogenase subunits B/C/D | some cnidarians, sponges |
| `rps3` | ribosomal protein S3 | cnidarians, sponges, some flatworms |
| `rpl16` | ribosomal protein L16 | some cnidarians, sponges |
| `tatC` | twin-arginine translocase subunit C | sponges (Amphimedon, Oscarella) |
| `polB` | DNA polymerase B | some demosponge mitogenomes |

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
