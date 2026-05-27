<style>
body {
    font-family: sans-serif;
    font-size: 11pt;
    line-height: 1.5;
    max-width: 960px;
    margin: 0 auto;
    padding: 1em 2em;
    color: #1a1a1a;
}
h1 {
    border-bottom: 3px solid #1a5276;
    padding-bottom: 6px;
    color: #1a5276;
}
h2 {
    border-bottom: 2px solid #aab7b8;
    padding-bottom: 4px;
    color: #2c3e50;
    margin-top: 2em;
}
h3 {
    border-left: 5px solid #2980b9;
    background-color: #eaf4fb;
    padding: 5px 10px;
    margin-top: 1.8em;
    color: #1a5276;
    font-family: monospace;
    font-size: 1.5em;
}
table {
    border-collapse: collapse;
    width: 100%;
    margin-bottom: 1em;
    font-size: 10pt;
}
th {
    background-color: #2c3e50;
    color: #ffffff;
    padding: 5px 10px;
    text-align: left;
}
td {
    padding: 4px 10px;
    border: 1px solid #d5d8dc;
}
tr:nth-child(even) td {
    background-color: #f4f6f7;
}
td strong {
    color: #c0392b;
    background-color: #fdecea;
    padding: 1px 4px;
    border-radius: 3px;
}
hr {
    border: none;
    border-top: 1px solid #d5d8dc;
    margin: 2em 0;
}
code {
    background-color: #f0f0f0;
    padding: 1px 4px;
    border-radius: 3px;
    font-size: 0.95em;
}
</style>
# Mitogenome Codon Review

**Generated:** 2026-05-27  |  **Records:** 30  |  **CDS analyzed:** 388

**Command:** `/home/dmacguig/Documents/GitHub/GenBankMitoReview/summarize_codons.py examples/asteroidea.gb -o examples/asteroidea --noncoding-stats`

**Input:** `examples/asteroidea.gb`

**Full per-CDS data table:** `examples/asteroidea.codons.tsv`

---

## Genetic Code Summary

**Table 9: Echinoderm Mitochondrial** — 30 record(s)
- Canonical start codons: `ATG` · `GTG`
- Canonical stop codons: `TAA` · `TAG`

---

## Codon Usage by Gene

*388 CDS analyzed across 13 gene group(s). Canonical status uses each record's `/transl_table` genetic code.*

### `atp6`

**30 CDS** · 30 complete · 0 partial

**Length (bp):** min=687 · max=693 · mean=692.3 · SD=1.8

**Start codons**

| Codon | n | % | Status |
|:------|--:|--:|:-------|
| `ATG` | 30 | 100.0% | ok |

Non-canonical start: 0/30 (0.0%)

**Stop codons**

| Codon | n | % | Status |
|:------|--:|--:|:-------|
| `TAA` | 28 | 93.3% | ok |
| `T` *(polyA)* | 1 | 3.3% | polyA |
| `TAG` | 1 | 3.3% | ok |

*polyA-completed (T/TA → TAA): 1/30 (3.3%)*

Non-canonical stop: 0/30 (0.0%)

---

### `atp8`

**29 CDS** · 29 complete · 0 partial

**Length (bp):** min=162 · max=177 · mean=167.2 · SD=3.2

**Start codons**

| Codon | n | % | Status |
|:------|--:|--:|:-------|
| `ATG` | 29 | 100.0% | ok |

Non-canonical start: 0/29 (0.0%)

**Stop codons**

| Codon | n | % | Status |
|:------|--:|--:|:-------|
| `TAA` | 28 | 96.6% | ok |
| `TAG` | 1 | 3.4% | ok |

Non-canonical stop: 0/29 (0.0%)

---

### `cob`

**30 CDS** · 30 complete · 0 partial

**Length (bp):** min=1138 · max=1158 · mean=1138.7 · SD=3.6

**Start codons**

| Codon | n | % | Status |
|:------|--:|--:|:-------|
| `ATG` | 30 | 100.0% | ok |

Non-canonical start: 0/30 (0.0%)

**Stop codons**

| Codon | n | % | Status |
|:------|--:|--:|:-------|
| `T` *(polyA)* | 29 | 96.7% | polyA |
| `TAA` | 1 | 3.3% | ok |

*polyA-completed (T/TA → TAA): 29/30 (96.7%)*

Non-canonical stop: 0/30 (0.0%)

---

### `cox1`

**30 CDS** · 30 complete · 0 partial

**Length (bp):** min=1551 · max=1557 · mean=1553.8 · SD=0.9

**Start codons**

| Codon | n | % | Status |
|:------|--:|--:|:-------|
| `ATG` | 30 | 100.0% | ok |

Non-canonical start: 0/30 (0.0%)

**Stop codons**

| Codon | n | % | Status |
|:------|--:|--:|:-------|
| `TAA` | 21 | 70.0% | ok |
| `TA` *(polyA)* | 5 | 16.7% | polyA |
| `TAG` | 3 | 10.0% | ok |
| `T` *(polyA)* | 1 | 3.3% | polyA |

*polyA-completed (T/TA → TAA): 6/30 (20.0%)*

Non-canonical stop: 0/30 (0.0%)

---

### `cox2`

**30 CDS** · 30 complete · 0 partial

**Length (bp):** min=688 · max=702 · mean=689.2 · SD=2.7

**Start codons**

| Codon | n | % | Status |
|:------|--:|--:|:-------|
| `ATG` | 30 | 100.0% | ok |

Non-canonical start: 0/30 (0.0%)

**Stop codons**

| Codon | n | % | Status |
|:------|--:|--:|:-------|
| `T` *(polyA)* | 19 | 63.3% | polyA |
| `TAA` | 9 | 30.0% | ok |
| `TAG` | 2 | 6.7% | ok |

*polyA-completed (T/TA → TAA): 19/30 (63.3%)*

Non-canonical stop: 0/30 (0.0%)

---

### `cox3`

**30 CDS** · 30 complete · 0 partial

**Length (bp):** min=780 · max=783 · mean=782.8 · SD=0.7

**Start codons**

| Codon | n | % | Status |
|:------|--:|--:|:-------|
| `ATG` | 30 | 100.0% | ok |

Non-canonical start: 0/30 (0.0%)

**Stop codons**

| Codon | n | % | Status |
|:------|--:|--:|:-------|
| `TAA` | 30 | 100.0% | ok |

Non-canonical stop: 0/30 (0.0%)

---

### `nad1`

**30 CDS** · 28 complete · 2 partial

**Length (bp):** min=948 · max=981 · mean=977.6 · SD=5.9

**Start codons**

| Codon | n | % | Status |
|:------|--:|--:|:-------|
| `GTG` | 17 | 60.7% | ok |
| `ATG` | 11 | 39.3% | ok |

Non-canonical start: 0/28 (0.0%)

**Stop codons**

| Codon | n | % | Status |
|:------|--:|--:|:-------|
| `TAG` | 15 | 53.6% | ok |
| `TAA` | 11 | 39.3% | ok |
| `T` *(polyA)* | 2 | 7.1% | polyA |

*polyA-completed (T/TA → TAA): 2/28 (7.1%)*

Non-canonical stop: 0/28 (0.0%)

*2 partial CDS (truncated/fuzzy ends — not assessed)*

| Record | Species | Start codon | Start flag | Stop codon | Stop flag |
|:-------|:--------|:------------|:-----------|:-----------|:----------|
| `NC_054231.1` | Euretaster insignis | `TTG` | partial | `TAG` | ok |
| `NC_085630.1` | Xyloplax princealberti | `TTG` | partial | `TAG` | ok |

---

### `nad2`

**30 CDS** · 30 complete · 0 partial

**Length (bp):** min=1062 · max=1095 · mean=1065.9 · SD=6.5

**Start codons**

| Codon | n | % | Status |
|:------|--:|--:|:-------|
| `ATG` | 15 | 50.0% | ok |
| `GTG` | 15 | 50.0% | ok |

Non-canonical start: 0/30 (0.0%)

**Stop codons**

| Codon | n | % | Status |
|:------|--:|--:|:-------|
| `TAA` | 23 | 76.7% | ok |
| `TAG` | 7 | 23.3% | ok |

Non-canonical stop: 0/30 (0.0%)

---

### `nad3`

**30 CDS** · 12 complete · 18 partial

**Length (bp):** min=333 · max=351 · mean=350.3 · SD=3.3

**Start codons**

| Codon | n | % | Status |
|:------|--:|--:|:-------|
| `ATG` | 8 | 66.7% | ok |
| `ATT` | 4 | 33.3% | **NON-CANONICAL** |

**Non-canonical start: 4/12 (33.3%)**

**Stop codons**

| Codon | n | % | Status |
|:------|--:|--:|:-------|
| `TAA` | 10 | 83.3% | ok |
| `TAG` | 2 | 16.7% | ok |

Non-canonical stop: 0/12 (0.0%)

*18 partial CDS (truncated/fuzzy ends — not assessed)*

| Record | Species | Start codon | Start flag | Stop codon | Stop flag |
|:-------|:--------|:------------|:-----------|:-----------|:----------|
| `NC_001627.1` | Patiria pectinifera | `ATT` | partial | `TAG` | ok |
| `NC_037943.1` | Echinaster brasiliensis | `ATT` | partial | `TAG` | ok |
| `NC_053361.1` | Crossaster papposus | `ATT` | partial | `TAG` | ok |
| `NC_054225.1` | Pentaceraster mammillatus | `ATT` | partial | `TAA` | ok |
| `NC_054226.1` | Protoreaster nodosus | `ATT` | partial | `TAA` | ok |
| `NC_054227.1` | Ophidiaster granifer | `ATT` | partial | `TAA` | ok |
| `NC_054228.1` | Iconaster longimanus | `ATT` | partial | `TAA` | ok |
| `NC_054229.1` | Culcita novaeguineae | `ATT` | partial | `TAA` | ok |
| `NC_054230.1` | Anthenea aspera | `ATT` | partial | `TAA` | ok |
| `NC_054231.1` | Euretaster insignis | `TTG` | partial | `TAA` | ok |
| `NC_063668.1` | Asthenactis papyraceus | `ATT` | partial | `TAA` | ok |
| `NC_063669.1` | Zoroaster ophiactis | `ATT` | partial | `TAA` | ok |
| `NC_081967.1` | Leptychaster arcticus | `ATT` | partial | `TAA` | ok |
| `NC_081981.1` | Crossaster japonicus | `ATT` | partial | `TAG` | ok |
| `NC_083191.1` | Poraniopsis inflata | `ATT` | partial | `TAA` | ok |
| `NC_085630.1` | Xyloplax princealberti | `ATT` | partial | `TAA` | ok |
| `NC_087784.1` | Henricia reniossa | `ATT` | partial | `TAA` | ok |
| `NC_087785.1` | Henricia sanguinolenta | `ATT` | partial | `TAA` | ok |

---

### `nad4`

**30 CDS** · 30 complete · 0 partial

**Length (bp):** min=1377 · max=1386 · mean=1382.5 · SD=1.6

**Start codons**

| Codon | n | % | Status |
|:------|--:|--:|:-------|
| `ATG` | 30 | 100.0% | ok |

Non-canonical start: 0/30 (0.0%)

**Stop codons**

| Codon | n | % | Status |
|:------|--:|--:|:-------|
| `TAA` | 23 | 76.7% | ok |
| `TAG` | 7 | 23.3% | ok |

Non-canonical stop: 0/30 (0.0%)

---

### `nad4l`

**30 CDS** · 8 complete · 22 partial

**Length (bp):** min=288 · max=321 · mean=297.3 · SD=4.7

**Start codons**

| Codon | n | % | Status |
|:------|--:|--:|:-------|
| `ATT` | 4 | 50.0% | **NON-CANONICAL** |
| `ATG` | 3 | 37.5% | ok |
| `ATC` | 1 | 12.5% | **NON-CANONICAL** |

**Non-canonical start: 5/8 (62.5%)**

**Stop codons**

| Codon | n | % | Status |
|:------|--:|--:|:-------|
| `TAA` | 8 | 100.0% | ok |

Non-canonical stop: 0/8 (0.0%)

*22 partial CDS (truncated/fuzzy ends — not assessed)*

| Record | Species | Start codon | Start flag | Stop codon | Stop flag |
|:-------|:--------|:------------|:-----------|:-----------|:----------|
| `NC_001627.1` | Patiria pectinifera | `ATT` | partial | `TAA` | ok |
| `NC_037943.1` | Echinaster brasiliensis | `ATT` | partial | `TAA` | ok |
| `NC_041450.1` | Styracaster yapensis | `ATT` | partial | `TAA` | ok |
| `NC_042741.1` | Pisaster ochraceus | `ATT` | partial | `TAA` | ok |
| `NC_053361.1` | Crossaster papposus | `GTC` | partial | `TAA` | ok |
| `NC_054225.1` | Pentaceraster mammillatus | `ATT` | partial | `TAA` | ok |
| `NC_054226.1` | Protoreaster nodosus | `ATT` | partial | `TAA` | ok |
| `NC_054227.1` | Ophidiaster granifer | `ATC` | partial | `TAA` | ok |
| `NC_054228.1` | Iconaster longimanus | `ATC` | partial | `TAA` | ok |
| `NC_054229.1` | Culcita novaeguineae | `ATC` | partial | `TAA` | ok |
| `NC_054230.1` | Anthenea aspera | `ATC` | partial | `TAA` | ok |
| `NC_054231.1` | Euretaster insignis | `ATG` | partial | `TAA` | ok |
| `NC_063669.1` | Zoroaster ophiactis | `ATT` | partial | `TAA` | ok |
| `NC_063787.1` | Coscinasterias acutispina | `ATT` | partial | `TAA` | ok |
| `NC_071876.1` | Ctenodiscus crispatus | `ATT` | partial | `TAA` | ok |
| `NC_081967.1` | Leptychaster arcticus | `ATT` | partial | `TAA` | ok |
| `NC_081981.1` | Crossaster japonicus | `GTT` | partial | `TAA` | ok |
| `NC_083190.1` | Sclerasterias satsumana | `ATT` | partial | `TAA` | ok |
| `NC_083191.1` | Poraniopsis inflata | `ATT` | partial | `TAA` | ok |
| `NC_085630.1` | Xyloplax princealberti | `TTG` | partial | `TAA` | ok |
| `NC_087784.1` | Henricia reniossa | `ATC` | partial | `TAA` | ok |
| `NC_087785.1` | Henricia sanguinolenta | `ATC` | partial | `TAA` | ok |

---

### `nad5`

**29 CDS** · 29 complete · 0 partial

**Length (bp):** min=1842 · max=1932 · mean=1903.0 · SD=18.3

**Start codons**

| Codon | n | % | Status |
|:------|--:|--:|:-------|
| `ATG` | 28 | 96.6% | ok |
| `GTG` | 1 | 3.4% | ok |

Non-canonical start: 0/29 (0.0%)

**Stop codons**

| Codon | n | % | Status |
|:------|--:|--:|:-------|
| `TAA` | 26 | 89.7% | ok |
| `TAG` | 3 | 10.3% | ok |

Non-canonical stop: 0/29 (0.0%)

---

### `nad6`

**30 CDS** · 30 complete · 0 partial

**Length (bp):** min=474 · max=492 · mean=488.9 · SD=3.0

**Start codons**

| Codon | n | % | Status |
|:------|--:|--:|:-------|
| `ATG` | 27 | 90.0% | ok |
| `GTG` | 3 | 10.0% | ok |

Non-canonical start: 0/30 (0.0%)

**Stop codons**

| Codon | n | % | Status |
|:------|--:|--:|:-------|
| `TAA` | 15 | 50.0% | ok |
| `TAG` | 15 | 50.0% | ok |

Non-canonical stop: 0/30 (0.0%)

---

## 1. CDS Length Distribution

![CDS length distribution](/home/dmacguig/Documents/GitHub/GenBankMitoReview/examples/asteroidea.fig_geneLen_ridgeline.png)

## 2. Start Codon Frequency Heatmap

![Start codon heatmap](/home/dmacguig/Documents/GitHub/GenBankMitoReview/examples/asteroidea.fig_start_heatmap.png)

## 3. Stop Codon Frequency Heatmap

![Stop codon heatmap](/home/dmacguig/Documents/GitHub/GenBankMitoReview/examples/asteroidea.fig_stop_heatmap.png)

## 4. Start and Stop Codon Distributions

![Codon distributions](/home/dmacguig/Documents/GitHub/GenBankMitoReview/examples/asteroidea.fig_stacked_bars.png)


---

## Non-coding Feature Summary

*707 non-coding features annotated across 30 record(s).*

### tRNA

**Count per genome** (n=30 records): min=19 · max=22 · mean=21.7

| Gene | n | min (bp) | max (bp) | mean (bp) | median (bp) |
|:-----|--:|--------:|--------:|----------:|------------:|
| tRNA-Ala | 30 | 68 | 73 | 71.1 | 72.0 |
| tRNA-Arg | 30 | 46 | 72 | 68.4 | 69.0 |
| tRNA-Asn | 30 | 72 | 73 | 72.0 | 72.0 |
| tRNA-Asp | 29 | 69 | 75 | 70.4 | 70.0 |
| tRNA-Cys | 30 | 66 | 73 | 69.3 | 69.0 |
| tRNA-Gln | 30 | 71 | 74 | 71.8 | 72.0 |
| tRNA-Glu | 30 | 67 | 73 | 69.2 | 69.0 |
| tRNA-Gly | 30 | 67 | 72 | 69.3 | 69.0 |
| tRNA-His | 30 | 67 | 72 | 69.2 | 69.0 |
| tRNA-Ile | 30 | 59 | 74 | 71.4 | 72.0 |
| tRNA-Leu | 60 | 71 | 74 | 72.4 | 72.0 |
| tRNA-Lys | 28 | 57 | 76 | 72.3 | 72.0 |
| tRNA-Met | 30 | 68 | 75 | 72.0 | 72.0 |
| tRNA-Phe | 30 | 50 | 76 | 71.1 | 72.0 |
| tRNA-Pro | 30 | 69 | 76 | 70.8 | 71.0 |
| tRNA-Ser | 57 | 67 | 73 | 69.8 | 70.0 |
| tRNA-Thr | 30 | 69 | 76 | 71.7 | 71.5 |
| tRNA-Trp | 27 | 69 | 73 | 69.6 | 69.0 |
| tRNA-Tyr | 30 | 69 | 72 | 70.7 | 71.0 |
| tRNA-Val | 30 | 62 | 75 | 70.6 | 71.0 |

**tRNA copy number per genome** (genes present >1x in at least one record)

| Gene | 1 copy | 2 copies | 3+ copies |
|:-----|-------:|---------:|----------:|
| tRNA-Leu | 0 | 30 | 0 |
| tRNA-Ser | 3 | 27 | 0 |

---

### rRNA

**Count per genome** (n=27 records): min=2 · max=2 · mean=2.0

| Gene | n | min (bp) | max (bp) | mean (bp) | median (bp) |
|:-----|--:|--------:|--------:|----------:|------------:|
| 12S rRNA | 27 | 885 | 929 | 902.1 | 898.0 |
| 16S rRNA | 27 | 1531 | 1646 | 1596.1 | 1603.0 |

---

### Control Region

**Count per genome** (n=2 records): min=1 · max=1 · mean=1.0

| Gene | n | min (bp) | max (bp) | mean (bp) | median (bp) |
|:-----|--:|--------:|--------:|----------:|------------:|
| control region | 2 | 531 | 552 | 541.5 | 541.5 |

---

## 5. Non-coding Feature Counts

![Non-coding feature counts](/home/dmacguig/Documents/GitHub/GenBankMitoReview/examples/asteroidea.fig_noncoding_counts.png)

## 6. tRNA Length Distribution

![tRNA lengths by gene](/home/dmacguig/Documents/GitHub/GenBankMitoReview/examples/asteroidea.fig_trna_lengths.png)

## 7. rRNA and Control Region Lengths

![rRNA and control region lengths](/home/dmacguig/Documents/GitHub/GenBankMitoReview/examples/asteroidea.fig_rrna_cr_lengths.png)
