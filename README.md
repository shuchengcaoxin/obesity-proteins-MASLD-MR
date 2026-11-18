# Proteome-wide Mendelian randomization prioritizes proteins linking obesity traits with MASLD

Analysis code accompanying the manuscript.

> **Citation:** Cao S, Su C-Y, Tan D, Sasako T, Yoshiji S, Butler-Laporte G.
> *Proteome-wide Mendelian randomization prioritizes proteins linking obesity traits with MASLD.*
> Manuscript in preparation.

---

## Overview

We applied a two-step Mendelian randomization (MR) framework to prioritize circulating
proteins that link adiposity to metabolic dysfunction-associated steatotic liver
disease (MASLD):

1. **Step 1 MR:** effect of three obesity proxies (BMI, waist-to-hip ratio, body fat
   percentage) on plasma protein levels, in two proteomic platforms
   (deCODE/SomaScan and UKB-PPP/Olink).
2. **Reverse MR:** protein &rarr; obesity, to exclude reverse causation.
3. **Step 2 MR:** effect of the resulting obesity-driven proteins on MASLD risk,
   using *cis*-pQTLs as instruments.
4. **Colocalization:** SharePro, to test whether the protein and MASLD signals share
   a causal variant.

## Repository layout

| Directory | Purpose |
|---|---|
| `1.1.BMI_to_proteins/` | Step 1 MR: obesity proxy &rarr; plasma proteome (instrument prep, MR, result summary) |
| `1.2.BMI_reverse/` | Reverse MR: plasma protein &rarr; obesity proxy, used to flag reverse causation |
| `1.3.collect_BMI_driven_proteins/` | Filters step 1 output (Bonferroni, heterogeneity, pleiotropy, reverse causation) to the set of obesity-driven proteins |
| `2.1.pQTL_to_MASLD_GWAS_MVP_new_filtered/` | Step 2 MR: obesity-driven proteins &rarr; MASLD (MVP), collection and filtering of results |
| `3.colocalization/` | Preparation of SharePro input files for the colocalization analyses |

Scripts are numbered in the order they were run. `.R` files are the analysis steps;
`.sh` files are the corresponding HPC job submission scripts.

## Data sources

All analyses use publicly available GWAS summary statistics; no individual-level data
were accessed, and no data are included in this repository.

| Data | Source |
|---|---|
| BMI | GIANT consortium + UK Biobank (Yengo et al. 2018) |
| Waist-to-hip ratio | GIANT consortium + UK Biobank (Pulit et al. 2019) |
| Body fat percentage | MRC-IEU UK Biobank GWAS, OpenGWAS accession `ukb-b-8909` |
| Plasma proteome (deCODE) | https://www.decode.com/summarydata/ |
| Plasma proteome (UKB-PPP) | Synapse `syn51364943` (doi:10.7303/syn51364943) |
| MASLD | Million Veteran Program |
| eQTL | GTEx v8 |

The only individual-level genotypes used are from the open-access 1000 Genomes Project
European reference panel, for LD clumping and LD estimation.

## Notes on reproducing

- Paths in the scripts point to the HPC directories on which the analyses were run and
  will need to be adapted.
- *cis*-pQTL instruments were re-derived from the full summary statistics:
  - deCODE: genome-wide significant variants (*P* < 5&times;10<sup>-8</sup>) within
    &plusmn;1 Mb of the encoding gene's transcription start site, LD-clumped at
    *r*<sup>2</sup> < 0.001 within a 10-Mb window.
  - UKB-PPP: LD clumping per protein GWAS (*P* < 5&times;10<sup>-8</sup>,
    *r*<sup>2</sup> < 0.001, 1-Mb window), then variants within &plusmn;500 kb of the
    canonical transcription start site of the protein-coding gene.
- Obesity instruments, *cis*-pQTL clumping and the colocalization LD matrices all used
  the 1000 Genomes Project European reference panel.

## Software

`R` with `TwoSampleMR`, `ieugwasr`, `vroom`, `tidyverse`, `magrittr`;
`PLINK 1.9`; `SharePro`; `locuszoomr`.
