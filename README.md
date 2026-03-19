# CRISPR_interactions

This is the collection of scripts used for the **"CRISPR-Cas immune repertoires as an ecological record of bacterial interactions with mobile genetic elements in the human gut"** by Avershina et al, 2026 [preprint](https://www.biorxiv.org/content/10.64898/2026.03.18.712547v1.full.pdf+html)

---

## Overview

This repository contains code and data processing workflows used for the extended microbiome repository generation and bacteria-MGE interaction analysis using CRISPR-Cas. 


---

## Repository Structure

```cctyper_snakemake_wf/ ```: Snakemake workflow for CRISPR-Cas data generation using [CRISPRCasTyper](https://github.com/Russel88/CRISPRCasTyper)

```scapp_snakemake_wf/ ```:  Snakemake workflow for plasmid data generation using [SCAPP](https://github.com/Shamir-Lab/SCAPP)

```data_analysis/ ```: Scripts used for processing of generated data

---

## Data Availability

* Raw sequencing data are available from The Federated European Genome-phenome Archive (FEGA): accession EGAS50000000170
* Prokaryotic metagenome-assembled genomes were described in [Birkeland et al, 2025](https://www.medrxiv.org/content/10.1101/2025.10.06.25336873v1) and generated using [ATLAS](https://github.com/metagenome-atlas/atlas)
* Viral data were described in [Istvan et al, 2026](https://www.nature.com/articles/s41467-024-46033-0) and generated using [VirMake](https://github.com/Rounge-lab/VirMake)
* Processed fastafiles and metadata are deposited in [Figshare](https://doi.org/10.6084/m9.figshare.31707688)

---

## Dependencies

* Python ≥ 3.9
* pandas
* numpy
* scikit-learn
* biopython
* stats
* matplotlib
* seaborn
* cctyper
* scapp

---

## Contact

Ekaterina Avershina<br> 
University of Oslo<br>
Email: ekateria@uio.no
<br>
<br>
Trine B Rounge<br> 
University of Oslo<br>
Email: t.b.rounge@uio.no
