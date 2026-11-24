# Genome-Wide Codon Usage and Translational Efficiency Analysis

A Python-based computational biology project investigating genome-wide codon usage patterns in *Drosophila melanogaster* and their potential relationship to translational efficiency.

## Project Overview

This project develops a reproducible Python pipeline for processing and analyzing annotated coding DNA sequences (CDS) from *Drosophila melanogaster*.

The initial analysis quantifies codon usage across the genome, generating structured datasets of codon counts and relative frequencies. Building on this foundation, the research investigates how codon usage patterns vary across amino acids and nucleotide positions, and how these patterns may relate to translational efficiency.

A central research objective is to investigate the relationship between codon composition, ribosomal occupancy, RNA expression, and protein abundance.

## Research Objectives

The project focuses on several interconnected research questions:

- **Genome-wide codon usage:** Which codons occur most frequently across annotated coding sequences?
- **Codon preferences:** What is the most frequently occurring codon for each amino acid, and how do codon preferences vary?
- **Nucleotide composition:** How do nucleotide patterns and positional distributions contribute to observed codon usage?
- **Translational efficiency:** Are codon usage patterns associated with differences in ribosomal occupancy, RNA expression, and protein abundance?
- **Positional patterns:** Can nucleotide and codon distributions be organized into positional bins to identify changes in sequence composition?

## Dataset

| Attribute | Description |
|---|---|
| Organism | *Drosophila melanogaster* |
| Genome assembly | GCA_000001215.4 |
| Data type | Coding DNA sequences (CDS) |
| Input format | FASTA |
| Records processed | 30,802 |
| Codon types | 64 canonical codons |

The initial dataset consists of CDS sequences obtained from the NCBI genome assembly.

Additional expression and ribosomal occupancy datasets will support downstream investigations of translational efficiency.

## Computational Methods

### 1. Sequence Processing and Quality Control

The Python pipeline processes genomic CDS records from FASTA format.

Current functionality includes:
- Parsing FASTA records and extracting sequence identifiers.
- Extracting gene annotations where available.
- Processing nucleotide sequences into codon triplets.
- Handling incomplete sequence triplets and excluding ambiguous codons from canonical codon counts.
- Tracking record-level and dataset-level quality-control metrics.

### 2. Codon Frequency Analysis

Codon usage is quantified across individual CDS records and aggregated at the genome-wide level.

Analyses include:
- Counting occurrences of the 64 canonical codons.
- Calculating genome-wide codon frequencies and relative frequencies.
- Annotating codons with their corresponding amino acids.
- Identifying the most frequently occurring codon for each amino acid.

### 3. Nucleotide and Positional Analysis

Further analysis investigates nucleotide composition and positional patterns within coding sequences.

Planned investigations include:
- Grouping codons and nucleotide positions into defined positional bins.
- Examining changes in nucleotide composition across positions.
- Identifying patterns in codon distributions and nucleotide frequencies.
- Visualizing positional changes in codon usage.

### 4. Translational Efficiency Analysis

A longer-term research objective is to investigate potential relationships between codon usage and translational efficiency.

The analysis will incorporate complementary biological measurements, including:

- **Ribosomal occupancy:** Ribosome association with RNA relative to RNA expression.
- **Protein-to-RNA ratios:** Relationships between protein abundance and corresponding RNA expression.
- **Codon usage patterns:** Associations between codon frequencies, nucleotide composition, and translational measurements.

These datasets will be integrated with the codon usage analysis to investigate potential relationships between sequence composition and protein production.

## Planned Visualizations

The project will develop visualizations to communicate genome-wide and positional patterns.

Planned figures include:

- **Codon frequency heatmap:** A 20-amino-acid by codon visualization showing codon usage patterns and the most frequently occurring codon for each amino acid.
- **Nucleotide composition plots:** Visualizations of nucleotide frequencies across positional bins.
- **Positional distribution plots:** Comparisons of codon and nucleotide patterns across defined sequence regions.
- **Translational efficiency comparisons:** Visualizations investigating relationships between codon usage, ribosomal occupancy, and protein-to-RNA ratios.

## Current Outputs

The initial computational workflow generates:

- Per-record codon usage counts in CSV format.
- Genome-wide codon counts and relative frequencies.
- An Excel analysis workbook with aggregated codon usage data.
- Quality-control summaries.
- Amino-acid annotations for the 64 canonical codons.

## Technology Stack

- **Programming:** Python
- **Data Analysis:** pandas
- **Data Processing:** FASTA sequence parsing and nucleotide sequence analysis
- **Development:** Visual Studio Code
- **Version Control:** Git and GitHub
- **Data Exploration:** Microsoft Excel

## Project Status

The initial genome-wide codon counting and frequency analysis pipeline has been developed. Current research is expanding this foundation toward positional nucleotide analysis, codon preference visualization, and integration of translational measurements.

Future development will focus on extending the computational workflow, generating research visualizations, and investigating the relationship between codon usage and translational efficiency.

---

**Research Project:** BI 401 Undergraduate Research  
**Institution:** Oregon State University  
**Organism:** *Drosophila melanogaster*
