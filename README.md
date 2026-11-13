# bi401-codon-usage
Python-based analysis of genome-wide codon usage in Drosophila melanogaster.
# Genome-Wide Codon Usage Analysis

## Project Overview
This project analyzes codon usage patterns in *Drosophila melanogaster* using Python. The goal is to characterize codon frequencies across annotated coding DNA sequences (CDS) and establish a reproducible computational workflow for downstream analysis of codon usage bias.

## Dataset
- **Organism:** *Drosophila melanogaster*
- **Genome assembly:** GCA_000001215.4
- **Data type:** Coding DNA sequences (CDS)
- **Records processed:** 30,802

## Methods
The Python pipeline:
- Parses CDS sequences from FASTA format.
- Extracts gene annotations where available.
- Counts codon occurrences for each CDS record.
- Generates a structured CSV containing per-record codon counts for the 64 canonical codons.

## Quality Control
Initial quality-control checks include:
- Number of CDS records processed.
- Unique record identifiers.
- Missing gene annotations.
- Handling of incomplete sequence triplets and canonical codon counts.

## Outputs
- Per-record codon usage counts (CSV).
- Genome-wide codon counts and relative frequencies (Excel analysis workbook).
- Quality-control summary.

## Future Directions
Planned winter-term work will investigate codon usage bias, GC and GC3 composition, variation across genes, and positional patterns, with the goal of developing a scientific report and research poster.

## Tools
- Python
- Visual Studio Code
- GitHub
- Microsoft Excel
