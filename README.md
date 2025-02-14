# Emergomyces_constructs_2025

This code is to generate codon-adjusted fluorescent proteins for Emergomyces africanus and Histoplasma capsulatum.

Edward Wallace, Edward.Wallace@ed.ac.uk, 2025

Funded by the Wellcome Trust Bioimaging award.

# Contents

## data

Input data:

- native sequences for codon counts from FungiDB. Precisely, Histoplasma capsulatum or Emergomyces africanus CDSs with annotation GO:0003735 "structural constituent of ribosome". Downloaded by EW in February 2025.
- protein sequences for fluorescent proteins from fpbase.org. Downloaded by EW February 2025.

## designs

Coding sequences for fluorescent proteins, adjusted to codon frequencies of Histoplasma capsulatum ribosomal proteins, moderate GC content, and to avoid select restriction enzyme sites. See source code for details.

## src

Source code in python format. Uses the DNAchisel package and Biopython tools.

- codon_adjust_FPs.py - using DNAchisel's DnaOptimizationProblem to design sequences encoding fluorescent proteins
- codon_count_functions.py - functions to calculate codon counts and frequencies


# How to run the code

To runs the code to design the proteins, from the root directory of of `Emergomyces_constructs_2025`, in a bash terminal run:

```bash
python src/codon_adjust_FPs.py
```

## Computational environment

You'll need the relevant python packages installed. Biopython and DNAchisel, etc.

It should be possible to install the packages using conda:

```bash
conda create --name codons2025 --file codons2025_env.txt
conda activate codons2025
```

Beware, I have not tested this.
