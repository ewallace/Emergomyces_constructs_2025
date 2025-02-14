# Emergomyces_constructs_2025

This code is to generate codon-adjusted fluorescent proteins for Emergomyces africanus and Histoplasma capsulatum.

Edward Wallace, Edward.Wallace@ed.ac.uk, 2025

Funded by the Wellcome Trust Bioimaging award.

# Contents

## data

Input data:

- native sequences for codon counts from FungiDB
- protein sequences for fluorescent proteins from fpbase.org

## designs

Coding sequences for fluorescent proteins, adjusted to codon frequencies of Histoplasma capsulatum ribosomal proteins, moderate GC content, and to avoid select restriction enzyme sites. See source code for details.

## src

Source code in python format. Uses the DNAchisel package and Biopython tools.

- codonadjust_FPs.py - using DNAchisel's DnaOptimizationProblem to make fluorescent proteins
- codoncount_functions.py - functions to calculate codon counts and frequencies
