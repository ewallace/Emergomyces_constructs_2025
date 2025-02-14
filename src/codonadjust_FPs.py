#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 12:06:49 2022

@author: Edward Wallace
"""
# This script uses DNAchisel to produces Histoplasma capsulatum codon-adjusted DNA sequences from an amino acid sqeuence
# Instruction for DNAchisel can be found at https://edinburgh-genome-foundry.github.io/DnaChisel/
# DNAchisel paper: https://doi.org/10.1093/bioinformatics/btaa558
# This script builds on anidulans_codon_adjust.py script from Domenico Modaffari

# import codon counting functions, EW adapted from CodonAdaptationIndex package
from codoncount_functions import *

# import relevant functions from DNAchisel
from dnachisel import (DnaOptimizationProblem,
                       EnforceTranslation,
                       CodonOptimize,
                       reverse_translate,
                       EnforceGCContent,
                       UniquifyAllKmers,
                       AvoidRareCodons,
                       AvoidPattern, #use this if you want to exclude restriction enzymes cut sites from DNA sequence
                       AvoidHairpins
                       )

# import SeqFeature, for labeling
from Bio import SeqFeature as sf

# import date, so we can time-stamp the output
from datetime import date

# Ensure reproducible output by setting the seed for the random number generator
# This is because the DnaOptimizationProblem uses some random processes in its optimization
from numpy import random


# Calculate the codon frequencies to use
HcCodonFreqs = CodonFrequencyByAA(
    fastaFileToStrings(
        "data/Histoplasma_capsulatum_ribosomalproteinCDS.fasta"
        ),
        outputtype = "freq"
    )

# Define constraints to avoid
# restriction enzyme sites, that rely on Bio.Restriction as called in dnachisel.AviodPattern
enzymes_to_avoid = ['HindIII', 'SwaI', 'AgeI', 'XhoI', 'SpeI', 'KpnI', 'NotI', 'EcoRI', 'BglII', 'SacI', 'SacII', 'BbsI', 'BsaI', 'BsmBI', 'AarI']
sites_to_avoid = [ AvoidPattern("%s_site" % enzyme)
    for enzyme in  enzymes_to_avoid ]

# also ensure translation, medium GC content, avoid rarest codons, avoid hairpins
other_constraints = [
    EnforceTranslation(),
    EnforceGCContent(mini=0.4, maxi=0.6, window=40), #GC content setting
    AvoidRareCodons(min_frequency = 0.1, codon_usage_table= HcCodonFreqs), #do not use codons below frequency of 0.1
    AvoidHairpins(stem_size=15) #avoid hairpins
    ] 

all_constraints = sites_to_avoid + other_constraints

# define objectives
objectives_match_uniquify = [
    CodonOptimize(method = "match_codon_usage", codon_usage_table = HcCodonFreqs), #adjust codons
    UniquifyAllKmers(k = 9) #try to avoid repeats in sequence, makes gene synthesys easier
    ]

# function that does all the adjustment/optimisation for a protein, and writes out to a record
def codon_adjust(proteinseq, 
        seed = None, 
        filepath = None, 
        record_id = None, 
        CDS_id = None,
        add_date = True,
        with_constraints = False,
        with_objectives = False):
    problem = DnaOptimizationProblem(
        sequence = reverse_translate(proteinseq),
        constraints = all_constraints,
        objectives = objectives_match_uniquify
        )
    
    if (seed is not None) :
        random.seed(seed)
    
    problem.resolve_constraints()
    problem.optimize()
    # convert to record
    record = problem.to_record(
        record_id = record_id,
        with_constraints = with_constraints,
        with_objectives = with_objectives)
    
    if (CDS_id is not None):
        # add the CDS feature into the record
        CDSfeature = sf.SeqFeature(
           location = sf.FeatureLocation(0, len(record.seq), strand = +1), 
           type = "CDS",
           id = CDS_id,
           qualifiers = {'gene':CDS_id})
        record.features = record.features + [ CDSfeature ] 
    
    if (add_date):
        # add date to record
        # note this doesn't write the date to the output file, don't know why
        record.annotations["date"] =  str(date.today())
    
    if (filepath is not None) :
        SeqIO.write(record, handle = filepath, format = "genbank")
    return record


# Apply to the proteins of interest. Write out to a record
codon_adjust(
    proteinseq = SeqIO.read("data/mCardinal.fasta", "fasta"),
    seed = 1,
    filepath = "designs/mCardinalHc_2025-02-14.gbk",
    record_id = "mCardinalHc",
    CDS_id = "mCardinal")


