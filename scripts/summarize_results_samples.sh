#!/bin/bash

run="minimum_reads_fix"
echo -e "sample\thp\thaplotype_name\t#YES\tmedian_of_mismatches\tmismatches_per_read" > /ifs/data/research/projects/michelle/17_prdm9_haplotypecheck/results/${run}/PRDM9_haplotypes.tsv
cat /ifs/data/research/projects/michelle/17_prdm9_haplotypecheck/results/${run}/*/*haplotypematch.tsv >> /ifs/data/research/projects/michelle/17_prdm9_haplotypecheck/results/${run}/PRDM9_haplotypes.tsv