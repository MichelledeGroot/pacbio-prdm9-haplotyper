#!/bin/bash

run="use_5_reads"
echo -e "sample\thp\thaplotype_name\t#YES\tmedian_of_mismatches\tmismatches_per_read\tmode\th1_reads\th2_reads" > /ifs/data/research/projects/michelle/17_prdm9_haplotypecheck/results/${run}/PRDM9_haplotypes.tsv
cat /ifs/data/research/projects/michelle/17_prdm9_haplotypecheck/results/${run}/*/*haplotypematch.tsv >> /ifs/data/research/projects/michelle/17_prdm9_haplotypecheck/results/${run}/PRDM9_haplotypes.tsv