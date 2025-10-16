# Directory to check PRMD9 haplotypes

Catalog: https://pmc.ncbi.nlm.nih.gov/articles/PMC8600002/

Other paper that does this:
https://www.nature.com/articles/s41586-025-09540-8


## Goal
Compare PRMD9 haplotype of 15q patients to "expected" haplotype distribution. Do we find anything that is unexpected?

## Strategy
1. Retrieve PRMD9 sequence per haplotype
Do this using samtools.
Get Region using samtools : "chr5:23526573-23527864"

2. minimap align to catalog (supp table 3)

