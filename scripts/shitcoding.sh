
INPUT_BAM="/ifs/data/research/revio/work/DNA18-10561/GRCh38_20240122/DNA18-10561.haplotagged.bam"
OUTPUT_BAM_PRDM9="/ifs/data/research/projects/michelle/17_prmd9_haplotypecheck/data/genotype/DNA18-10561.PRDM9.bam"

# Get GeneRegion using samtools
REGION="chr5:23526573-23527864"
samtools view -b -h ${INPUT_BAM} ${REGION} > ${OUTPUT_BAM_PRDM9}

#split into hap
samtools view -d HP:1 -o DNA18-10561.HP1.PRDM9.bam DNA18-10561.PRDM9.bam
samtools view -d HP:2 -o DNA18-10561.HP2.PRDM9.bam DNA18-10561.PRDM9.bam

#bedtools to find reads spanning whole region
bedtools intersect -f 1.0 -wb -a prdm9_target.bed -b DNA18-10561.HP1.PRDM9.bam  | cut -f7 | head -n 5 > spanning_readnames_H1.txt
bedtools intersect -f 1.0 -wb -a prdm9_target.bed -b DNA18-10561.HP2.PRDM9.bam  | cut -f7 | head -n 5 > spanning_readnames_H2.txt

# samtools and extract reads
samtools fasta DNA18-10561.HP1.PRDM9.bam > DNA18-10561.HP1.PRDM9.fasta
grep -f spanning_readnames_H1.txt DNA18-10561.HP1.PRDM9.fasta -A 1 | grep -v "\-\-" > targetreads_h1.fasta

# minimaap 

# make for loop to increment head lines
for (int i = 1; i < 6; i++) {

}

minimap2 -c -x map-hifi --secondary=no --paf-no-hit <(head -n 2 targetreads_h1.fasta | tail -n 2) ../prdm9_ref.fa  > read1.H1_prdm9motifs.paf
minimap2 -c -x map-hifi --secondary=no --paf-no-hit <(head -n 4 targetreads_h1.fasta | tail -n 2) ../prdm9_ref.fa > read2.H1_prdm9motifs.paf


awk '!seen[$1]++' read1.H1_prdm9motifs.paf > read1.H1_prdm9motifs_unique.paf
awk '!seen[$1]++' read2.H1_prdm9motifs.paf > read2.H1_prdm9motifs_unique.paf

# removing duplicates from PAF

# formt output
awk 'BEGIN {FS=OFS="\t"} ($2==$11) {print $1, "YES", $13} ($2!=$11) {print $1, "NO", $13}' read1.H1_prdm9motifs_unique.paf | sed 's/NM:i://g' | sort > read1.H1.prdm9matches.tsv
awk 'BEGIN {FS=OFS="\t"} ($2==$11) {print $1, "YES", $13} ($2!=$11) {print $1, "NO", $13}' read2.H1_prdm9motifs_unique.paf | sed 's/NM:i://g' | sort > read2.H1.prdm9matches.tsv


# merging per hapltype
#aggregate NM for reads together + how many "YES"


# paste command for combining paf outputs (processed) (one combined per haplotypes)

#output is super main witg 2 lines per sample H1 H2
# sample hap results confidence (#YES +NM)

# input is full bam file (haplotagged), output is the 2 lines for main result
#runtime less than half a minute

#folders
# input scripts results README.md