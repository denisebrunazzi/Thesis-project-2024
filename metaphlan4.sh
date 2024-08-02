# MetaPhlAn4 = computatioonal tool which identifies the microbiota species or strains present in a sample, based on specific unique genomic sequences

conda activate metaphlan4

mkdir metaphlan/

# Generates a .tsv file containing the abundances of specific species and strains inside the sample.
for i in $(cut -f7 PRJNA561510_tsv.txt | sed 1d | cut -f1 -d\; | cut -f6 -d/)  # PRJNA561510_tsv.txt contains the FASTA sequences
do
metaphlan sequences/${i}_1.fastq,sequences/${i}_2.fastq \
--nproc 10 \   # Uses 10 processing to be faster
--input_type fastq \
-o metaphlan/${i}_metagenome.tsv \
-s metaphlan/${i}_metagenome.sam \
--bowtie2out metaphlan/${i}.bowtie2.bz2  # specify this when you have multiple sequences to peed up the procedure
done

# The .tsv files are merged together in a table containing the sequence reads in the columns, and all the identified bacterial species and subspecies in the rows.
merge_metaphlan_tables.py SRR125577*.tsv > Wibowo_metaphlan4.tsv
