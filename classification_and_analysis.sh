# GTDB-Tk = software toolkit for assigning objective taxonomic classifications to bacterial and archaeal genomes based on the Genome Database Taxonomy

conda activate prokka 

# Names of the samples used in our study
ID="SRR12557704 SRR12557705 SRR12557706 SRR12557707 SRR12557708 SRR12557709 SRR12557710 SRR12557711 SRR12557712 SRR12557713 SRR12557714 
SRR12557715 SRR12557716 SRR12557717 SRR12557718 SRR12557719 SRR12557720 SRR12557721 SRR12557722 SRR12557723 SRR12557724 SRR12557725 
SRR12557726 SRR12557727 SRR12557728 SRR12557729 SRR12557730 SRR12557731 SRR12557733 SRR12557734"

for i in ${ID}
do
gtdbtk classify_wf --genome_dir ${i}/metawrap_50_5_bins/ \
--out_dir ${i}/GTDBTK/ \
-x fa \
--cpus 10 \
--mash_db /mnt/fat1/databases/gtdbtk/release214/mash/gtdb_ref_sketch.msh
done

# Select the bins containing Faecalibacterium prausnitzii
cat *GTDBTK/gtdbtk.bac120.summary.tsv | grep "Faecalibacterium"

# Select the files containing the selected bins'
for i in ${ID}
do
grep "Faecalibacterium" ${i}/GTDBTK/gtdbtk.bac120.summary.tsv
echo "${i}"
done

cp  SRR12557728/metawrap_50_5_bins/bin.28.fa ../../Fprau_genomes/

for i in $(ls bin* |cut -f1,2 -d\.); do mv ${i}.fa Wibowo_${i}.fna; done


# Analysis of the genome panel - ANI pairwise distances 

conda activate pyani
average_nucleotide_identity.py -i Fprau_genomes -o pyani_output_fprau \
-m ANIb \  # Uses BLASTN+ for 1020nt alignments
-g --workers 10 -v  


# Analysis of the genome panel - Annotation with Prokka

conda activate prokka

mkdir prokka 
for i in $(ls ../Fprau_genomes | rev |cut -f2- -d\. |rev) 
do 
prokka ../Fprau_genomes/${i}.fna -o ./${i} --fast --locustag ${i} --cpus 10 \
--force  # Force overwriting existing output folder
done 

# Group together the GFF files
mkdir GFF
for i in $(ls Fprau_genomes | rev |cut -f2- -d\. |rev) 
do 
cp prokka/${i}/*.gff GFF/${i}.gff
done


# Analysis of the genome panel - Pangenome by ROARY

conda activate roary

mkdir roary
cd roary
roary -i 90 \  # Minimum percentage identity BLASTP
-cd 90 \  # Percentage of genomes in which a gene has to be present in order to classify as “core”
-e \  # PRANK for multiFASTA gene alignment
-p 10 \ # Number of threads
../GFF/*.gff --group_limit 100000 -v
