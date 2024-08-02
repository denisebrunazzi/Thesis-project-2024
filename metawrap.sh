# MetaWRAP = modular pipeline for shotgun metagenomic data analysis, used for the assembly of high-quality draft bacterial genomes.

conda activate metawrap-env

# Names of the samples used in our study
ID="SRR12557704 SRR12557705 SRR12557706 SRR12557707 SRR12557708 SRR12557709 SRR12557710 SRR12557711 SRR12557712 SRR12557713 SRR12557714 
SRR12557715 SRR12557716 SRR12557717 SRR12557718 SRR12557719 SRR12557720 SRR12557721 SRR12557722 SRR12557723 SRR12557724 SRR12557725 
SRR12557726 SRR12557727 SRR12557728 SRR12557729 SRR12557730 SRR12557731 SRR12557733 SRR12557734"

# Assembly 
mkdir assembly/
for i in ${ID}
do
metawrap assembly -1 sequences/${i}_1.fastq -2 sequences/${i}_2.fastq \ # Requires two FASTA files resulting from the shotgun sequencing of each sample
-m 200 \  # Memory in GB (default=10)
-t 10 \  # Number of threads (default=1)
--megahit \  # Memory efficient, faster, and scales well with large datasets
-o assembly/${i}_assembly_megahit
done

# Binning
mkdir binning/
for i in ${ID}
do
metawrap binning -o binning/${i} \
-t 30 \
-a assembly/${i}_assembly_megahit/final_assembly.fasta \
--metabat2 --maxbin2 --concoct sequences/${i}_1.fastq sequences/${i}_2.fastq  # Uses three diffeerent metagenomic binning softwares
done

# Bin refinement
mkdir bin_refinement/
for i in ${ID}
do
metawrap bin_refinement -o bin_refinement/${i} \
-t 30 \
-A binning/${i}/maxbin2_bins -B binning/${i}/metabat2_bins -C binning/${i}/concoct_bins \ # Each binning softwares is assigned a letter
-c 50 \  # Completeness
-x 5  # Contamination
done
