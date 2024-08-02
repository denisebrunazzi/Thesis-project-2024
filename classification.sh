# GTDB-Tk = software toolkit for assigning objective taxonomic classifications to bacterial and archaeal genomes based on the Genome Database Taxonomy

conda activate /mnt/mini1/work/daniel/miniconda3/envs/prokka 

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
