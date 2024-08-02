# PhyloPhlAn = pipeline for large-scale phylogenetic profiling of genomes and metagenome.

conda activate phylophlan

# Custom configuration file
phylophlan_write_config_file \
-o custom_config_fprau_nt.cfg \
-d a \  # Type of database (uses 136 universal marker genes)
--db_aa diamond --map_dna diamond --map_aa diamond \  # BLAST Diamond for the mapping step
--msa mafft \  # MAFFT for the mul9ple sequence alignment
--trim trimal \  # trimAl for trimming
--tree1 fasttree \  # FastTree to build the first tree
--tree2 raxml \  # RAxML to refine the final tree
--overwrite

# Tree building
phylophlan \
-i ../Fprau_genomes/ \
-d phylophlan \
--diversity low \
-f custom_config_fprau_nt.cfg \
--nproc 10 \
--fast
