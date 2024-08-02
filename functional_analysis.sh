# EggNOG = mapper used for fast functional annotation, assigns func9ons to the bacterial genes identified, based on the KEGG Orthology (KO) database
# It requires the output files of Prokka (the ORFs)

conda activate eggnog

# Mapping
cd prokka/
for i in $(ls ../Fprau_genomes | rev |cut -f2- -d\. |rev) 
do 
emapper.py -i ${i}/*.faa --output ${i}_out --output_dir ../eggnog/ --temp_dir ../temp/ --cpu 10 --data_dir /mnt/fat1/databases/eggnog-mapper/data/ \
-d bact  # Target database for sequence searches
done

# The annotated files produced by EggNOG are scan to look for KOs
for i in $(ls Fprau_genomes | rev | cut -f2- -d\. | rev)
do
cat eggnog/${i}_out.emapper.annotations | cut -f12 | grep -v "#" | grep -v "KEGG_ko" | grep -v "-" >> ./all_KOs_list.txt
done
tr , '\n' < all_KOs_list.txt > all_KOs_list_nocommas.txt

# All KOs identified in the genome panel
cat all_KOs_list_nocommas.txt | uniq | sort | uniq > uniq_KOs_list.txt

# Construction of a file for each genome, containing infor regarding the presence/absence of each KO
mkdir KEGG_KOs
mkdir annotations

for i in $(ls Fprau_genomes | rev | cut -f2- -d\. | rev)
do
cat eggnog/${i}_out.emapper.annotations | cut -f12 | grep -v "#" | grep -v "KEGG_ko" | grep -v "-" > annotations/${i}.annotations
done

for i in $(ls Fprau_genomes | rev |cut -f2- -d\. |rev)
do
  NAME=${i%%.fna*};
  echo "KEGG_ko ${NAME}" > KEGG_KOs/KOs__${NAME}.txt ; # stampare a schermata
  for s in $(cat uniq_KOs_list.txt)
  do
    if grep -w "${s}" annotations/${i}.annotations ;
    then var="1" ; 
    echo "${s}   ${var}" >> KEGG_KOs/KOs__${NAME}.txt ;
    else var="0" ; 
    echo "${s}   ${var}" >> KEGG_KOs/KOs__${NAME}.txt ;
    fi
  done
done

# Presence/absence table of the KOs
cd KEGG_KOs/
for i in $(ls *.txt)
do
	cat ${i} | cut -f2 > ${i}_column_2.txt # il simbolo ">" indica che sovrascrive, mentre ">>" indica che aggiunge 
done
paste uniq_KOs_list.txt KEGG_KOs/*column_2.txt > KOs_binary_table.tsv
