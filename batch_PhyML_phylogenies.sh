# cd "4_species_orthologs_longer100aa_aligned"

for a in *.fasta; do
    perl fasta2phylip.pl $a
done
 

for i in *.phy; do
    sed -i 's/Atri_/Atri*_/g' $i
#    sed -i 's/Amur*_/Amur_/g' $i
done


for gene in *.phy
do
	echo "Building tree for "$gene"...."
	PhyML-3.1_linux64 -i $gene -d aa -b -4 -m WAG -f m -v e -c 4 -a e -o tlr --quiet --no_memory_check
	echo "ML tree built for "$gene"."
	echo "---"
done

rename 's/.phy_phyml_tree.txt/_tree.tre/g' *
