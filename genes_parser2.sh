#!/bin/sh

#  genes_parser.sh
#  
#
#  Created by Damien HINSINGER on 26/01/2018.
#

echo "genes_parser script, by Damien"
echo ""
echo "*********************************************************************************************************"
#echo "How to use this script :"
#echo ""
#echo "- From Geneious, export both a fasta file (sequence) and a gff file (annotations, without the sequence)"
#echo "  for each species of interest, and put all files in a folder"
#echo "- Copy this script into this folder, copy also the perl script 'gff2perl_modified.pl' in the same folder"
#echo "- Type 'sh gff2alignments.sh' in a terminal window (being in the folder)"
#echo "- Wait and enjoy a coffee"
#echo "- The results are in the 'Genes_unaligned' folder"
#echo "*********************************************************************************************************"
#echo ;echo;echo;echo

# cd fixed_alignment_file
rm -r genes_parsed_2
mkdir genes_parsed_2

cat all.muscle.pep.new | while read linefic1
	do
#		echo "line"
		if ! echo $linefic1 | egrep -q '^ *>'
		then
			# that is a sequence
			# echo $linefic1
			species=$(echo $linefic1 | cut -d" " -f1 | cut -f1 | rev | cut -d"_" -f1 | rev )
			
			sequence=$(echo $linefic1 | cut -d" " -f2)
			echo ">"$species"_"$linefic1 >>"genes_parsed_2/"$gene".pep_aligned.fasta"
			echo $sequence >>"genes_parsed_2/"$gene".pep_aligned.fasta"
			#echo $sequence >> "genes_parsed/"$gene".pep_aligned.fasta"
#			if [ "$gene" = "$gene_old" ]
			#echo $species
#			then
#				echo "$gene showed several exons - they were merged"
#			else
#			# the sequence's name
        	else

			# Gene name found
			gene=$(echo $linefic1 | cut -d">" -f2)
			echo "Dealing with gene "$gene
						fi
	done
echo "#####################################################################"









cd genes_parsed_2

echo "Converting fasta files to phylip for phyML"
for gene in *.pep_aligned.fasta
do
	echo "converting "$gene
	perl ../fasta2phylip.pl $gene
	echo $gene" converted..."
done
echo "#######################################"


for i in *.phy; do
    sed -i 's/Atr_/Atr*_/g' $i
#    sed -i 's/Amur*_/Amur_/g' $i
done


# run the ML tree using PhyML

for gene in *.pep_aligned.fasta.phy
do
	echo "Building tree for "$gene"...."
	PhyML-3.1_linux64 -i $gene -d aa -b -4 -m WAG -f m -v e -c 4 -a e -o tlr --quiet --no_memory_check
	echo "ML tree built for "$gene"."
	echo "---"
done


## Still to do :

# move the stat files to an other folder

# format the tree files to go in R


# cd ../



