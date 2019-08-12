#!/bin/sh

File_list=$(ls *genes-per-scaffold-count.txt)

echo "scaffold_ID scaffold_length nb_gene_background nb_gene_AnnMonocots nb_gene_basalAnnona nb_gene_Dicots" > genes_count_per_scaffold_all_topology.txt

cat Amur_contigs_lengths_total.txt | while read line
do
	scaffold_ID=$(echo $line | cut -d" " -f1)
	scaffold_length=$(echo $line | cut -d" " -f2)
	scaffold_IDg=$scaffold_ID"[[:space:]]"

	nb_gene_background=0
	nb_gene_AnnMonocots=0
	nb_gene_basalAnnona=0
	nb_gene_Dicots=0

#	echo $scaffold_IDg
	nb_gene_background=$(grep $scaffold_IDg orthologs_background_SwissProt_annotations.txt_genes-per-scaffold-count.txt |cut -d" " -f2)
	if [ -z "$nb_gene_background" ]
	then
		nb_gene_background=0
	fi
	nb_gene_AnnMonocots=$(grep $scaffold_IDg AnnMonocots_BP70.txt_SwissProt_annotations.txt_genes-per-scaffold-count.txt |cut -d" " -f2)
	if [ -z "$nb_gene_AnnMonocots" ];
	then
		nb_gene_AnnMonocots=0
	fi
	nb_gene_basalAnnona=$(grep $scaffold_IDg basal_Annona_BP70.txt_SwissProt_annotations.txt_genes-per-scaffold-count.txt |cut -d" " -f2)
	if [ -z "$nb_gene_basalAnnona" ];
	then
		nb_gene_basalAnnona=0
	fi
	nb_gene_Dicots=$(grep $scaffold_IDg Dicots_BP70.txt_SwissProt_annotations.txt_genes-per-scaffold-count.txt |cut -d" " -f2)
	if [ -z "$nb_gene_Dicots" ];
	then
		nb_gene_Dicots=0
	fi

	echo $scaffold_ID" "$scaffold_length" "$nb_gene_background" "$nb_gene_AnnMonocots" "$nb_gene_basalAnnona" "$nb_gene_Dicots >> genes_count_per_scaffold_all_topology.txt
done
