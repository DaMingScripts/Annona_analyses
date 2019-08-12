#! #!/usr/bin/env python

import re
import os
import sys
from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline


#### look in a fasta file for a specific sequence
def getSeqRec(seq_id,species_file):
	record=""
        fasta=SeqIO.parse(species_file,"fasta")
        for record in fasta:
                if seq_id==record.id:
			return record
##


#### Create new directory
def _mkdir(newdir):
    """works the way a good mkdir should :)
        - already exists, silently complete
        - regular file in the way, raise an exception
        - parent directory(ies) does not exist, make them as well
    """
    if os.path.isdir(newdir):
        pass
    elif os.path.isfile(newdir):
        raise OSError("a file with the same name as the desired " \
                      "dir, '%s', already exists." % newdir)
    else:
        head, tail = os.path.split(newdir)
        if head and not os.path.isdir(head):
            _mkdir(head)
        #print "_mkdir %s" % repr(newdir)
        if tail:
            os.mkdir(newdir)
##


#### We start organising some folders to put the data
_mkdir("pairwise_orthologs")
_mkdir("4_species_orthologs")
_mkdir("pairwise_orthologs_aligned")
_mkdir("4_species_orthologs_aligned")


orthomcl_output = open("mclOutput", 'rU')

OutFile4orth = open ("temp.4orth.txt",  'w')

group=0

#### Parsing the orthoMCL output file to retrieve the species and genes names
for Line in orthomcl_output:
	genes=Line.split()
	#### First for the orthologs among all the 4 species
	if len(genes)==4:
		gene1=genes[0].split('|')
		sp_gene1=gene1[0]
		seq_gene1=gene1[1]
		if sp_gene1=="Amur":
			file1="Amuricata_CDS_nc.fasta"
		elif sp_gene1=="Atri":
			file1="Amborella_trichopoda_CDS_nc.fasta"
		elif sp_gene1=="Atha":
			file1="Arabidopsis_thaliana_CDS_nc.fasta"
		elif sp_gene1=="Osat":
			file1="Oryza_sativa_CDS_nc.fasta"

		gene2=genes[1].split('|')	
		sp_gene2=gene2[0]
		seq_gene2=gene2[1]
		if sp_gene2=="Amur":
			file2="Amuricata_CDS_nc.fasta"
		elif sp_gene2=="Atri":
			file2="Amborella_trichopoda_CDS_nc.fasta"
		elif sp_gene2=="Atha":
			file2="Arabidopsis_thaliana_CDS_nc.fasta"
		elif sp_gene2=="Osat":
			file2="Oryza_sativa_CDS_nc.fasta"
		gene3=genes[2].split('|')
		sp_gene3=gene3[0]
		seq_gene3=gene3[1]
		if sp_gene3=="Amur":
			file3="Amuricata_CDS_nc.fasta"
		elif sp_gene3=="Atri":
			file3="Amborella_trichopoda_CDS_nc.fasta"
		elif sp_gene3=="Atha":
			file3="Arabidopsis_thaliana_CDS_nc.fasta"
		elif sp_gene3=="Osat":
			file3="Oryza_sativa_CDS_nc.fasta"

		gene4=genes[3].split('|')	
		sp_gene4=gene4[0]
		seq_gene4=gene4[1]
		if sp_gene4=="Amur":
			file4="Amuricata_CDS_nc.fasta"
		elif sp_gene4=="Atri":
			file4="Amborella_trichopoda_CDS_nc.fasta"
		elif sp_gene4=="Atha":
			file4="Arabidopsis_thaliana_CDS_nc.fasta"
		elif sp_gene4=="Osat":
			file4="Oryza_sativa_CDS_nc.fasta"

		#### if all the species are different, we have 1 ortholog per species
		if (sp_gene1 != sp_gene2) and (sp_gene4 != sp_gene3) and (sp_gene2 != sp_gene4) and (sp_gene1 != sp_gene3) and (sp_gene1 != sp_gene4) and (sp_gene2 != sp_gene3):
			unaligned_file="4_species_orthologs/group" + str(group)+".fasta"
			#### write gene_name 1 and sequence 1
			OutFile4orth = open (unaligned_file,  'w')
			seq_rec1=getSeqRec(seq_gene1,file1)
			OutFile4orth.write(">" + sp_gene1 + "_" + seq_rec1.id + "\n" + str(seq_rec1.seq+"\n"))
			OutFile4orth.close()

			#### write gene_name 2 and sequence 2
			OutFile4orth = open (unaligned_file,  'a')
			seq_rec2=getSeqRec(seq_gene2,file2)
			OutFile4orth.write(">" + sp_gene2 + "_" + seq_rec2.id + "\n" + str(seq_rec2.seq+"\n"))
			OutFile4orth.close()

			#### write gene_name 3 and sequence 3
			OutFile4orth = open (unaligned_file,  'a')
			seq_rec3=getSeqRec(seq_gene3,file3)
			OutFile4orth.write(">" + sp_gene3 + "_" + seq_rec3.id + "\n" + str(seq_rec3.seq+"\n"))
			OutFile4orth.close()

			#### write gene_name 4 and sequence 4
			OutFile4orth = open (unaligned_file,  'a')
			seq_rec4=getSeqRec(seq_gene4,file4)
			OutFile4orth.write(">" + sp_gene4 + "_" + seq_rec4.id + "\n" + str(seq_rec4.seq+"\n"))
			OutFile4orth.close()



			#### Alignment of each pair of ortholog using mafft, with direction adjustment
			mafft_cline = MafftCommandline(input=unaligned_file)
			mafft_cline.adjustdirection = True
			stdout, stderr = mafft_cline()
			aligned_file="4_species_orthologs_aligned/group" + str(group)+".fasta"
			outfile = open (aligned_file, 'w')
			outfile.write(stdout)
			outfile.close()

			
	#### Then for each pair of orthologs
	elif len(genes)==2:

		gene1=genes[0].split('|')
		sp_gene1=gene1[0]
		seq_gene1=gene1[1]
		#i=0
		if sp_gene1=="Amur":
			file1="Amuricata_CDS_nc.fasta"
		elif sp_gene1=="Atri":
			file1="Amborella_trichopoda_CDS_nc.fasta"
		elif sp_gene1=="Atha":
			file1="Arabidopsis_thaliana_CDS_nc.fasta"
		elif sp_gene1=="Osat":
			file1="Oryza_sativa_CDS_nc.fasta"

		gene2=genes[1].split('|')	
		sp_gene2=gene2[0]
		seq_gene2=gene2[1]
		if sp_gene2=="Amur":
			file2="Amuricata_CDS_nc.fasta"
		elif sp_gene2=="Atri":
			file2="Amborella_trichopoda_CDS_nc.fasta"
		elif sp_gene2=="Atha":
			file2="Arabidopsis_thaliana_CDS_nc.fasta"
		elif sp_gene2=="Osat":
			file2="Oryza_sativa_CDS_nc.fasta"

		#### if the species are different, it means it's an ortholog
		if (sp_gene1 != sp_gene2):
			unaligned_file="pairwise_orthologs/group" + str(group)+".fasta"

			OutFile2orth = open (unaligned_file,  'w')
			seq_rec1=getSeqRec(seq_gene1,file1)
  			to_write1=">" + sp_gene1 + "_" + seq_rec1.id + "\n"
			OutFile2orth.write(to_write1)
			OutFile2orth.close()

			OutFile2orth = open (unaligned_file,  'a')
 			to_write2=seq_rec1.seq+"\n"
			OutFile2orth.write(str(to_write2))
			OutFile2orth.close()

			OutFile2orth = open (unaligned_file,  'a')
			seq_rec2=""
			seq_rec2=getSeqRec(seq_gene2,file2)
			OutFile2orth.write(">" + sp_gene2 + "_"+ seq_rec2.id + "\n")
			OutFile2orth.close()
			OutFile2orth = open (unaligned_file,  'a')
 			to_write4=seq_rec2.seq+"\n"
			OutFile2orth.write(str(to_write4))
			OutFile2orth.close()

			#### Alignment of each pair of ortholog using mafft, with direction adjustment
			mafft_cline = MafftCommandline(input=unaligned_file)
			mafft_cline.adjustdirection = True
			stdout, stderr = mafft_cline()
			aligned_file="pairwise_orthologs_aligned/group" + str(group)+".fasta"
			outfile = open (aligned_file, 'w')
			outfile.write(stdout)
			outfile.close()



			
	group +=1



OutFile4orth.close()

orthomcl_output.close()

###############
