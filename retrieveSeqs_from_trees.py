#! #!/usr/bin/env python

import re
import os
import sys
from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline




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

#########################


#### We start organising some folders to put the data
_mkdir("genes_from_orthologs_analysis_BP0")
_mkdir("genes_from_orthologs_analysis_BP50")




for tree_names in ("AnnMonocots_BP0.txt", "Dicots_BP0.txt", "basal_Annona_BP0.txt", "AnnMonocots_BP50.txt", "Dicots_BP50.txt", "basal_Annona_BP50.txt"):
	File1 = open(tree_names, 'rU')
	Outfile_name=tree_names+"_seqs_names.txt"
	Outfile_name2=tree_names+"_SwissProt_annotations.txt"
	OutFile = open (Outfile_name,  'w')
	OutFile2 = open (Outfile_name2,  'w')
	for Line in File1:
		group=Line.split("_")	
		# print group
		#### First for the orthologs among all the 4 species
		file_seq="4_species_orthologs_aligned_aa/"+str(group[0])+".fasta"
		for line in open(file_seq):
			if ">Amur_" in line:
				gene=line.split("_",1)
				#print gene
				OutFile.write(gene[1])
				gene4grep0=gene[1].split("\n")
				gene4grep=gene4grep0[0]+"\t"
				#print gene4grep
				for annotation in open("../Annotations_orthologs_Atha-Osat-Atri-Amur/blast2go_gaf_20180306_0030_annotations_GO-slim.txt"):
					if gene4grep in annotation:
						#print annotation
						OutFile2.write(annotation)
	File1.close()
	OutFile.close()
	OutFile2.close()

###############
