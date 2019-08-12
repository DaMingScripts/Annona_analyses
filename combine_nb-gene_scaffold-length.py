#! #!/usr/bin/env python

import re
import os
import sys

Usage = """
Combine the number of genes found in each scaffold for each case with the scaffolds lengths table.
Usage :
	combine_nb-gene_scaffold-length.py *genes-per-scaffold-count.txt 
	(input files are the result files of the script script_extract_length_of_contigs.py)
"""


FileList =sys.argv[1:]
for InfileName in FileList:
	print InfileName


FileNum=0
MasterList =[]
for InfileName in FileList:
		Infile = open(InfileName, 'r')
		LineNumber = 0
		RecordNum = 0
		for Line in Infile:
			Line=Line.strip('\n')
			if FileNum == 0:
				MasterList.append(Line)
			else:
				ElementList = Line.split(' ')
				print ElementList[1]
				print MasterList[RecordNum]
				MasterList[RecordNum] += (' '+ElementList[1])
				RecordNum += 1
			LineNumber +=1
		Infile.close()
		FileNum +=1
