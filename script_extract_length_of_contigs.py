from Bio import SeqIO

input_file="../Annotations/01.genome/genome.fa"
output_file = open ("Amur_contigs_lengths_total.txt",  'w')
fasta_sequences = SeqIO.parse(open(input_file),'fasta')
for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        seq_len=len(sequence)
	to_export=name+" "+str(seq_len)+"\n"
        output_file.write(to_export)

