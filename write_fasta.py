from util import *
import random
import itertools
import time
from pyfasta import Fasta
import sys

seqs = Fasta('data/hg19.fa')
crnas = read_bed('data/hsa_hg19_Rybak2015.bed')
negs = read_bed('tmp/negatives.bed')

# Write CRNA and negative sample content to fasta file
# for further processing
crnaFile = open('tmp/crnas.fa', 'w')
for gene in crnas:
	print(">%s" % gene.gene_name, file=crnaFile)
	end = min(gene.end, gene.start + 2047)
	print(seqs[gene.chr_n][gene.start:end] + '\n', file=crnaFile)
crnaFile.close()

negFile = open('tmp/negs.fa', 'w')
for gene in negs:
	print(">%s" % gene.gene_name, file=negFile)
	end = min(gene.end, gene.start + 2047)
	print(seqs[gene.chr_n][gene.start:end] + '\n', file=negFile)
negFile.close()
