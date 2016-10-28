from intervaltree import IntervalTree
from util import *
import random
import itertools
import time
import sys
import numpy as np
from pyfasta import Fasta

tick("Loading files...")

seqs = Fasta('data/hg19.fa')
crnas = read_bed('data/hsa_hg19_Rybak2015.bed')
exons = read_bed('data/all_exons.bed')

tick("Finding free ranges...")

# Constants for scanline points
CRNA, EXON = 0, 1
START, END = 0, 1

# Group by chromosome
points = {}
free_ranges = []
free_exons = []

# Prepare scanline points
for gene in crnas:
	points.setdefault(gene.chr_n, []).extend([
		(CRNA, START, gene.start, gene),
		(CRNA, END, gene.end, gene)])

for gene in exons:
	points.setdefault(gene.chr_n, []).extend([
		(EXON, START, gene.start, gene),
		(EXON, END, gene.end, gene)])

for c in points.keys():
	# Sort by position, then ends before starts (False < True)
	points[c] = sorted(points[c], \
		key=lambda p: (p[2], p[1] == START))

	# Add cRNA start at the end for scanline algorithm
	points[c].append((CRNA, START, points[c][-1][2], None))

	# Scanline to find free ranges
	current_crnas = 0
	current_exons, current_range = [], []

	for geneType, pointType, _, gene in points[c]:
		if geneType == CRNA and pointType == START:
			# Add working range to list, close all open exons
			current_crnas += 1
			current_exons = []

			if len(current_range) > 0:
				free_ranges.append(( current_range[0], current_range[-1] ))
				current_range = []
		elif geneType == CRNA and pointType == END:
			# Decrease current cRNA count
			current_crnas -= 1
		elif geneType == EXON and pointType == START and current_crnas == 0:
			# Add exon to list if no cRNA is covering it
			current_exons.append(gene)
		elif geneType == EXON and pointType == END and current_crnas == 0:
			# Add exon to working range if no cRNA is covering it
			if gene in current_exons:
				current_range.append(gene)
				free_exons.append(gene)
				current_exons.remove(gene)

tick("Generating negative examples...")

# Generate negatives from multiple exons, using the same length and 
# chromosome distribution as the positives
bedFile = open('tmp/negatives.bed', 'w')

not_crnas = []
ok, failed = 0, 0
for i, crna in enumerate(crnas):
	# Choose length
	length = crna.end - crna.start

	start, end = crna.start, crna.end

	# Find spot with exactly that length
	r = None
	for i in range(0, len(free_ranges)):
		r = random.choice(free_ranges)
		if r[1].end - r[0].start >= length \
		and r[1].start - r[0].end < length:
			break

	# Skip if no range with that length can be found
	if r is None:
		continue

	start = random.randint(max(r[0].start, r[1].start - length),
		min(r[0].end, r[1].end - length))
	end = start + length

	# Choose strand
	chr_n = r[0].chr_n
	strand = random.choice('+-')

	# Write to file
	print("%s\t%9d\t%9d\t%s" % (chr_n, start, end, strand), file=bedFile)

	if i % 1000 == 0:
		print('.', end='')
		sys.stdout.flush()


bedFile.close()

# Write unoverlapped exons to file (primitive negative samples)
exonFile = open('tmp/free_exons.bed', 'w')
for exon in free_exons:
	print("%s\t%9d\t%9d\t%s" % (exon.chr_n, exon.start, exon.end, exon.strand), file=exonFile)

exonFile.close()

tick()