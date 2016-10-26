from intervaltree import IntervalTree
from parse import *
import random
import itertools
import time
import numpy as np
from pyfasta import Fasta

print("Loading files...")
t0 = time.time()

seqs = Fasta('data/hg19.fa')
crnas = read_bed('data/hsa_hg19_Rybak2015.bed')
exons = read_bed('data/all_exons.bed')

t1 = time.time()
print("%.2fs" % (t1 - t0))
print("Finding free ranges...")
t0 = time.time()

# Constants for scanline points
CRNA, EXON = 0, 1
START, END = 0, 1

# Group by chromosome
chromosomes = [ 'chr%d' % i for i in range(1, 23) ] + [ 'chrX', 'chrY' ]
points = {}
free_ranges = {}
free_exons = []

# Only use useful chromosomes
for c in chromosomes:
	points[c], free_ranges[c] = [], []

for crna in crnas:
	if crna.chr_n in chromosomes:
		points[crna.chr_n].append((CRNA, START, crna.start, crna))
		points[crna.chr_n].append((CRNA, END, crna.end, crna))

for exon in exons:
	if exon.chr_n in chromosomes:
		points[exon.chr_n].append((EXON, START, exon.start, exon))
		points[exon.chr_n].append((EXON, END, exon.end, exon))

# Sort by position, then ends before starts (False < True)
for c in chromosomes:
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
				free_ranges[c].append(( current_range[0], current_range[-1] ))
				current_range = []
		elif geneType == CRNA and pointType == END:
			# Decrease current cRNA count
			current_crnas -= 1
		elif geneType == EXON and pointType == START and current_crnas == 0:
			# Add exon to list if no cRNA is covering it
			current_exons.append(gene)
			free_exons.append(gene)
		elif geneType == EXON and pointType == END and current_crnas == 0:
			# Add exon to working range if no cRNA is covering it
			if gene in current_exons:
				current_range.append(gene)
				current_exons.remove(gene)

t1 = time.time()
print("%.2fs" % (t1 - t0))
print("Generating negative examples...")
t0 = time.time()

bedFile = open('tmp/negatives.bed', 'w')

not_crnas = []
ok, failed = 0, 0
for crna in crnas:
	# Choose length
	length = crna.end - crna.start

	# Use same chromosome distribution FIXME?
	start, end = crna.start, crna.end
	chr_n = crna.chr_n

	# Only use useful chromosomes
	if chr_n not in chromosomes: continue

	# Find spot with exactly that length
	suitable_ranges = [ 
		r for r in free_ranges[chr_n] 
		if r[1].end - r[0].start >= length 
		and r[1].start - r[0].end < length ]

	if len(suitable_ranges) == 0:
		failed += 1
		continue
	else:
		ok += 1

	r = random.choice(suitable_ranges)
	start = random.randint(max(r[0].start, r[1].start - length),
		min(r[0].end, r[1].end - length))
	end = start + length

	# Choose strand
	strand = random.choice('+-')

	# Write to file
	print("%s\t%9d\t%9d\t%s" % (chr_n, start, end, strand), file=bedFile)

bedFile.close()

# Write unoverlapped exons to file (primitive negative samples)
exonFile = open('tmp/free_exons.bed', 'w')
for exon in free_exons:
	print("%s\t%9d\t%9d\t%s" % (exon.chr_n, exon.start, exon.end, exon.strand), file=exonFile)

exonFile.close()

print("%d %d" % (ok, failed))

t1 = time.time()
print("%.2fs" % (t1 - t0))