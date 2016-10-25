from intervaltree import Interval, IntervalTree
import matplotlib.pyplot as plt

# Parse exons
ep = []
exons = IntervalTree()

bedFile = open('../data/all_exons.bed', 'r')
circRnasPositions = {}
circRnaBlockPositions = {}

for line in bedFile:
    chr_n, start, end, _ = line.split(maxsplit=3)

    if chr_n == 'chr1':
        start, end = int(start), int(end)
        ep.append(start)
        ep.append(end)

        if end > start:
            exons[start:end] = True

ep = set(ep)
bedFile.close()

# Parse circRNAs from BED file
bedFile = open('../data/hsa_hg19_Rybak2015.bed', 'r')
cstarts = []
cends = []
csizes = []
cblocknos = []
crbp = []

for line in bedFile:
    chr_n, start, end, gene_name, _, strand, _, _, _, _, sizes, starts = line.split(maxsplit=12)

    if chr_n == 'chr1' and strand == '-':
        start, end = int(start), int(end)
        starts = [int(n) for n in starts.split(',')]
        sizes = [int(n) for n in sizes.split(',')]
        cstarts.append(start)
        cends.append(end)
        csizes.append(end-start)
        cblocknos.append(len(sizes))

        for i in range(1, len(starts) - 1):
            crbp.append(start + starts[i])
            crbp.append(start + starts[i] + sizes[i])

plt.hist(csizes, range(0, 100000, 1000))
plt.show()

plt.hist(cblocknos, range(0, 25))
plt.show()

bedFile.close()

crp = set(cstarts + cends)
crbp = set(crbp)

# Show sequence
seq = [ ('E', p) for p in ep ] + \
      [ ('R', p) for p in crp ] + \
      [ ('I', p) for p in crbp ]
seq.sort(key=lambda x: x[1])
seq.append((None, -1))

filteredSeq = []

lastPos = -1
keys = []
for key, pos in seq:
    if pos == lastPos:
        keys.append(key)
    else:
        if len(keys) > 0:
            filteredSeq.append((keys, lastPos))

            if len(filteredSeq) >= 3 and all(k == ['E'] for k, p in filteredSeq[-3:]):
                del filteredSeq[-2]

        lastPos = pos
        keys = [key]

#for keys, pos in filteredSeq[1:1000]:
#    print("%7d %s" % (pos, ''.join(keys)))

print("circRna boundaries: %d\ncircRna block boundaries: %d\nExon boundaries: %d" % (len(crp), len(crbp), len(ep)))
print("circRNA boundaries not found in exons: %d" % len(crp - ep))
print("circRNA block boundaries not found in exons: %d" % len(crbp - ep))
print("circRNA block boundaries not found in exons, ±1: %d" % len([ True for n in crbp if len(set([n,n+1]) & ep) == 0 ]))
print("Exons not found in circRNA (block) boundaries: %d" % len(ep - crp -crbp))
print("circRNA starts outside exon: %d" % len([ True for n in cstarts if exons[n] == [] ]))
print("circRNA ends outside exon: %d" % len([ True for n in cends if exons[n] == [] ]))