import pysam
import collections

i = pysam.Samfile("../../mbf_align/tests/sample_data/rnaseq_spliced.bam")
input = list(i.fetch(until_eof=True))
op_sam = pysam.Samfile("spliced_reads.bam", "wb", template=i)
op_blocks = open("spliced_reads.blocks", "w")
counts = collections.Counter()
count = 0
for r in input:
    if counts[len(r.cigar)] > 10:
        continue
    counts[len(r.cigar)] += 1
    op_sam.write(r)
    op_blocks.write("let blocks = it.next().unwrap().unwrap().blocks();\n")
    ii = 0
    for start, stop in r.blocks:
        op_blocks.write("//" + r.cigarstring + " - %i \n"% count)
        op_blocks.write("assert!(blocks[%i] == (%i, %i));\n" % (ii, start, stop))
        ii += 1
    op_blocks.write("\n")
    count += 1
    if count > 100:
        print('leaving because enough')
        break
