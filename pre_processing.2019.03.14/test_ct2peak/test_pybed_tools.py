import pybedtools

peaks = pybedtools.BedTool("peak.bed")
reads = pybedtools.BedTool("foo.bam")
for line in peaks.intersect(reads, wa=True, wb=True):
  elems = str(line).split()
  print(elems)
