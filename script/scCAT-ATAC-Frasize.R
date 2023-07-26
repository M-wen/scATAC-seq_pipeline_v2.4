### Get the parameters
parser = argparse::ArgumentParser(description="Script to plot fragment size")
parser$add_argument('-B','--bam', help='bam file')
parser$add_argument('-L','--lab', help='the label')
parser$add_argument('-O','--out', help='the out put pdf')
args = parser$parse_args()

###
library(ATACseqQC)
source("/jdfssz1/ST_SUPERCELLS/P21Z10200N0090/Automated/USER/songyumo/pipeline/scATAC//script/fragSizeDist.R")
library(Rsamtools)
svg(args$out,width = 6,height=5)
bamfile=args$bam
bamfile.labels=args$lab
fragSize <- fragSizeDist(bamfile, bamfile.labels,color="blue")
dev.off()
