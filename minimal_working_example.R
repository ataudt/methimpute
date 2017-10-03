setwd('~/methimpute/')
# compileAttributes()
load_all()

### Run contexts separately ###
fit.on.chrom=NULL; transDist=Inf; eps=1; max.time=Inf; max.iter=Inf; count.cutoff=500; verbosity=1; initial.params=NULL; include.intermediate=FALSE; update='independent'; min.reads=0; num.threads=2+include.intermediate



### Import Poplar methylomes ###
file <- '~/work_ERIBA/test/methimpute_files/allc_methyl-13-1-1.tsv'
file <- '~/work_ERIBA/test/methimpute_files/test.tsv'
methylome <- importMethylpy(file)




### Make cytosine positions ###
file <- '~/work_ERIBA/test/methimpute_files/TAIR10_chr_all.fa'
cs <- extractCytosinesFromFASTA(file)
ms <- get(load('~/work_ERIBA/test/methimpute_files/cyt_positions_TAIR10.RData'))

cs1 <- cs[seqnames(cs) == 'chr1']
ms1 <- ms[seqnames(ms) == 'chr1']
