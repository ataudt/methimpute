setwd('~/methimpute/')
# compileAttributes()
load_all()

### Make cytosine positions ###
file <- '~/work_ERIBA/test/methimpute_files/TAIR10_chr_all.fa'
cs <- extractCytosinesFromFASTA(file)
ms <- get(load('~/work_ERIBA/test/methimpute_files/cyt_positions_TAIR10.RData'))

cs1 <- cs[seqnames(cs) == 'chr1']
ms1 <- ms[seqnames(ms) == 'chr1']
