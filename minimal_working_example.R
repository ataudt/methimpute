setwd('~/methimpute/')
# compileAttributes()
load_all()

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
