# Set up R error handling to go to stderr
options(show.error.messages=FALSE,
        error=function() {
          cat(geterrmessage(), file=stderr())
          q(save="no", status=1, runLast=FALSE)
          }
        )


library(methimpute)


filepath <- system.file("extdata", "arabidopsis_sequence.fa.gz", package="methimpute")
cytosines <- extractCytosinesFromFASTA(filepath, contexts = c('CG', 'HCCG', 'GCCG', 'CWG', 'CHH'))

## Get an example file in BSSeeker format
file <- system.file("extdata","arabidopsis_bsseeker.txt.gz", package="methimpute")
data(arabidopsis_chromosomes)
bsseeker.data <- importBSSeeker(file, chrom.lengths=arabidopsis_chromosomes)
seqlengths(bsseeker.data) <- seqlengths(cytosines)

data <- inflateMethylome(bsseeker.data, cytosines)
model <- callMethylation(data, verbosity=5, num.threads=1)
