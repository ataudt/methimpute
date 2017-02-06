setwd('~/methimpute/')
compileAttributes()
load_all()

### Test convergence of different settings ###
load('rene-data-chr1.RData')
eps <- 0.01
# data <- data[data@seqnames == 1]
data <- data[1:1e5]

## Get to work
distcor <- distanceCorrelation(data)
load_all(); fit <- estimateTransDist(distcor)
transDist <- fit$transDist
load_all(); model <- callMethylation(data, eps=1, max.iter=2, num.threads = 1, verbosity=1, fit.on.chrom=1, transDist = Inf)

