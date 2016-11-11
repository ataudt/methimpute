setwd('~/popmeth/')
compileAttributes()
load_all()

### Binning ###
files <- list.files('/media/aaron/seagate2/DATA/icb_work/popmeth/test_methylomes/', full.names=TRUE, pattern="gz$")
chrom.lengths <- GenomeInfoDb::seqlengths(Rsamtools::BamFile('~/work_ERIBA/test/popmeth_files/sample_Paw-13_1_map_paired_rmd_and_rmo.bam'))
names(chrom.lengths) <- sub('^chr', '', names(chrom.lengths))

for (file in files) {
    filename <- paste0(file,'_binsize100.RData')
    if (!file.exists(filename)) {
        data <- importRene(file, chrom.lengths = chrom.lengths)
        for (binsize in rev(c(100, 1000, 10000))) {
            filename.binsize <- paste0(file,'_binsize', binsize, '.RData')
            if (!file.exists(filename.binsize)) {
                bins <- binCounts(data, binsize=binsize)
                save(bins, file=filename.binsize)
                rm(bins)
            }
        }
        rm (data)
    }
}

### Methylation calling ###
savedir <- '~/work_ERIBA/test/popmeth_files/test_methylomes'
if (!file.exists(savedir)) { dir.create(savedir) }
files <- list.files('/media/aaron/seagate2/DATA/icb_work/popmeth/test_methylomes/', full.names=TRUE, pattern='binsize100\\.')
file <- files[1]
for (file in files) {
    filename <- file.path(savedir, paste0(basename(file)))
    if (!file.exists(filename)) {
        bins <- get(load(file))
        model <- callMethylation(bins, ID=basename(file), fit.on.chrom = 1, eps=10)
        save(model, file=filename)
    }
}
# plotHistogram(model$mhmm, binwidth=100)
# plotBoxplot(model$mhmm)
# plotHistogram(model$phmm, binwidth=100)
# plotBoxplot(model$phmm)
# plotScatter(model, datapoints=1000)

### Draw heatmap ###
files <- list.files(savedir, full.names=TRUE)
models <- loadFromFiles(files)
for (i1 in 1:length(models)) {
    # class(models[[i1]]) <- 'binnedMethylome'
    # seqlengths(models[[i1]]$segments) <- seqlengths(models[[i1]]$data)[seqlevels(models[[i1]]$segments)]
    # seqlevels(models[[i1]]$segments) <- seqlevels(models[[i1]]$data)
}
bins <- models[[1]]$data
consensus <- makeConsensus(models)
range <- trim(resize(reduce(consensus[which(consensus$variance > 2)]), fix='center', width=100e3))[20]
heatmapGenomewide(models, file='test_range.pdf', range=range)

