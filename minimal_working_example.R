setwd('~/popmeth/')
compileAttributes()
load_all()

### Test convergence of different settings ###
load('rene-data-chr1.RData')

## 2-state HMM, min coverage = 0
model <- callMethylationBinomial(data, min.reads=0, include.heterozygosity = FALSE)
# CONVERGENCE FINE AT ITERATION 21

## 3-state HMM, min coverage = 0
model <- callMethylationBinomial(data, min.reads=0, include.heterozygosity = TRUE)
# NEGATIVE CHANGE AT ITERATION 16

## 2-state HMM, min coverage = 3
model <- callMethylationBinomial(data, min.reads=3, include.heterozygosity = FALSE)
# NEGATIVE CHANGE AT ITERATION 9

## 3-state HMM, min coverage = 3
model <- callMethylationBinomial(data, min.reads=3, include.heterozygosity = TRUE)
# NEGATIVE CHANGE AT ITERATION 7

## 2-state HMM, min coverage = 0, context-specific
model <- callMethylationBinomialContext(data, min.reads=0, include.heterozygosity = FALSE)
# NEGATIVE CHANGE AT ITERATION 18

## 3-state HMM, min coverage = 0, context-specific
model <- callMethylationBinomialContext(data, min.reads=0, include.heterozygosity = TRUE)
# NEGATIVE CHANGE AT ITERATION 12

## 2-state HMM, min coverage = 3, context-specific
model <- callMethylationBinomialContext(data, min.reads=3, include.heterozygosity = FALSE)
# CONVERGENCE FINE AT ITERATION 33

## 3-state HMM, min coverage = 3, context-specific
model <- callMethylationBinomialContext(data, min.reads=3, include.heterozygosity = TRUE)
# NEGATIVE CHANGE AT ITERATION 16
