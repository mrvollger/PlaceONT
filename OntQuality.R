#!/usr/bin/env Rscript
#lib="/home/mrvollger/anaconda3/lib/R/library"
#.libPaths(lib)
#.libPaths( c( .libPaths(), lib ))
#library(Cairo) 
library(ggplot2)
library(plyr)
require(gridExtra)
require(scales)
library(RColorBrewer)
library(reshape2)
#install.packages("tidyr")
library(tidyr)
library(dplyr)
library(splitstackshape)
library(stringr)
library(data.table)
#library(networkD3)
library(bedr)
library(evaluate)
library(stringdist)
library(GenomicRanges)
require(gridExtra)
library(ggExtra)
library(Hmisc)

suppressPackageStartupMessages(library("argparse"))




# create parser object
tsv = "~/Desktop/data/genomeWide/CHM1_V4/LocalAssemblies/000324F.0.78599/ont/ontToRef.bam.tsv"
#tsv = "~/Desktop/data/genomeWide/CHM13/LocalAssemblies/PlaceOnt/SDAtoOnt.bam.tsv"
tsv = "~/Desktop/data/genomeWide/Mitchell_CHM1_V2/LocalAssemblies/SDAtoSDA.bam.tsv"
tsv = "~/Desktop/data/genomeWide/Clint/LocalAssemblies/CmpToClone/SDA.to.clone.tsv"
tsv = "~/Desktop/data/genomeWide/Clint/LocalAssemblies/wCmp//SDA.to.clone.tsv"

parser <- ArgumentParser()
parser$add_argument("-n", "--ont", default=tsv, help="Input tsv file")
args <- parser$parse_args()
args



ont = fread(args$ont)
colnames(ont)
ontl <- melt(ont, measure.vars = c("perID_by_matches", "perID_by_events", "perID_by_all")) 
ontl


mycdf = ecdf(ont$perID_by_matches)
1-mycdf(99.8)
ggplot(ontl) + stat_ecdf(aes(x=value, color=variable)) + coord_cartesian(xlim=c(99,100)) +geom_vline(aes(xintercept=99.8)) + geom_histogram(aes(x=value, fill=variable), bins=1000 )


text = data.table(ontl %>% 
  group_by(variable) %>%
  summarise(mean=mean(value), median=median(value)))


ggplot() + geom_violin(data=ontl, aes(y=value, x=variable), draw_quantiles=c(.25,.50,.75)) +
  geom_text(data=text, aes(y=mean, x=variable, label=mean))+
  geom_text(data=text, aes(y=median, x=variable, label=median))+
  coord_cartesian(ylim = c(85, 100), expand=F)

ont$perID_by_all


all  = ontl[ontl$variable == "perID_by_all", ]
all

p1 = ggplot() + geom_histogram(data=all, aes(value), bins=1000) + 
  xlab("Best %ID Match") + ylab("Count") +
  coord_cartesian(xlim=c(90,100)) + theme(text = element_text(size=20))

p2 = p1 +  coord_cartesian(xlim=c(97.5,100)) 

p1
p2





#
# LHR overlap 
#
tsv = "~/Desktop/data/genomeWide/CHM13_V2/LocalAssemblies/PlaceOnt/LHR.tsv"
lhr = fread(tsv)
colnames(lhr) = c("correct", "total", "frac", "LHR", "rname", "qname", "dove", "match", "mismatch", "insertion","deletion", "indelsOverThree", "largestIndel", "perID", "rlen", "qlen")
lhr$NumOfPSV <- cut2(lhr$total, g=10)
lhr$logdove <- cut2((lhr$dove+1), g=2)

#
# read in true values
#
bed = "~/Desktop/data/genomeWide/CHM13_V2/LocalAssemblies/PlaceOnt/intersect.bed"
bed = fread(bed)
tpairs = paste(bed$V4, bed$V10, sep="__")
dpairs = paste(lhr$qname, lhr$rname, sep = "__")
Y = as.factor(dpairs %in% tpairs)

lhr$real = Y

mmax = 300
lhr$LHR[lhr$LHR > mmax ] = mmax
lhr$LHR[lhr$LHR < -mmax ] = -mmax
p3 = ggplot() + geom_histogram(data=lhr, aes(LHR), bins=mmax) + 
  geom_vline(xintercept = 5, color="darkred", linetype=2) +
  xlab("Log Likelihood Ratio") + ylab("Count") + theme(text = element_text(size=20))
p4 = ggplot() + geom_histogram(data=lhr, aes(frac), bins=mmax) + 
  geom_vline(xintercept = 0.65, color="darkred", linetype=2) +
  xlab("Fraction of PSVs in Alignment") + ylab("Count") + theme(text = element_text(size=20))
grid.arrange(p4, p3, nrow=2)

p5 = ggplot() + geom_point(data=lhr, aes(x=LHR, y=frac, color = Y, shape = Y)) +
  geom_vline(aes(xintercept = 0), color = "red", linetype = 2, size = 1) + theme(legend.position = "bottom") +
  scale_color_manual(values=c("#ff0000","#000000")) 
p5 = ggMarginal(p5, type="histogram", bins = 200, groupFill=T, alpha = 1); p5



bed = "~/Desktop/data/genomeWide/CHM13_V2/LocalAssemblies/PlaceOnt/intersect.bed"
bed = fread(bed)
tpairs = paste(bed$V4, bed$V10, sep="__")
dpairs = paste(lhr$qname, lhr$rname, sep = "__")
Y = as.factor(dpairs %in% tpairs)

X = lhr[,c("correct", "total", "frac", "LHR", "dove", "match", "mismatch", "insertion","deletion", "indelsOverThree", "largestIndel", "perID", "rlen", "qlen")]
X

require(randomForest)
set.seed(101)
model.rf = randomForest(Y~., data=X, trees=100000, importance=TRUE)
print(model.rf)
importance(model.rf)

#install.packages("e1071")
library(e1071)
model.svm = svm(Y~., data=X)
pred <- predict(model.svm, X)
table(pred, Y)
