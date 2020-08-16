## Monte Carlo calculation of posterior probability ##

jobid <- as.numeric(commandArgs(trailingOnly=TRUE))
if(length(jobid)==0){jobid <- 1}

pdff <- function(filename){pdf(file=paste0(filename,'.pdf'),paper='a4r',height=11.7,width=16.5)}
#library('ggplot2')
library('RColorBrewer')
#library('cowplot')
#library('png')
#library('plot3D')
library('foreach')
library('doFuture')
registerDoFuture()
plan(multiprocess, workers=6)
library('doRNG')
library('LaplacesDemon')
library('RNetCDF')
#library('Rmpfr')
options(bitmapType='cairo')
mypurpleblue <- '#4477AA'
myblue <- '#66CCEE'
mygreen <- '#228833'
myyellow <- '#CCBB44'
myred <- '#EE6677'
myredpurple <- '#AA3377'
mygrey <- '#BBBBBB'
mycolours <- c(myblue, myred, mygreen, myyellow, myredpurple, mypurpleblue, mygrey, 'black')
palette(mycolours)
barpalette <- colorRampPalette(c(mypurpleblue,'white',myredpurple),space='Lab')
barpalettepos <- colorRampPalette(c('white','black'),space='Lab')
#dev.off()

set.seed(181225+jobid)

scriptname <- 'mc6'

#################################################################
######################## EQUATIONS SETUP ########################
#################################################################

generalname <- paste0('result_f_dirich')
datafilename <- 'data1'

dataabsfreqs <- as.matrix(read.csv(paste0(datafilename,'.csv'),header=FALSE,sep=','))
## format: 1 row = 1 stimulus, 1 column = 1 response

nmcsamples <- 10000

nstimuli <- nrow(dataabsfreqs)
nresponses <- ncol(dataabsfreqs)
datasize <- rowSums(dataabsfreqs)
datarelfreqs <- dataabsfreqs/datasize

## Mutual info: argument is matrix of *relative* conditional frequencies:
## 1 row = 1 stimulus, 1 column = 1 response
## output is in base=number of stimuli (2 stimuli -> bits)
MI <- function(condfreqs){
    nstim <- nrow(condfreqs)
    sum(condfreqs*log(t(t(condfreqs)/colMeans(condfreqs)), base=nstim), na.rm=TRUE)/nstim
}

dataMI <- MI(datarelfreqs)
## Dirichlet parameters for F-uniform prior
dirichpriorfreqs <- foreach(i=1:nstimuli, .combine=rbind)%do%{rep(1/nresponses, nresponses)}
dirichpriorsize <- foreach(i=1:nstimuli, .combine=c)%do%{nresponses}*0.1

dirichpostparams <- dirichpriorsize * dirichpriorfreqs + dataabsfreqs

Fsamples <- aperm(simplify2array(foreach(i=1:nstimuli)%do%{
    rdirichlet(n=nmcsamples, alpha=dirichpostparams[i,]) }) , perm=c(1,3,2))
## format: 1 row = 1 MCsample, 1 column = 1 stimulus, 1 3rdDim = 1 response

MIsamples <- foreach(i=1:nmcsamples, .combine=c)%do%{MI(Fsamples[i,,])}


histcells <- seq(from=0, to=1, length.out=round(25/diff(range(MIsamples))))
pdff(paste0(generalname,dirichpriorsize[1]))
hist(MIsamples, breaks=histcells, freq=FALSE, xlab='long-run MI/bit',
     main=paste0('Dirichlet prior with uniform reference distribution and prior size = ', dirichpriorsize[1], '\nsample MI = ', signif(dataMI,3), ' bit'))
matpoints(x=dataMI, y=0, type='p', pch=17, cex=2, col=myred)
dev.off()

#####
## Forcasts using only the mutual info of the sample, rather than the sample frequencies

nmcsamplesbis <- 40000

## Dirichlet parameters for F-uniform prior
dirichpriorfreqs <- foreach(i=1:nstimuli, .combine=rbind)%do%{rep(1/nresponses, nresponses)}
dirichpriorsize <- foreach(i=1:nstimuli, .combine=c)%do%{nresponses}*10

dirichpriorparams <- dirichpriorsize * dirichpriorfreqs

Fpriorsamples <- aperm(simplify2array(foreach(i=1:nstimuli)%do%{
    rdirichlet(n=nmcsamplesbis, alpha=dirichpriorparams[i,]) }) , perm=c(1,3,2))
## format: 1 row = 1 MCsample, 1 column = 1 stimulus, 1 3rdDim = 1 response

doubleMIsamples <- foreach(i=1:nmcsamplesbis, .combine=rbind)%dorng%{
    c(MI(foreach(s=1:nstimuli, .combine=rbind)%do%{
        tabulate(sample(x=1:nresponses, size=datasize[s], replace=TRUE,
                        Fpriorsamples[i,s,]), nbins=nresponses)/datasize[s]
    })
    , MI(Fpriorsamples[i,,]))}


miwidth <- 0.01 # width around the mutual info
MIsamplesbis <- doubleMIsamples[(doubleMIsamples[,1] > dataMI - miwidth) & (doubleMIsamples[,1] < dataMI + miwidth),2]

histcells <- seq(from=0, to=1, length.out=round(20/diff(range(MIsamplesbis))))
pdff(paste0('onlyMI_',generalname,dirichpriorsize[1]))
hist(MIsamplesbis, breaks=histcells, freq=FALSE, xlab='long-run MI/bit',
     main=paste0('Dirichlet prior with uniform reference distribution and prior size = ', dirichpriorsize[1], '\nsample MI = ', signif(dataMI,3), ' bit'))
matpoints(x=dataMI, y=0, type='p', pch=17, cex=2, col=myred)
dev.off()


mcdataMI <- doubleMIsamples[,1]
histcells <- seq(from=0, to=1, length.out=round(25/diff(range(mcdataMI))))
pdff(paste0('test',dirichpriorsize[1]))
hist(mcdataMI, breaks=histcells, freq=FALSE, xlab='long-run MI/bit',
     main=paste0('Dirichlet prior with uniform reference distribution and prior size = ', dirichpriorsize[1], '\nsample MI = ', signif(dataMI,3), ' bit'))
matpoints(x=dataMI, y=0, type='p', pch=17, cex=2, col=myred)
dev.off()


plan(sequential)
