## Monte Carlo calculation of posterior probability ##

ntasks <- as.numeric(commandArgs(trailingOnly=TRUE))

pdff <- function(filename){pdf(file=paste0(filename,'.pdf'),paper='a4r',height=11.7,width=16.5)}
#library('ggplot2')
library('RColorBrewer')
#library('cowplot')
#library('png')
#library('plot3D')
library('foreach')
library('doFuture')

registerDoFuture()
print(paste0('available workers: ', availableCores()))
if(length(ntasks)==0){ntasks <- availableCores()}
plan(multiprocess, workers=ntasks-1)
print(paste0('number of workers: ', nbrOfWorkers()))

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
dev.off()

set.seed(181225+0)

#################################################################
######################## EQUATIONS SETUP ########################
#################################################################

generalname <- paste0('result_f_dirich')
datafilename <- 'data1'

## read example sample
dataabsfreqs <- as.matrix(read.csv(paste0(datafilename,'.csv'),header=FALSE,sep=','))
## format: 1 row = 1 stimulus, 1 column = 1 response

nmcsamples <- 1e5

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

## MI of the example sample
dataMI <- MI(datarelfreqs)

## Dirichlet parameters for F-uniform prior
priorsize <- nresponses * 10

dirichpriorfreqs <- foreach(i=1:nstimuli, .combine=rbind)%do%{rep(1/nresponses, nresponses)}
dirichpriorsize <- foreach(i=1:nstimuli, .combine=c)%do%{priorsize}

dirichpostparams <- dirichpriorsize * dirichpriorfreqs + dataabsfreqs

## generate posterior samples of long-run frequency distrs
Fsamples <- aperm(simplify2array(foreach(i=1:nstimuli)%do%{
    rdirichlet(n=nmcsamples, alpha=dirichpostparams[i,]) }) , perm=c(1,3,2))
## format: 1 row = 1 MCsample, 1 column = 1 stimulus, 1 3rdDim = 1 response

## generate posterior samples of corresponding long-run MI
MIsamples <- foreach(i=1:nmcsamples, .combine=c)%dopar%{MI(Fsamples[i,,])}

## generate posterior samples of corresponding log-likelihoods
#LLLsamples <- foreach(i=1:nmcsamples, .combine=c)%:%foreach(s=1:nstimuli, .combine='+')%dopar%{dmultinom(dataabsfreqs[s,], prob=Fsamples[i,s,], log=TRUE)}
LLLsamples <- foreach(i=1:nmcsamples, .combine=c)%:%foreach(s=1:nstimuli, .combine='+')%dopar%{lfactorial(datasize[s]) + sum(dataabsfreqs[s,]*log(Fsamples[i,s,])-lfactorial(dataabsfreqs[s,]))}
    
## save posterior samples
#saveRDS(Fsamples, paste0('Fsamples_', generalname, priorsize,'.rds'))
saveRDS(MIsamples, paste0('MIsamples_', generalname, priorsize,'.rds'))
saveRDS(LLLsamples, paste0('LLLsamples_', generalname, priorsize,'.rds'))

write.table(cbind(LLLsamples, MIsamples),file= paste0('LLL-MI_samples_', generalname, priorsize,'.csv'),row.names=FALSE,col.names=FALSE, sep=',')

## find and same long-run frequencies with min & max MI
#minmi <- which(MIsamples==min(MIsamples))
#maxmi <- which(MIsamples==max(MIsamples))

## write.table(Fsamples[minmi,,],file= paste0('minMI_Fsamples_', generalname, priorsize,'.csv'),row.names=FALSE,col.names=FALSE, sep=',')
## write.table(Fsamples[maxmi,,],file= paste0('maxMI_Fsamples_', generalname, priorsize,'.csv'),row.names=FALSE,col.names=FALSE, sep=',')


exit()

## credibility interval & median
quants <- unname(quantile(MIsamples, probs=c(0.025, 0.5, 0.975)))

## histogram of posterior for long-run MI
histcells <- seq(from=0, to=1, length.out=round(75/diff(range(MIsamples))))
pdff(paste0(generalname,priorsize))
hist(MIsamples, breaks=histcells, freq=FALSE, xlab='long-run MI/bit',
     main=paste0('Dirichlet prior with uniform reference distribution and prior size = ', priorsize, '\n95% probability range: (', signif(quants[1],3), ', ', signif(quants[3],3), ') bit\nmedian: ', signif(quants[2],3),
                 ' bit\nsample MI = ', signif(dataMI,3), ' bit'))
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
pdff(paste0('onlyMI_',generalname,priorsize))
hist(MIsamplesbis, breaks=histcells, freq=FALSE, xlab='long-run MI/bit',
     main=paste0('Dirichlet prior with uniform reference distribution and prior size = ', priorsize, '\nsample MI = ', signif(dataMI,3), ' bit'))
matpoints(x=dataMI, y=0, type='p', pch=17, cex=2, col=myred)
dev.off()


mcdataMI <- doubleMIsamples[,1]
histcells <- seq(from=0, to=1, length.out=round(25/diff(range(mcdataMI))))
pdff(paste0('test',priorsize))
hist(mcdataMI, breaks=histcells, freq=FALSE, xlab='long-run MI/bit',
     main=paste0('Dirichlet prior with uniform reference distribution and prior size = ', priorsize, '\nsample MI = ', signif(dataMI,3), ' bit'))
matpoints(x=dataMI, y=0, type='p', pch=17, cex=2, col=myred)
dev.off()


plan(sequential)
