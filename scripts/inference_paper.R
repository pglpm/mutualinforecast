## Author: Battistin, Gonzalo Cogno, Porta Mana
## Last-Updated: 2021-03-20T10:07:17+0100
################
## Script for:
## - outputting samples of prior & posterior distributions
## - calculating posteriors
## Uses Dirichlet prior
################

#### Custom setup ####
## Colour-blind friendly palettes, from https://personal.sron.nl/~pault/
## (consider using khroma package instead)
library('RColorBrewer')
mypurpleblue <- '#4477AA'
myblue <- '#66CCEE'
mygreen <- '#228833'
myyellow <- '#CCBB44'
myred <- '#EE6677'
myredpurple <- '#AA3377'
mygrey <- '#BBBBBB'
mypalette <- c(myblue, myred, mygreen, myyellow, myredpurple, mypurpleblue, mygrey, 'black')
palette(mypalette)
barpalette <- colorRampPalette(c(mypurpleblue,'white',myredpurple),space='Lab')
barpalettepos <- colorRampPalette(c('white','black'),space='Lab')
#dev.off()
####
library('data.table')
library('khroma')
library('ggplot2')
library('ggthemes')
theme_set(theme_bw(base_size=18))
scale_colour_discrete <- scale_colour_bright
#library('cowplot')
library('png')
library('foreach')
library('doRNG')
library('LaplacesDemon') # used for Dirichlet generator
library('ash')
options(bitmapType='cairo')
pdff <- function(filename){pdf(file=paste0(filename,'.pdf'),paper='a4r',height=11.7,width=16.5)} # to output in pdf format
pngf <- function(filename,res=300){png(file=paste0(filename,'.png'),height=11.7*1.2,width=16.5,units='in',res=res,pointsize=36)} # to output in pdf format
#### End custom setup ####

#######################################
#### FUNCTION TO CALCULATE MUTUAL INFO FROM FREQUENCY PAIRS
## freqs[,S] = response freqs for stimulus S: one column per stimulus
## assumes all stimuli equally probable
mutualinfo <- function(freqs,base=2){##in bits by default
    stimulusFreqs <- 1/ncol(freqs)
    jointFreqs <- freqs * stimulusFreqs
    responseFreqs <- rowSums(jointFreqs)
    sum(jointFreqs *
        log2(jointFreqs/outer(responseFreqs,rep(stimulusFreqs,ncol(freqs)))),
        na.rm=TRUE)/log2(base)
}
## function to normalize absolute frequencies
normalize <- function(freqs){freqs/sum(freqs)}

fullDataFile  <- 'SpikeCounts_and_direction.csv'
sampleIndexFile  <- 'index_mat_160.csv'
chunk <- 1
maxSpikes <- 15
priorMeanSpikes <- 0.2 # 5Hz * (40Hz/1000s)
priorWeight <- 20
priorBaseDistr <- normalize(foreach(i=0:maxSpikes, .combine=c)%do%{dgeom(x=i, prob=1/(priorMeanSpikes+1), log=FALSE)})
nPostSamples <- 1000
nPlotSamples <- 200
maxX <- 8
set.seed(149)
##
## load full recording
fullData  <- as.data.table(t(read.csv(fullDataFile,header=FALSE,sep=',')))
colnames(fullData) <- c('nspikes', 'stimulus')
stimuli <- unlist(unique(fullData[,2]))
nStimuli <- length(stimuli)
## frequencies of full recording
fullFreqs <- foreach(stim=stimuli, .combine=cbind)%do%{
    tabulate(unlist(subset(fullData, stimulus==stim)$nspikes)+1, nbins=maxSpikes+1)
}
fullMI <- mutualinfo(normalize(fullFreqs))
## frequencies of sample
chunkIndices <- as.matrix(read.csv(sampleIndexFile,header=FALSE,sep=','))[chunk,]
sampleData <- fullData[chunkIndices,]
sampleFreqs <- foreach(stim=stimuli, .combine=cbind)%do%{
    tabulate(unlist(subset(sampleData, stimulus==stim)$nspikes)+1, nbins=maxSpikes+1)
}
nSamples <- sum(sampleFreqs)
sampleMI <- mutualinfo(normalize(sampleFreqs))
##
##
## posterior superdistribution samples
postAlphas <- sampleFreqs + priorWeight*priorBaseDistr
postSamplePairs <- sapply(1:nStimuli, function(x){
    t(rdirichlet(n=nPostSamples, alpha=postAlphas[,x]))},
    simplify='array')
## posterior MI samples
    postMISamples <- foreach(i=1:nPostSamples, .combine=c)%do%{
        mutualinfo(postSamplePairs[,i,])
    }
    postMIDistr <- hist(postMISamples, breaks=seq(0,1,by=0.02), plot=FALSE)
    postMIQuantiles <- quantile(x=postMISamples, probs=c(0.025,0.5,0.975))
##
##
## prior superdistribution samples
priorSamplePairs <- sapply(1:nStimuli, function(x){
    t(rdirichlet(n=nPostSamples, alpha=priorWeight*priorBaseDistr))},
    simplify='array')
## prior MI samples
    priorMISamples <- foreach(i=1:nPostSamples, .combine=c)%do%{
        mutualinfo(priorSamplePairs[,i,])
    }
    priorMIDistr <- hist(priorMISamples, breaks=seq(0,1,by=0.02), plot=FALSE)
    priorMIQuantiles <- quantile(x=priorMISamples, probs=c(0.025,0.5,0.975))
##
##
##
## Plots
##
## plot prior samples
pdff(paste0('plots_samples',nSamples,'_chunk',chunk,'_weight',priorWeight))
##
matplot(x=0:maxSpikes, y=priorSamplePairs[,1:nPlotSamples,1],
        type='l', lty=1, lwd=2, col=paste0(mypurpleblue,'22'), ylim=c(-1,1),  xlim=c(0,maxX),
        xlab='spikes/bin', ylab='freq', cex.lab=2, cex.axis=2)
for(i in 2:nStimuli){
matplot(x=0:maxSpikes, y=-priorSamplePairs[,1:nPlotSamples,i],
        type='l', lty=1, lwd=1, col=paste0(myredpurple,'22'),
        add=TRUE)
}
matplot(x=0:maxSpikes, y=normalize(sampleFreqs[,1]),
        type='l', lty=3, lwd=4, col='black',
        add=TRUE)
matplot(x=0:maxSpikes, y=-normalize(sampleFreqs[,2]),
        type='l', lty=3, lwd=5, col='black',
        add=TRUE)
title(paste0('(',nSamples,' data samples,',
             ' chunk ', chunk,
             ', prior weight = ', priorWeight, ')',
             '\nprior superdistr'), cex.main=2)
legend('topright',c('sample freqs'),lty=c(3,2),lwd=c(4,3),cex=1.5)
## matplot(x=0:maxSpikes, y=normalize(fullFreqs[,1]),
##         type='l', lty=2, lwd=4, col='black',
##         xlab='spikes', ylab='freq',add=TRUE)
## matplot(x=0:maxSpikes, y=-normalize(fullFreqs[,2]),
##         type='l', lty=2, lwd=4, col='black',
##         xlab='spikes', ylab='freq',add=TRUE)
##
##
matplot(x=0:maxSpikes, y=postSamplePairs[,1:nPlotSamples,1],
        type='l', lty=1, lwd=2, col=paste0(mypurpleblue,'22'), ylim=c(-1,1),  xlim=c(0,maxX),
        xlab='spikes/bin', ylab='freq', cex.lab=2, cex.axis=2)
for(i in 2:nStimuli){
matplot(x=0:maxSpikes, y=-postSamplePairs[,1:nPlotSamples,i],
        type='l', lty=1, lwd=1, col=paste0(myredpurple,'22'),
        add=TRUE)
}
matplot(x=0:maxSpikes, y=normalize(sampleFreqs[,1]),
        type='l', lty=3, lwd=4, col='black',
        add=TRUE)
matplot(x=0:maxSpikes, y=-normalize(sampleFreqs[,2]),
        type='l', lty=3, lwd=5, col='black',
        add=TRUE)
matplot(x=0:maxSpikes, y=normalize(fullFreqs[,1]),
        type='l', lty=2, lwd=4, col='black',
        add=TRUE)
matplot(x=0:maxSpikes, y=-normalize(fullFreqs[,2]),
        type='l', lty=2, lwd=4, col='black',
        add=TRUE)
legend('topright',c('sample freqs','long-run freqs'),lty=c(3,2),lwd=c(4,3),cex=1.5)
title('posterior superdistr', cex.main=2)
##
##
## Prior MI
matplot(x=priorMIDistr$mids, y=priorMIDistr$density,
        type='h', lty=1, lwd=15, col=paste0(mypurpleblue,'88'), xlim=c(0,1),
        xlab='MI/bit', ylab='prob dens', cex.lab=2, cex.axis=2)
for(q in priorMIQuantiles){
    matlines(x=rep(q,2),y=c(-1,1/2)*max(priorMIDistr$density), lty=2, lwd=6, col=mygreen)
    }
    matlines(x=rep(sampleMI,2),y=c(-1,2/3)*max(priorMIDistr$density), lty=1, lwd=6, col=myredpurple)
title('prior MI distr', cex.main=2)
legend('topright',c('sample MI'),lty=1,col=myredpurple,lwd=4,cex=1.5)
##
## Posterior MI
matplot(x=postMIDistr$mids, y=postMIDistr$density,
        type='h', lty=1, lwd=15, col=paste0(mypurpleblue,'88'), xlim=c(0,1),
        xlab='MI/bit', ylab='prob dens', cex.lab=2, cex.axis=2)
for(q in postMIQuantiles){
    matlines(x=rep(q,2),y=c(-1,1/2)*max(postMIDistr$density), lty=2, lwd=6, col=mygreen)
    }
    matlines(x=rep(sampleMI,2),y=c(-1,2/3)*max(postMIDistr$density), lty=1, lwd=6, col=myredpurple)
    matlines(x=rep(fullMI,2),y=c(-1,2/3)*max(postMIDistr$density), lty=1, lwd=6, col=myyellow)
title('posterior MI distr', cex.main=2)
legend('topright',c('sample MI', 'long-run MI'),lty=1,col=c(myredpurple,myyellow),lwd=4,cex=1.5)
dev.off()


exit <- function() { invokeRestart("abort") }
exit()
