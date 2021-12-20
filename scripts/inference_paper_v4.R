## Author: Battistin, Gonzalo Cogno, Porta Mana
## Last-Updated: 2021-12-19T17:09:30+0100
################
## Script for:
## - outputting samples of prior & posterior distributions
## - calculating posteriors
## Uses Dirichlet prior
################

if(file.exists("/cluster/home/pglpm/R")){
    .libPaths(c("/cluster/home/pglpm/R",.libPaths()))
}
#### Custom setup ####
## Colour-blind friendly palettes, from https://personal.sron.nl/~pault/
## library('khroma')
## palette(colour('bright')())
## scale_colour_discrete <- scale_colour_bright
## palette(colour('muted')())
library('data.table')
## library('ggplot2')
## library('ggthemes')
## theme_set(theme_bw(base_size=18))
#library('cowplot')
library('png')
library('foreach')
library('doFuture')
library('doRNG')
registerDoFuture()
print('availableCores:')
print(availableCores())
print('availableCores-multicore:')
print(availableCores('multicore'))
if(file.exists("/cluster/home/pglpm/R")){
    plan(multicore, workers=availableCores()-1)
}else{
    plan(multisession, workers=6)
}
##library('ash')
## library('LaplacesDemon')
## library('extraDistr')
## library('mvtnorm')
## options(bitmapType='cairo')
## pdff <- function(filename){pdf(file=paste0(filename,'.pdf'),paper='a4r',height=11.7,width=16.5)} # to output in pdf format
## pngf <- function(filename,res=300){png(file=paste0(filename,'.png'),height=11.7*1.2,width=16.5,units='in',res=res,pointsize=36)} # to output in pdf format
library('nimble')
#### End custom setup ####

#######################################
#### FUNCTION TO CALCULATE MUTUAL INFO FROM JOINT DISTRIBUTION
## freqs[S,B] = freq spike count B and stimulus S (one ROW per stimulus)
## The function calculates the conditional frequencies of B|S
## and constructs a new joint distribution with equal marginals for S
## Note: don't need to normalize input to mutualinfo
mutualinfo <- function(jointFreqs,base=2){##in bits by default
    stimulusFreqs <- 1/nrow(jointFreqs)
    ## (conditional freqs B|S) * (new freq S)
    jointFreqs <- (jointFreqs/rowSums(jointFreqs)) * stimulusFreqs
    sum(jointFreqs *
        log2(jointFreqs/outer(rowSums(jointFreqs), colSums(jointFreqs))),
        na.rm=TRUE)/log2(base)
}
## function to normalize absolute frequencies
normalize <- function(freqs){freqs/sum(freqs)}
normalizerows <- function(freqs){freqs/rowSums(freqs)}
normalizecols <- function(freqs){t(t(freqs)/colSums(freqs))}

longrunDataFile  <- 'SpikeCounts_and_direction.csv'
sampleIndexFile  <- 'index_mat_160.csv'
plan(sequential)
maxSpikes <- 12
maxSpikes1 <- maxSpikes + 1
T <- 100 ## prior weight
priorMeanSpikes <- 0.4 # 0.2 = 5Hz * (40Hz/1000s)
#priorBaseDistr <- normalize(dpois(x=0:maxSpikes, lambda=priorMeanSpikes, log=FALSE))
priorBaseDistr <- normalize(dgeom(x=0:maxSpikes, prob=1/(priorMeanSpikes+1), log=FALSE))
## priorBaseDistr <- normalize(rep(1,maxSpikes1))
nDraws <- 2000
nPlotSamples <- 100
maxX <- 8
##
## load full recording
longrunData  <- as.data.table(t(read.csv(longrunDataFile,header=FALSE,sep=',')))
colnames(longrunData) <- c('nspikes', 'stimulus')
stimuli <- unique(longrunData[,stimulus])
nStimuli <- length(stimuli)
## longrunData <- longrunData[nspikes<=maxSpikes] # for debug, REMEMBER TO REMOVE
## frequencies of full recording
longrunFreqs <- foreach(stim=stimuli, .combine=rbind)%do%{
    tabulate(longrunData[stimulus==stim,nspikes]+1, nbins=maxSpikes1)
}
##priorBaseDistr <- normalize(c(colSums(longrunFreqs))+0.1)
rownames(longrunFreqs) <- nameStimulus <- paste0('stimulus',stimuli)
colnames(longrunFreqs) <- nameNspikes <- paste0('nspikes',0:maxSpikes)
longrunMI <- c(bit=mutualinfo(longrunFreqs))
## frequencies of sample
gc()
plan(sequential)
plan(multisession, workers = 6L)
nores <- foreach(chunk=0:2, .inorder=F, .packages=c('data.table','extraDistr'))%dorng%{
    pflag <- 0
    if(chunk==0){chunk <- 1
        pflag <- 1}
    chunkIndices <- as.matrix(read.csv(sampleIndexFile,header=FALSE,sep=','))[chunk,]
    sampleData <- longrunData[chunkIndices,]
    ##print(str(sampleData))
    sampleFreqs <- foreach(stim=stimuli, .combine=rbind)%do%{
        tabulate(sampleData[stimulus==stim,nspikes]+1, nbins=maxSpikes1)
    } 
    dimnames(sampleFreqs) <- dimnames(longrunFreqs)
    nSamples <- sum(sampleFreqs)
    if(pflag==0){
        sampleMI <- c(bit=mutualinfo(sampleFreqs))
    } else {
        sampleMI <- c(bit=-2)
        chunk <- 0}
    ##
    ##
    ## Alphas for Dirichlet distribution
    ##
    priorAlphas <- normalize(matrix(rep(priorBaseDistr, each=nstimuli),nrow=nstimuli))
    dimnames(priorAlphas) <- list(nameStimulus, nameNspikes)
    dAlphas <- T*priorAlphas + if(pflag==0){sampleFreqs}else{0}
    ## Generate samples
    set.seed(147+chunk)
    mcmcrun <- rdirichlet(n=nDraws, alpha=c(dAlphas))
    dim(mcmcrun) <- c(nDraws,nstimuli,maxSpikes1)
    dimnames(mcmcrun) <- list(NULL, nameStimulus, nameNspikes)
    postMISamples <- apply(mcmcrun,1,mutualinfo)
    postMIDistr <- hist(postMISamples, breaks=seq(0,1,by=0.02), plot=F)
    postMIQuantiles <- quantile(x=postMISamples, probs=c(0.025,0.5,0.975))
    condfreqSamples <- t(apply(mcmcrun,1,normalizerows))
    dim(condfreqSamples) <- c(nDraws, nStimuli, maxSpikes1)
    dimnames(condfreqSamples) <- dimnames(mcmcrun)
    ##
    ##
    ## Plot draws
    pdff(paste0('DDplots_samples',nSamples,'_chunk',chunk))
    ##
    matplot(x=0:maxSpikes, y=t(condfreqSamples[round(seq(1,nDraws,length.out=nPlotSamples)),1,]),
            type='l', lty=1, lwd=2, col=paste0(mypurpleblue,'22'), ylim=c(-1,1),  xlim=c(0,maxX),
            xlab='spikes/bin', ylab='freq', cex.lab=2, cex.axis=2)
    for(i in 2:nStimuli){
        matplot(x=0:maxSpikes, y=-t(condfreqSamples[round(seq(1,nDraws,length.out=nPlotSamples)),i,]),
                type='l', lty=1, lwd=1, col=paste0(myredpurple,'22'),
                add=TRUE)
    }
    matplot(x=0:maxSpikes, y=normalize(longrunFreqs[1,]),
            type='l', lty=1, lwd=2, col='black',
            add=TRUE)
    matplot(x=0:maxSpikes, y=-normalize(longrunFreqs[2,]),
            type='l', lty=1, lwd=2, col='black',
            add=TRUE)
    if(pflag==0){matplot(x=0:maxSpikes, y=normalize(sampleFreqs[1,]),
                         type='l', lty=2, lwd=5, col=myyellow,
                         add=TRUE)
                         matplot(x=0:maxSpikes, y=-normalize(sampleFreqs[2,]),
                                 type='l', lty=2, lwd=5, col=myyellow,
                                 add=TRUE)}
    title(paste0('(',nSamples,' data samples,',
                 ' chunk ', chunk,
                 ', prior weight = ', T, ')',
                 '\n superdistr ',chunk), cex.main=2)
    legend('topright',c('long-run freqs','sample freqs'),lty=c(1,2),lwd=c(2,5),col=c('black',myyellow),cex=1.5)
    ## Posterior MI
    matplot(x=postMIDistr$mids, y=postMIDistr$density,
            type='h', lty=1, lwd=15, col=paste0(mypurpleblue,'88'), xlim=c(0,1),
            xlab='MI/bit', ylab='prob dens', cex.lab=2, cex.axis=2)
    for(q in postMIQuantiles){
        matlines(x=rep(q,2),y=c(-1,1/2)*max(postMIDistr$density), lty=2, lwd=6, col=mygreen)
    }
    matlines(x=rep(sampleMI,2),y=c(-1,2/3)*max(postMIDistr$density), lty=4, lwd=6, col=myyellow)
    matlines(x=rep(longrunMI,2),y=c(-1,2/3)*max(postMIDistr$density), lty=1, lwd=6, col=myredpurple)
    title('posterior MI distr', cex.main=2)
    legend('topright',c('sample MI', 'long-run MI'),lty=1,col=c(myyellow,myredpurple),lwd=4,cex=1.5)
    dev.off()
    ##
    NULL
}
