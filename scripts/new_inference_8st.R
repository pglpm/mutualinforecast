## Author: Battistin, Gonzalo Cogno, Porta Mana
## Last-Updated: 2021-12-20T07:32:39+0100
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
## library('nimble')
#### End custom setup ####

#######################################
#### FUNCTION TO CALCULATE MUTUAL INFO FROM JOINT DISTRIBUTION
## freqs[S,B] = freq spike count B and stimulus S (one ROW per stimulus)
## The function calculates the conditional frequencies of B|S
## and if requested constructs
## a new joint distribution with equal marginals for S
## Note: don't need to normalize input to mutualinfo
mutualinfo <- function(jointFreqs, equalstim=FALSE, base=2L){##in bits by default
    if(equalstim){
        jointFreqs <- (jointFreqs/rowSums(jointFreqs))/nrow(jointFreqs)
    }else{
        jointFreqs <- jointFreqs/sum(jointFreqs)
    }
    jointFreqs[is.na(jointFreqs)] <- 0
    sum(jointFreqs *
        log2(jointFreqs/outer(rowSums(jointFreqs), colSums(jointFreqs))),
        na.rm=TRUE)/log2(base)
}
entropy <- function(freqs, base=2L){
    freqs <- cbind(freqs)
    freqs <- t(t(freqs)/colSums(freqs, na.rm=T))
    colSums(freqs*log2(1/freqs), na.rm=T)/log2(base)
}
## function to normalize absolute frequencies
normalize <- function(freqs){freqs/sum(freqs, na.rm=T)}
normalizerows <- function(freqs){freqs/rowSums(freqs, na.rm=T)}
normalizecols <- function(freqs){t(t(freqs)/colSums(freqs, na.rm=T))}

equalstim <- FALSE
longrunDataFile  <- 'SpikeCounts_and_direction.csv'
#sampleIndexFile  <- 'index_mat_80.csv'
maxSpikes <- 12
maxSpikes1 <- maxSpikes + 1
T <- 10 ## prior weight
priorMeanSpikes <- 0.4 # 0.2 = 5Hz * (40Hz/1000s)
#priorBaseDistr <- normalize(dpois(x=0:maxSpikes, lambda=priorMeanSpikes, log=FALSE))
priorBaseDistr <- normalize(dgeom(x=0:maxSpikes, prob=1/(priorMeanSpikes+1), log=FALSE))
## priorBaseDistr <- normalize(rep(1,maxSpikes1))
##
## load full recording
longrunData  <- as.data.table(t(read.csv(longrunDataFile,header=FALSE,sep=',')))
colnames(longrunData) <- c('nspikes', 'stimulus')
longrunData$stimulus_lag2 <- NA
longrunData$stimulus_lag2[3:nrow(longrunData)] <- 1L*longrunData$stimulus[3:nrow(longrunData)] + 2L*longrunData$stimulus[2:(nrow(longrunData)-1)] + 4L*longrunData$stimulus[1:(nrow(longrunData)-2)]
##
stimuli <- sort(unique(longrunData[2:nrow(longrunData),stimulus_lag2],na.rm=T))
nStimuli <- length(stimuli)
## longrunData <- longrunData[nspikes<=maxSpikes] # for debug, REMEMBER TO REMOVE
## frequencies of full recording
longrunFreqs <- t(sapply(stimuli, function(stim){
    tabulate(longrunData[stimulus_lag2==stim,nspikes]+1L, nbins=maxSpikes1)
}))
##priorBaseDistr <- normalize(c(colSums(longrunFreqs))+0.1)
rownames(longrunFreqs) <- nameStimulus <- paste0('stimulus',stimuli)
colnames(longrunFreqs) <- nameNspikes <- paste0('nspikes',0:maxSpikes)
longrunMI <- c(bit=mutualinfo(longrunFreqs, equalstim=equalstim))

nDraws <- 2^14
nPlotSamples <- 100
maxX <- 8
priorAlphas <- NULL
for(i in 1:nStimuli){priorAlphas <- rbind(priorAlphas, priorBaseDistr)}
dimnames(priorAlphas) <- list(nameStimulus, nameNspikes)
T <- 32 ## prior weight
startr <- 3 ## 3 or 0 for sequence with max uniformity of stimuli
set.seed(147)
##
##pdff(paste0('MIhistogram_4stimEQ_prior',T,'_start',(if(startr==2){2}else{'OPT'})))
pdff('_justtest')
for(lsample in c(2^(5:13), nrow(longrunData)-2)){
    if(startr==0){
        entseq <- foreach(i=2L:(nrow(longrunData)-lsample+2L), .combine=c, .inorder=TRUE)%dopar%{
            entropy(
                tabulate(longrunData[i-1L+(1L:lsample),stimulus_lag2]+1L, nbins=nStimuli)
            )
        }
        entorder <- order(entseq, decreasing=T)+1L
        startbin <- entorder[1]
    }else{
        startbin <- 3
    }
    ##
    sampleData <- longrunData[startbin-1+(1:lsample),]
    sampleFreqs <- t(sapply(stimuli, function(stim){
        tabulate(sampleData[stimulus_lag2==stim,nspikes]+1L, nbins=maxSpikes1)
    }))
    dimnames(sampleFreqs) <- dimnames(longrunFreqs)
    sampleMI <- c(bit=mutualinfo(sampleFreqs, equalstim=equalstim))
    ##
    ## Alphas for Dirichlet distribution
    dAlphas <- sampleFreqs + T*priorAlphas
    ## Generate samples
    if(FALSE){
        mcmcrun <- LaplacesDemon::rdirichlet(n=nStimuli*nDraws, alpha=dAlphas)
        dim(mcmcrun) <- c(nStimuli, nDraws, maxSpikes1)
        mcmcrun <- aperm(mcmcrun, c(2,1,3))
    }else{
        mcmcrun <- LaplacesDemon::rdirichlet(n=nDraws, alpha=c(dAlphas))
        dim(mcmcrun) <- c(nDraws, nStimuli, maxSpikes1)
        }
    dimnames(mcmcrun) <- list(NULL, nameStimulus, nameNspikes)
    postMISamples <- apply(mcmcrun, 1, function(x){mutualinfo(x, equalstim=equalstim)})
    ##
    postMIDistr <- thist(postMISamples)
    ##
    tplot(x=postMIDistr$breaks,y=postMIDistr$density,
          xlim=c(0,max(postMIDistr$breaks,longrunMI,sampleMI)),
          xlab='long-run MI/bit', ylab='probability density',
          main=paste0(lsample, ' observed bins; start observation at bin ', startbin),
          cex.main=1.25)
    abline(v=longrunMI, col=2, lty=1, lwd=2)
    abline(v=sampleMI, col=3, lty=2, lwd=2)
    postMIQuantiles <- quant(x=postMISamples, probs=c(2.5,50,97.5)/100)
    abline(v=postMIQuantiles, col=c(7,7,7), lty=c(4,4,4), lwd=c(3,3,3))
    legend('topright', legend=c(names(postMIQuantiles),'sample','true'), col=c(7,7,7,3,2), lty=c(4,4,4,2,1), lwd=3, bty='n', cex=1.5)
}
dev.off()



condfreqSamples <- t(apply(mcmcrun,1,normalizerows))
dim(condfreqSamples) <- c(nDraws, nStimuli, maxSpikes1)
dimnames(condfreqSamples) <- dimnames(mcmcrun)
    ##

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


## frequencies of sample
gc()
nores <- foreach(chunk=0:2, .inorder=F, .packages=c('data.table'))%dorng%{
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
