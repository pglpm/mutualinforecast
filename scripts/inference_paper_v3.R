## Author: Battistin, Gonzalo Cogno, Porta Mana
## Last-Updated: 2021-07-26T15:42:08+0200
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
library('doFuture')
registerDoFuture()
library('doRNG')
library('ash')
library('LaplacesDemon') # used for Dirichlet generator
library('extraDistr')
options(bitmapType='cairo')
pdff <- function(filename){pdf(file=paste0(filename,'.pdf'),paper='a4r',height=11.7,width=16.5)} # to output in pdf format
pngf <- function(filename,res=300){png(file=paste0(filename,'.png'),height=11.7*1.2,width=16.5,units='in',res=res,pointsize=36)} # to output in png format
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
normalizerows <- function(freqs){freqs/rowSums(freqs)}
normalizecols <- function(freqs){t(t(freqs)/colSums(freqs))}

longrunDataFile  <- 'SpikeCounts_and_direction.csv'
sampleIndexFile  <- 'index_mat_160.csv'
plan(sequential)
maxSpikes <- 12
maxSpikes1 <- maxSpikes + 1
priorMeanSpikes <- 0.2 # 5Hz * (40Hz/1000s)
priorWeight <- 20
priorBaseDistr <- normalize(foreach(i=0:maxSpikes, .combine=c)%do%{dgeom(x=i, prob=1/(priorMeanSpikes+1), log=FALSE)})
nPostSamples <- 1000
nPlotSamples <- 100
maxX <- 8
set.seed(149)
##
## load full recording
longrunData  <- as.data.table(t(read.csv(longrunDataFile,header=FALSE,sep=',')))
colnames(longrunData) <- c('nspikes', 'stimulus')
stimuli <- unique(longrunData[,stimulus])
nStimuli <- length(stimuli)
## longrunData <- longrunData[nspikes<=maxSpikes] # for debug, REMEMBER TO REMOVE
## frequencies of full recording
longrunFreqs <- foreach(stim=stimuli, .combine=cbind)%do%{
    tabulate(longrunData[stimulus==stim,nspikes]+1, nbins=maxSpikes1)
}
colnames(longrunFreqs) <- nameStimulus <- paste0('stimulus',stimuli)
rownames(longrunFreqs) <- nameNspikes <- paste0('nspikes',0:maxSpikes)
longrunMI <- mutualinfo(normalize(longrunFreqs))
names(longrunMI) <- 'bit'
## frequencies of sample
gc()
plan(sequential)
plan(multisession, workers = 6L)
nores <- foreach(chunk=0:2)%dorng%{
    pflag <- 0
    if(chunk==0){chunk <- 1
        pflag <- 1}
chunkIndices <- as.matrix(read.csv(sampleIndexFile,header=FALSE,sep=','))[chunk,]
sampleData <- longrunData[chunkIndices,]
sampleFreqs <- foreach(stim=stimuli, .combine=cbind)%do%{
    tabulate(samplrownames(longrunFreqs)eData[stimulus==stim,nspikes]+1, nbins=maxSpikes1)
}
dimnames(sampleFreqs) <- dimnames(longrunFreqs)
nSamples <- sum(sampleFreqs)
if(pflag==0){sampleMI <- mutualinfo(normalize(sampleFreqs))} else {sampleMI <- 2}
names(sampleMI) <- 'bit'
##
##
## posterior superdistribution samples
##

priorMeanSpikes <- 0.2 # 5Hz * (40Hz/1000s)
priorBaseDistr <- normalize(foreach(i=0:maxSpikes, .combine=c)%do%{dgeom(x=i, prob=1/(priorMeanSpikes+1), log=FALSE)})
T <- 500 ## prior weight
priorAlphas <- matrix(rep(priorBaseDistr, nstimuli),nrow=nstimuli)
dimnames(priorAlphas) <- list(nameStimulus, nameNspikes)
dAlphas <- T*priorAlphas - 1 + if(pflag==0){t(sampleFreqs)}else{0}
smoothness <- 1000
smoothm <- sqrt(smoothness)*diff(diag(maxSpikes1),differences=2)
##
outfile <- paste0('_LDoutput',chunk)
thinning <- 2 #10
mcadapt <- 1000 #10e3
mcburnin <- mcadapt
mciterations <- 2000 #mcburnin + 10e3
mcstatus <- 100 #1e3
parnames <- paste0('F',rep(0:maxSpikes,each=nStimuli),'_',rep(stimuli,times=maxSpikes1))
stepwidth <- 1 #c(rep(nn/100,length(inu)),rep(nnyR/100,length(iR)))
nsteps <- 100 #c(rep(nn/100,length(inu)),rep(nnyR/100,length(iR)))
nprev <- 0
covm <- list()
B <- list()
##
PGF <- function(data){
    c(maxSpikes1 * t(normalize(t(rdirichlet(2,(dAlphas+1)*100+1)))))
}
##
mydata <- list(y=1, PGF=PGF,
               parm.names=parnames,
               mon.names=c('lpf','MI')
               ##nn=nn,
               ##L=L,
               ##T=T,
               ##hh=hh,
               ##ft=ft,
               ##rpos=rpos,
               ##inu=inu,
               ##iR=iR,
               ##smoothm=smoothm,
               ##nnyR=***
               )
##
logprob <- function(parm,data){
    parm <- interval(parm,0,Inf)
    FF <- parm
    dim(FF) <- c(nstimuli, maxSpikes1)
    sFF <- rowSums(FF)
    tFF <- sum(sFF)
    mi <- mutualinfo(t(FF/tFF))
    ##
##
    lpf <- sum(dAlphas*log(FF/tFF), na.rm=TRUE) - sum(tcrossprod(FF/sFF, smoothm)^2)
##
    LP <- lpf - 1e5*(tFF-maxSpikes1)^2
##
    list(LP=LP, Dev=-2*LP, Monitor=c(lpf,mi), yhat=1, parm=parm)
}
##
Initial.Values <- PGF(mydata)
## print('running Monte Carlo...')
mcoutput <- LaplacesDemon(logprob, mydata, Initial.Values,
                      Covar=NULL,
                      ##Covar=covm,
                      Thinning=thinning,
                      Iterations=mciterations, Status=mcstatus,
                      Debug=list(DB.chol=FALSE, DB.eigen=FALSE,
                                 DB.MCSE=FALSE, DB.Model=FALSE),
                      LogFile=paste0('_LDlog',outfile), #Type="MPI", #Packages=c('MASS'),
                      ##Algorithm="RDMH"#, Specs=list(B=list(1:d,d1:d2,d3:dnp))
                      ##Algorithm="Slice", Specs=list(B=NULL, Bounds=c(0,1), m=Inf, Type="Continuous", w=0.001)
                      ##Algorithm="MALA", Specs=list(A=1e7, alpha.star=0.574, gamma=mcadapt, delta=1, epsilon=c(1e-6,1e-7))
                      ##Algorithm="pCN", Specs=list(beta=0.001)
                      Algorithm="AFSS", Specs=list(A=mcadapt, B=B, m=nsteps, n=nprev, w=stepwidth)
                      ##Algorithm="AIES", Specs=list(Nc=4*nparm, Z=NULL, beta=2, CPUs=1, Packages=NULL, Dyn.libs=NULL)                      
                      ##Algorithm="DRM", Specs=NULL
                      )






    
## Preparation as prior as pseudocount data
priorMeanSpikes <- 0.2 # 5Hz * (40Hz/1000s)
priorBaseDistr <- normalize(foreach(i=0:maxSpikes, .combine=c)%do%{dgeom(x=i, prob=1/(priorMeanSpikes+1), log=FALSE)})
priorWeight <- 1000
priorBaseData <- foreach(count=0:maxSpikes, .combine=c)%do%{
    rep(count, times=max(1,round(priorBaseDistr[count+1] * priorWeight)))
}
## priorData <- data.table(nspikes=rep(priorBaseData, times=nStimuli),
##                         stimulus=rep(stimuli, each=length(priorBaseData)))
    priorData <- rbind(
        data.table(nspikes=priorBaseData,
                        stimulus=rep(NA, times=length(priorBaseData))),
        data.table(nspikes=rep(NA,times=nStimuli),
                   stimulus=stimuli)
    )
##     datum <- priorData[,nspikes]
## datum <- datum[!is.na(datum)]
## for(level in setdiff(0:maxSpikes, as.numeric(names(table(datum))))){
##     datamcr <- rbind(priorData, data.table(nspikes=level), fill=TRUE)
## }
    if(pflag==0){datamcr <- rbind(priorData, sampleData)
    } else {
        datamcr <- priorData
    chunk <- 0}
## include all possible spike counts up to maxSpikes
##
## hyperparameters
hyper <- setHyperparams(aPhi=c(1/13,1/20))
##
covNames <- colnames(sampleData)
outfile <- paste0('_mcoutput',chunk)
mcmcrun <- profRegr(excludeY=TRUE, xModel='Discrete', nSweeps=160e3, nBurn=20e3, nFilter=80, data=as.data.frame(datamcr), nClusInit=80, covNames=covNames, discreteCovs=covNames, nProgress=5e3, seed=148, output=outfile, useHyperpriorR1=FALSE, useNormInvWishPrior=TRUE, alpha=-2, hyper=hyper)
##
## Save MCMC samples
## log-likelihood and log-posteriors
fc <- file(paste0(outfile,"_logPost.txt"))
logPost <- sapply(strsplit(readLines(fc), " +"), as.numeric)
rownames(logPost) <- c('log-post','log-likelihood','log-prior')
close(fc)
## Samples of numbers of clusters
fc <- file(paste0(outfile,'_nClusters.txt'))
nList <- sapply(strsplit(readLines(fc), " +"), as.integer)
close(fc)
## Samples of Dirichlet-process alpha
fc <- file(paste0(outfile,'_alpha.txt'))
alphaList <- sapply(strsplit(readLines(fc), " +"), as.numeric)
close(fc)
## Samples of cluster weights
fc <- file(paste0(outfile,'_psi.txt'))
psiList <- lapply(strsplit(readLines(fc), " +"), function(x){x <- as.numeric(x); x[x>=0]})
close(fc)
##  Samples of Dirichlet-distribution phis (discrete covariates)
fc <- file(paste0(outfile,'_phi.txt'))
phiList <- lapply(strsplit(readLines(fc), " +"), function(x){x <- as.numeric(x); x[x>=0]})
close(fc)
nCat <- mcmcrun$nCategories
cumcats <- c(0,cumsum(nCat))
for(i in 1:length(nList)){
    catgroups <- cumcats*nList[i]
    datum <- phiList[[i]]
    phiList[[i]] <- lapply(1:length(mcmcrun$nCategories), function(j){
        y <- datum[(catgroups[j]+1):catgroups[j+1]]
        dim(y) <- c(nList[i], mcmcrun$nCategories[j])
        rownames(y) <- paste0('cluster',1:nList[i])
        aperm(y)
    })
    names(phiList[[i]]) <- mcmcrun$covNames
}
        if(mcmcrun$alpha>0){rgalpha <- c(0, mcmcrun$alpha+1)
    alphaList <- rep(mcmcrun$alpha,length(nList))}
## Save the samples above
MCMCdata <- list(nList=nList, alphaList=alphaList, psiList=psiList, phiList=phiList, logPost=logPost)
##
save.image(file=paste0('_MCdata_chunk',chunk,'.RData'))
##
##
## Diagnostic plots
pdff(paste0('mcsummary',chunk))
matplot(MCMCdata$logPost[2,], type='l',ylim=range(MCMCdata$logPost[2,],na.rm=T,finite=T),ylab='log-likelihood')
matplot(MCMCdata$nList,type='l',ylim=range(MCMCdata$nList,na.rm=T,finite=T),ylab='no. clusters')
matplot(MCMCdata$alphaList,type='l',ylim=range(MCMCdata$alphaList,na.rm=T,finite=T),ylab='alpha')
for(i in c(1,3)){matplot(MCMCdata$logPost[i,],type='l',ylim=range(MCMCdata$logPost[i,],na.rm=T,finite=T))}
dev.off()
##
##
## Samples of conditional frequencies
## dimensions: (sample, nspikes, stimulus)
condfreqSamples <- foreach(sample=seq_along(MCMCdata$nList), .combine=cbind, .inorder=FALSE)%do%{
    ## weights <- normalizerows(t(t(MCMCdata$phiList[[sample]]$stimulus) * MCMCdata$psiList[[sample]]))
    tcrossprod(
        normalizecols(MCMCdata$phiList[[sample]]$nspikes),
        normalizerows(t(t(MCMCdata$phiList[[sample]]$stimulus) * MCMCdata$psiList[[sample]]))
    )
}
dim(condfreqSamples) <- c(maxSpikes1, nStimuli, length(MCMCdata$nList))
condfreqSamples <- aperm(condfreqSamples, c(3,1,2))
dimnames(condfreqSamples) <- list(NULL, rownames(longrunFreqs), colnames(longrunFreqs))
##
postMISamples <- foreach(sample=seq_along(MCMCdata$nList), .combine=c)%do%{
    mutualinfo(condfreqSamples[sample,,])
}
postMIDistr <- hist(postMISamples, breaks=seq(0,1,by=0.02), plot=F)
postMIQuantiles <- quantile(x=postMISamples, probs=c(0.025,0.5,0.975))
##
## plot posterior samples
pdff(paste0('MCplots_samples',nSamples,'_chunk',chunk))
##
matplot(x=0:maxSpikes, y=t(condfreqSamples[round(seq(1,length(MCMCdata$nList),length.out=nPlotSamples)),,1]),
        type='l', lty=1, lwd=2, col=paste0(mypurpleblue,'22'), ylim=c(-1,1),  xlim=c(0,maxX),
        xlab='spikes/bin', ylab='freq', cex.lab=2, cex.axis=2)
for(i in 2:nStimuli){
matplot(x=0:maxSpikes, y=-t(condfreqSamples[round(seq(1,length(MCMCdata$nList),length.out=nPlotSamples)),,1]),
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
## Posterior MI
matplot(x=postMIDistr$mids, y=postMIDistr$density,
        type='h', lty=1, lwd=15, col=paste0(mypurpleblue,'88'), xlim=c(0,1),
        xlab='MI/bit', ylab='prob dens', cex.lab=2, cex.axis=2)
for(q in postMIQuantiles){
    matlines(x=rep(q,2),y=c(-1,1/2)*max(postMIDistr$density), lty=2, lwd=6, col=mygreen)
    }
    matlines(x=rep(sampleMI,2),y=c(-1,2/3)*max(postMIDistr$density), lty=1, lwd=6, col=myredpurple)
    matlines(x=rep(longrunMI,2),y=c(-1,2/3)*max(postMIDistr$density), lty=1, lwd=6, col=myyellow)
title('posterior MI distr', cex.main=2)
legend('topright',c('sample MI', 'long-run MI'),lty=1,col=c(myredpurple,myyellow),lwd=4,cex=1.5)
dev.off()
##
NULL
}












## Prior superdistribution samples
##
## preparation of dataset
nPlotSamples <- 100
outfile <- paste0('_priormcoutput')
covNames <- colnames(sampleData)
## include all possible spike counts up to maxSpikes
priorMeanSpikes <- 0.2 # 5Hz * (40Hz/1000s)
priorBaseDistr <- normalize(foreach(i=0:maxSpikes, .combine=c)%do%{dgeom(x=i, prob=1/(priorMeanSpikes+1), log=FALSE)})
priorWeight <- 100
priorBaseData <- foreach(count=0:maxSpikes, .combine=c)%do%{
    rep(count, times=round(priorBaseDistr[count+1] * priorWeight))
}
## priorData <- data.table(nspikes=rep(priorBaseData, times=nStimuli),
##                         stimulus=rep(stimuli, each=length(priorBaseData)))
priorData <- data.table(nspikes=priorBaseData,
                        stimulus=rep(NA, times=length(priorBaseData)))
datum <- priorData[,nspikes]
for(level in setdiff(0:maxSpikes, as.numeric(names(table(datum))))){
    priorData <- rbind(priorData, data.table(nspikes=level), fill=TRUE)
}
for(level in stimuli){
    priorData <- rbind(priorData, data.table(stimulus=level), fill=TRUE)
}
## hyperparameters
hyper <- setHyperparams(aPhi=c(0.1,10))
##
mcmcrun <- profRegr(excludeY=TRUE, xModel='Discrete', nSweeps=20e3, nBurn=20e3, nFilter=20, data=as.data.frame(priorData), nClusInit=80, covNames=covNames, discreteCovs=covNames, nProgress=1e3, seed=148, output=outfile, useHyperpriorR1=FALSE, useNormInvWishPrior=TRUE, alpha=4, hyper=hyper)
##
## Save MCMC samples
## log-likelihood and log-posteriors
fc <- file(paste0(outfile,"_logPost.txt"))
logPost <- sapply(strsplit(readLines(fc), " +"), as.numeric)
rownames(logPost) <- c('log-post','log-likelihood','log-prior')
close(fc)
## Samples of numbers of clusters
fc <- file(paste0(outfile,'_nClusters.txt'))
nList <- sapply(strsplit(readLines(fc), " +"), as.integer)
close(fc)
## Samples of Dirichlet-process alpha
fc <- file(paste0(outfile,'_alpha.txt'))
alphaList <- sapply(strsplit(readLines(fc), " +"), as.numeric)
close(fc)
## Samples of cluster weights
fc <- file(paste0(outfile,'_psi.txt'))
psiList <- lapply(strsplit(readLines(fc), " +"), function(x){x <- as.numeric(x); x[x>=0]})
close(fc)
##  Samples of Dirichlet-distribution phis (discrete covariates)
fc <- file(paste0(outfile,'_phi.txt'))
phiList <- lapply(strsplit(readLines(fc), " +"), function(x){x <- as.numeric(x); x[x>=0]})
close(fc)
nCat <- mcmcrun$nCategories
cumcats <- c(0,cumsum(nCat))
for(i in 1:length(nList)){
    catgroups <- cumcats*nList[i]
    datum <- phiList[[i]]
    phiList[[i]] <- lapply(1:length(mcmcrun$nCategories), function(j){
        y <- datum[(catgroups[j]+1):catgroups[j+1]]
        dim(y) <- c(nList[i], mcmcrun$nCategories[j])
        rownames(y) <- paste0('cluster',1:nList[i])
        aperm(y)
    })
    names(phiList[[i]]) <- mcmcrun$covNames
}
## Save the samples above
MCMCdata <- list(nList=nList, alphaList=alphaList, psiList=psiList, phiList=phiList, logPost=logPost)
##
save.image(file=paste0('_MCdata_prior_chunk',chunk,'.RData'))
##
##
## Diagnostic plots
pdff('mcsummary_prior')
matplot(MCMCdata$logPost[2,], type='l',ylim=range(MCMCdata$logPost[2,],na.rm=T,finite=T),ylab='log-likelihood')
matplot(MCMCdata$nList,type='l',ylim=range(MCMCdata$nList,na.rm=T,finite=T),ylab='no. clusters')
matplot(MCMCdata$alphaList,type='l',ylim=range(MCMCdata$alphaList,na.rm=T,finite=T),ylab='alpha')
for(i in c(1,3)){matplot(MCMCdata$logPost[i,],type='l',ylim=range(MCMCdata$logPost[i,],na.rm=T,finite=T))}
dev.off()
##
##
## Samples of conditional frequencies
## dimensions: (sample, nspikes, stimulus)
plan(sequential)
plan(multisession, workers = 6L)
condfreqSamples <- foreach(sample=seq_along(MCMCdata$nList), .combine=cbind, .inorder=FALSE)%dopar%{
    ## weights <- normalizerows(t(t(MCMCdata$phiList[[sample]]$stimulus) * MCMCdata$psiList[[sample]]))
    tcrossprod(
        normalizecols(MCMCdata$phiList[[sample]]$nspikes),
        normalizerows(t(t(MCMCdata$phiList[[sample]]$stimulus) * MCMCdata$psiList[[sample]]))
    )
}
dim(condfreqSamples) <- c(maxSpikes1, nStimuli, length(MCMCdata$nList))
condfreqSamples <- aperm(condfreqSamples, c(3,1,2))
dimnames(condfreqSamples) <- list(NULL, rownames(longrunFreqs), colnames(longrunFreqs))
##
plan(sequential)
plan(multisession, workers = 6L)
postMISamples <- foreach(sample=seq_along(MCMCdata$nList), .combine=c)%dopar%{
    mutualinfo(condfreqSamples[sample,,])
}
postMIDistr <- hist(postMISamples, breaks=seq(0,1,by=0.02), plot=F)
postMIQuantiles <- quantile(x=postMISamples, probs=c(0.025,0.5,0.975))
##
## plot posterior samples
pdff(paste0('MCplots_samples',nSamples,'_prior'))
##
matplot(x=0:maxSpikes, y=t(condfreqSamples[round(seq(1,length(MCMCdata$nList),length.out=nPlotSamples)),,1]),
        type='l', lty=1, lwd=2, col=paste0(mypurpleblue,'22'), ylim=c(-1,1),  xlim=c(0,maxX),
        xlab='spikes/bin', ylab='freq', cex.lab=2, cex.axis=2)
for(i in 2:nStimuli){
matplot(x=0:maxSpikes, y=-t(condfreqSamples[round(seq(1,length(MCMCdata$nList),length.out=nPlotSamples)),,1]),
        type='l', lty=1, lwd=1, col=paste0(myredpurple,'22'),
        add=TRUE)
}
for(i in 2:nStimuli){
matplot(x=0:maxSpikes, y=normalize(sampleFreqs[,i]),
        type='l', lty=3, lwd=4, col='black',
        add=TRUE)
}
title(paste0('(',nSamples,' data samples,',
             ', prior', ')',
             '\nprior superdistr'), cex.main=2)
legend('topright',c('sample freqs'),lty=c(3,2),lwd=c(4,3),cex=1.5)
## Posterior MI
matplot(x=postMIDistr$mids, y=postMIDistr$density,
        type='h', lty=1, lwd=15, col=paste0(mypurpleblue,'88'), xlim=c(0,1),
        xlab='MI/bit', ylab='prob dens', cex.lab=2, cex.axis=2)
for(q in postMIQuantiles){
    matlines(x=rep(q,2),y=c(-1,1/2)*max(postMIDistr$density), lty=2, lwd=6, col=mygreen)
    }
    matlines(x=rep(sampleMI,2),y=c(-1,2/3)*max(postMIDistr$density), lty=1, lwd=6, col=myredpurple)
    matlines(x=rep(longrunMI,2),y=c(-1,2/3)*max(postMIDistr$density), lty=1, lwd=6, col=myyellow)
title('posterior MI distr', cex.main=2)
legend('topright',c('sample MI', 'long-run MI'),lty=1,col=c(myredpurple,myyellow),lwd=4,cex=1.5)
dev.off()

















































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


## matplot(x=0:maxSpikes, y=normalize(longrunFreqs[,1]),
##         type='l', lty=2, lwd=4, col='black',
##         xlab='spikes', ylab='freq',add=TRUE)
## matplot(x=0:maxSpikes, y=-normalize(longrunFreqs[,2]),
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
matplot(x=0:maxSpikes, y=normalize(longrunFreqs[,1]),
        type='l', lty=2, lwd=4, col='black',
        add=TRUE)
matplot(x=0:maxSpikes, y=-normalize(longrunFreqs[,2]),
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
    matlines(x=rep(longrunMI,2),y=c(-1,2/3)*max(postMIDistr$density), lty=1, lwd=6, col=myyellow)
title('posterior MI distr', cex.main=2)
legend('topright',c('sample MI', 'long-run MI'),lty=1,col=c(myredpurple,myyellow),lwd=4,cex=1.5)
dev.off()


exit <- function() { invokeRestart("abort") }
exit()
