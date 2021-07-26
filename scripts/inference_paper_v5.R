## Author: Battistin, Gonzalo Cogno, Porta Mana
## Last-Updated: 2021-07-26T21:02:08+0200
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
maxSpikes2 <- 2 * maxSpikes1
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
stimulusVals <- unique(longrunData[,stimulus])
nStimuli <- length(stimulusVals)
nspikesVals <- 0:maxSpikes
## longrunData <- longrunData[nspikes<=maxSpikes] # for debug, REMEMBER TO REMOVE
## frequencies of full recording
longrunFreqs <- foreach(stim=stimulusVals, .combine=rbind)%do%{
    tabulate(longrunData[stimulus==stim,nspikes]+1, nbins=maxSpikes1)
}
##priorBaseDistr <- normalize(c(colSums(longrunFreqs))+0.1)
rownames(longrunFreqs) <- nameStimulus <- paste0('stimulus',stimulusVals)
colnames(longrunFreqs) <- nameNspikes <- paste0('nspikes',nspikesVals)
longrunMI <- c(bit=mutualinfo(longrunFreqs))
## frequencies of sample
gc()
plan(sequential)
plan(multisession, workers = 6L)
nores <- foreach(chunk=0:2, .inorder=F, .packages=c('data.table','LaplacesDemon','extraDistr'))%dorng%{
    pflag <- 0
    if(chunk==0){chunk <- 1
        pflag <- 1}
    chunkIndices <- as.matrix(read.csv(sampleIndexFile,header=FALSE,sep=','))[chunk,]
    sampleData <- longrunData[chunkIndices,]
    ##print(str(sampleData))
    sampleFreqs <- foreach(stim=stimulusVals, .combine=rbind)%do%{
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
    dataAlphas <- c(if(pflag==0){sampleFreqs}else{0})
    ## Generate samples
    set.seed(147+chunk)
    smoothness <- 1000
    smoothm <- sqrt(smoothness)*diff(diag(maxSpikes1),differences=4)
    meangamma <- 0.5
    sdgamma <- 0.5
    shapegamma <- (meangamma/sdgamma)^2
    rategamma <- meangamma/sdgamma^2
    ##
    outfile <- paste0('_LDoutput',chunk)
    thinning <- 2 #10
    mcadapt <- 1e3 #10e3
    mcburnin <- mcadapt
    mciterations <- mcburnin + 2e3
    mcstatus <- 100 #1e3
    Fnames <- paste0('F',rep(nspikesVals,each=nStimuli),'_',rep(stimulusVals,times=maxSpikes1))
    Mnames <- paste0('mean',stimuli)
    stepwidth <- 1 #c(rep(nn/100,length(inu)),rep(nnyR/100,length(iR)))
    nsteps <- 100 #c(rep(nn/100,length(inu)),rep(nnyR/100,length(iR)))
    nprev <- 0
    covm <- list()
    B <- list()
    meansInd <- 1:2
    nspikesVals2 <- rep(nspikesVals,each=nStimuli)
    ##
    PGF <- function(data){
        c(rgamma(n=2,shape=shapegamma, rate=rategamma),
        maxSpikes2 * rdirichlet(n=1, alpha=c(dataAlphas)+1))
    }
    ##
    mydata <- list(y=1, PGF=PGF,
                   parm.names=c(Mnames,Fnames),
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
        means <- parm[meansInd]
        FF <- parm[-meansInd]
        dim(FF) <- c(nstimuli, maxSpikes1)
        sFF <- rowSums(FF)
        tFF <- sum(sFF)
        mi <- mutualinfo(FF/tFF)
        ##
        ##
        lpf <- ddirichlet(x=c(FF)/tFF,
                          alpha = dataAlphas +
                              T * dgeom(nspikesVals2,prob=1/(means+1), log=FALSE),
                          log=TRUE) + 
            sum(dgamma(means, shape=shapegamma, rate=rategamma, log=TRUE)) -
            sum(tcrossprod(FF/sFF, smoothm)^2)
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
    ##
    ##
    mcmcrun <- mcoutput$Posterior1[-round((1:mcburnin)/thinning),-meansInd]
    nDraws <- nrow(mcmcrun)
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
    pdff(paste0('LDplots_samples',nSamples,'_chunk',chunk))
    ##
    matplot(x=nspikesVals, y=t(condfreqSamples[round(seq(1,nDraws,length.out=nPlotSamples)),1,]),
            type='l', lty=1, lwd=2, col=paste0(mypurpleblue,'22'), ylim=c(-1,1),  xlim=c(0,maxX),
            xlab='spikes/bin', ylab='freq', cex.lab=2, cex.axis=2)
    for(i in 2:nStimuli){
        matplot(x=nspikesVals, y=-t(condfreqSamples[round(seq(1,nDraws,length.out=nPlotSamples)),i,]),
                type='l', lty=1, lwd=1, col=paste0(myredpurple,'22'),
                add=TRUE)
    }
    matplot(x=nspikesVals, y=normalize(longrunFreqs[1,]),
            type='l', lty=1, lwd=2, col='black',
            add=TRUE)
    matplot(x=nspikesVals, y=-normalize(longrunFreqs[2,]),
            type='l', lty=1, lwd=2, col='black',
            add=TRUE)
    if(pflag==0){matplot(x=nspikesVals, y=normalize(sampleFreqs[1,]),
                         type='l', lty=2, lwd=5, col=myyellow,
                         add=TRUE)
                         matplot(x=nspikesVals, y=-normalize(sampleFreqs[2,]),
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
