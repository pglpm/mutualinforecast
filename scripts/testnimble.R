## Author: PGL Porta Mana
## Last-Updated: 2021-08-01T08:43:09+0200
################
## Script to test Nimble

library('data.table')
library('khroma')
library('ggplot2')
library('ggthemes')
theme_set(theme_bw(base_size=18))
scale_colour_discrete <- scale_colour_bright
#library('cowplot')
library('png')
library('foreach')
## library('doFuture')
## registerDoFuture()
## library('doRNG')
options(bitmapType='cairo')
pdff <- function(filename){pdf(file=paste0(filename,'.pdf'),paper='a4r',height=11.7,width=16.5)} # to output in pdf format
pngf <- function(filename,res=300){png(file=paste0(filename,'.png'),height=11.7*1.2,width=16.5,units='in',res=res,pointsize=36)} # to output in png format
library('nimble')

deregisterDistributions(c('dmynorm','rmynorm'))
dmynorm <- nimbleFunction(
    run = function(x=double(1), mean=double(1), prec=double(2), log=integer(0, default=0)){ lp <- -sum(asRow(x-mean) %*% prec %*% asCol(x-mean))/2
        if(log) return(lp)
        else return(exp(lp))
        returnType(double(0))
    })
assign('dmynorm',dmynorm,envir=.GlobalEnv)
#assign('dlogsmoothmean', dlogsmoothmean, envir = .GlobalEnv)
#Cdlogsmoothmean <- compileNimble(dlogsmoothmean)
##
code <- nimbleCode({
    mean ~ dnorm(mean=0, sd=0.1)
    sd ~ dgamma(shape=2, rate=1)
    means[1:2] <- c(mean,mean)
    prec[1:2,1:2] <- diag(c(1/sd^2,1/(sd/100)^2))
    
    x[1:2] ~ dmnorm(mean=means[1:2], prec=prec[1:2,1:2])
    y[1:2] ~ dmynorm(mean=means[1:2], prec=prec[1:2,1:2])
})
##
constants <- list()
##
inits <- list(mean=0, sd=1)
##
modeldata <- list()
##
model <- nimbleModel(code=code, name='model', constants=constants, inits=inits, data=modeldata)
Cmodel <- compileNimble(model, showCompilerOutput = TRUE, resetFunctions = TRUE)
##
confmodel <- configureMCMC(Cmodel,nodes=NULL)
confmodel$addSampler(target=c('mean'), type='posterior_predictive')
confmodel$addSampler(target=c('sd'), type='posterior_predictive')
confmodel$addSampler(target=c('x[1:2]','y[1:2]'), type='AF_slice', control=list(sliceAdaptFactorMaxIter=1000, sliceAdaptFactorInterval=100, sliceAdaptWidthMaxIter=500, sliceMaxSteps=100, maxContractions=100))
confmodel$setMonitors(c('x', 'y'))
##confmodel$setMonitors(c('x'))
confmodel
##
mcmc <- buildMCMC(confmodel)
Cmcmc <- compileNimble(mcmc, project= model, resetFunctions = TRUE)
test <- runMCMC(mcmc=Cmcmc,niter=2000,nburnin=1000)


##################################################
library('nimble')
##
dmymnorm <- nimbleFunction(
    run = function(x=double(1), mean=double(1), prec=double(2), log=integer(0, default=0)){ lp <- -sum(asRow(x-mean) %*% prec %*% asCol(x-mean))/2
        if(log) return(lp)
        else return(exp(lp))
        returnType(double(0))
    })
assign('dmymnorm',dmymnorm,envir=.GlobalEnv)
##
code <- nimbleCode({
    mean ~ dnorm(mean=0, sd=1)
    sd ~ dgamma(shape=2, rate=1)
    means[1:2] <- c(mean, mean)
    prec[1:2,1:2] <- diag(c(1/sd^2, 1/(sd/10)^2))
    ##    
    x[1:2] ~ dmymnorm(mean=means[1:2], prec=prec[1:2,1:2])
})
##
constants <- list()
inits <- list(mean=0, sd=1)
modeldata <- list()
##
model <- nimbleModel(code=code, name='model', constants=constants, inits=inits, data=modeldata)
Cmodel <- compileNimble(model, showCompilerOutput = TRUE, resetFunctions = TRUE)
##
confmodel <- configureMCMC(Cmodel)
confmodel$addSampler(target='x', type='AF_slice', control=list(sliceAdaptFactorMaxIter=1000, sliceAdaptFactorInterval=100, sliceAdaptWidthMaxIter=500, sliceMaxSteps=100, maxContractions=100))
confmodel$addMonitors('x')
confmodel
## ===== Monitors =====
## thin = 1: mean, sd, x
## ===== Samplers =====
## posterior_predictive_branch sampler (2)
##   - mean
##   - sd
## AF_slice sampler (1)
##   - x
##
mcmc <- buildMCMC(confmodel)
samples <- runMCMC(mcmc=mcmc,niter=2000,nburnin=1000)
## Warning: running an uncompiled MCMC algorithm, use compileNimble() for faster execution.
## running chain 1...
## Error: user-defined distribution dmymnorm provided without random generation function


confmodel$setMonitors(c('X','logProb_X'))
    confmodel


confmodel <- configureMCMC(Cmodel)
    ## confmodel$addSampler(target='X', type='AF_slice', control=list(sliceAdaptFactorMaxIter=20000, sliceAdaptFactorInterval=1000, sliceAdaptWidthMaxIter=1000, sliceMaxSteps=100, maxContractions=100))
    ## for(i in 1:nStimuli){
    ##     confmodel$addSampler(target=paste0('X[',i,',]'), type='RW_', control=list(sliceAdaptFactorMaxIter=10000, sliceAdaptFactorInterval=500, sliceAdaptWidthMaxIter=500, sliceMaxSteps=100, maxContractions=100))
        ## confmodel$addSampler(target=paste0('X[',i,',]'), type='AF_slice', control=list(sliceAdaptFactorMaxIter=10000, sliceAdaptFactorInterval=500, sliceAdaptWidthMaxIter=500, sliceMaxSteps=100, maxContractions=100))
    ## }
    mcsamples <- runMCMC(Cmcmcsampler, nburnin=20000, niter=20000+50000, thin=10, setSeed=147)
    nDraws <- nrow(mcsamples)
    llsamples <- mcsamples[,-(1:(maxS1*2)),drop=F]
##
condfreqSamples <- mcsamples[,1:(maxS1*nStimuli)]
    ## condfreqSamples <- t(apply(mcsamples[,1:(maxS1*2)],1,function(x){
    ##     dim(x) <- c(2,maxS1)
    ##     Nx2f(x)}))
    dim(condfreqSamples) <- c(nDraws,nStimuli,maxS1)
    dimnames(condfreqSamples) <- list(NULL, nameStimulus, nameNspikes)
    ##
    MIsamples <- apply(condfreqSamples,1,Nmutualinfo)
    MIDistr <- hist(MIsamples, breaks=seq(0,1,by=0.02), plot=F)
    MIQuantiles <- quantile(x=MIsamples, probs=c(0.025,0.5,0.975))
    ##
    meanSsamples <- t(apply(condfreqSamples,1,function(x){x %*% (0:maxS)}))
    ##
    ## PLOTS
    nPlotSamples <- 100
    maxX <- maxS
    if(chunk==0){psign <- 1}else{psign <- -1}
    pdff(paste0('testsummaryNv2random_',chunk))
    matplot(x=nspikesVals, y=t(condfreqSamples[round(seq(1,nDraws,length.out=nPlotSamples)),1,]),
            type='l', lty=1, lwd=2, col=paste0(mygrey,'66'), ylim=c(min(0,psign),1),  xlim=c(0,maxX), xlab='spikes/bin', ylab='freq', cex.lab=2, cex.axis=2)
        matplot(x=nspikesVals, y=psign*t(condfreqSamples[round(seq(1,nDraws,length.out=nPlotSamples)),2,]),
                type='l', lty=1, lwd=1, col=paste0(mygrey,'66'), add=TRUE)
    ##
    if(chunk>0){matplot(x=nspikesVals, y=t(normalizerows(sampleFreqs)*c(1,psign)),
                         type='l', lty=2, lwd=5, col=myyellow, add=TRUE)}
    ##
    matplot(x=nspikesVals, y=t(normalizerows(longrunFreqs)*c(1,psign)),
            type='l', lty=4, lwd=4, col='black', add=TRUE)
    ##
    title(paste0(nSamples,' data samples,',
                 ' chunk =', chunk,
                 ', prior weight = ', sum(prioralphas),
                 '\n superdistr ',chunk), cex.main=2)
    legend('topright',c('long-run freqs','sample freqs'),lty=c(1,2),lwd=c(2,5),col=c('black',myyellow),cex=1.5)
    ##
    ##
    matplot(x=MIDistr$mids, y=MIDistr$density,
            type='h', lty=1, lwd=15, col=paste0(mypurpleblue,'88'), xlim=c(0,1),
            xlab='MI/bit', ylab='prob dens', cex.lab=2, cex.axis=2)
    for(q in MIQuantiles){
        matlines(x=rep(q,2),y=c(-1,1/2)*max(MIDistr$density), lty=2, lwd=6, col=mygreen)
    }
    matlines(x=rep(sampleMI,2),y=c(-1,2/3)*max(MIDistr$density), lty=4, lwd=6, col=myyellow)
    matlines(x=rep(longrunMI,2),y=c(-1,2/3)*max(MIDistr$density), lty=1, lwd=6, col=myredpurple)
    title('predicted MI distr', cex.main=2)
    legend('topright',c('sample MI', 'long-run MI'),lty=1,col=c(myyellow,myredpurple),lwd=4,cex=1.5)
    ##
    ## Diagnostics
    hist(meanSsamples[,1], xlim=c(0,3),ylab='mean spikecountsy')
    hist(meanSsamples[,2], xlim=c(0,3),ylab='mean spikecounts')
    matplot((MIsamples),type='l', lty=1,ylab='MI samples')
    matplot((llsamples[,]),type='l', lty=1,ylab='log-posterior')
    matplot((mcsamples[,1]),type='l', lty=1,ylab='samples of first freq')
    dev.off()
print(paste0('gamma parms: ',c(shapegamma,rategamma)))
NULL
##
##
## plan(sequential)
## plan(multisession, workers = 3L)
## clusterExport(cl=mycluster, c('runcode'))
## alloutput <- parLapply(cl = mycluster, X = 0:2, fun = function(chunk){runcode(chunk)})
## stopCluster(mycluster)

