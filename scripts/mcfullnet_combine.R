## prob. distribution for joint frequencies R and long-run frequencies nu
inputname <- commandArgs(trailingOnly=TRUE)
set.seed(181225)

#### Read parameters ####
source(paste0('_sc_params_',inputname,'.R'))

#### Read Monte Carlo outputs ####
mcoutputp <- foreach(jobid=1:nchains)%do%{
    readRDS(file=paste0('_mcoutput',jobid,'_',inputname,'.rds'))
}
print('...finished reading MC output')

elim <- (1:(mcburnin/thinning)) ## what to retain if no stat samples

samples <- foreach(i=1:nchains,.combine='rbind')%do%{
        if(anyNA(mcoutputp[[i]]$Posterior2) || is.null(dim(mcoutputp[[i]]$Posterior2))){
            print('Using posterior 1')
            mcoutputp[[i]]$Posterior1[-elim,]
        }else{
            print('Using posterior 2')
            mcoutputp[[i]]$Posterior2
        }
    }

lpx <- -(foreach(i=1:nchains,.combine='c')%do%{mcoutputp[[i]]$Deviance})/2
if(mcburnin>0){
    lpfi <- foreach(i=1:nchains,.combine='c')%do%{mcoutputp[[i]]$Monitor[elim]}}
lpff <- foreach(i=1:nchains,.combine='c')%do%{mcoutputp[[i]]$Monitor[-elim]}



#samplesf <- t(apply(samples,1,function(v){v0 <- c(-sum(v),c(v/sum(v)}))
##isamples <- invlogit(samples)


## transformation from parameters to frequencies and joint frequencies
renormR <- function(Rvec){
    r <- rzeros
    r[rpos] <- Rvec
    (t(t(r)*(ft/colSums(r))))[rpos]
}

FfromR <- function(Rvec){
    r <- rzeros
    r[rpos] <- Rvec
    rowSums(r)
}

nsamples <- nrow(samples)

samplesnu <- samples[,inu]
samplesnu <- samplesnu/rowSums(samplesnu)

samplesR <- t(apply(samples[,iR],1,renormR))

samplesF <- t(apply(samplesR,1,FfromR))

namedata <- paste0(fileid,'_ss',nsamples)

##file.rename(filetemp,paste0('_mcoutput',namedata,'.rds'))

## calculate, mean, quantiles, std
sstatnu <- apply(samplesnu,2,function(v){c(mean(v),quantile(v,c(50,16,84)/100),sd(v)*(nsamples+1)/nsamples)})

sstatR <- apply(samplesR,2,function(v){c(mean(v),quantile(v,c(50,16,84)/100),sd(v)*(nsamples+1)/nsamples)})

sstatF <- apply(samplesF,2,function(v){c(mean(v),quantile(v,c(50,16,84)/100),sd(v)*(nsamples+1)/nsamples)})


print('writing samples to file...')

## mean, quantiles, std
write.table(sstatnu,file=paste0('_statsnu',namedata,'.csv'),row.names=FALSE,col.names=FALSE, sep=',')
write.table(sstatR,file=paste0('_statsR',namedata,'.csv'),row.names=FALSE,col.names=FALSE, sep=',')
write.table(sstatF,file=paste0('_statsF',namedata,'.csv'),row.names=FALSE,col.names=FALSE, sep=',')

## 1000 samples
write.table(samplesnu[round(seq(1,nsamples,length.out=min(1000,nsamples))),],file=paste0('_nusamples',namedata,'.csv'),row.names=FALSE,col.names=FALSE, sep=',')

write.table(samplesF[round(seq(1,nsamples,length.out=min(1000,nsamples))),],file=paste0('_Fsamples',namedata,'.csv'),row.names=FALSE,col.names=FALSE, sep=',')

## log-probabilities of X
write.table(lpx,file=paste0('_lpx',namedata,'.csv'),row.names=FALSE,col.names=FALSE, sep=',')

## log-probabilities of F
if(mcburnin>0){write.table(lpfi,file=paste0('_lpfi',namedata,'.csv'),row.names=FALSE,col.names=FALSE, sep=',')}
write.table(lpff,file=paste0('_lpff',namedata,'.csv'),row.names=FALSE,col.names=FALSE, sep=',')

saveRDS(samplesF,file=paste0('_Fallsamples',namedata,'.rds'))
saveRDS(samplesnu,file=paste0('_nuallsamples',namedata,'.rds'))

print('plots...')

pdff(paste0('_plotF',namedata))
matplot(x=(0:nnt)/nn,y=t(rbind(sstatF[-5,],Fme))*nn,type='l',lty=c(1,1,1,1,5),lwd=c(3,3,1,1,3),pch='.',col=c('blue','green','yellow3','yellow3','red'),xlab='A',ylab='F')
dev.off()

pdff(paste0('_plotnu',namedata))
matplot(x=(0:nnt)/nn,y=t(rbind(sstatnu[-5,],Fme))*nn,type='l',lty=c(1,1,1,1,5),lwd=c(3,3,1,1,3),pch='.',col=c('blue','green','yellow3','yellow3','red'),xlab='A',ylab='F')
dev.off()

if(mcburnin>0){pdff(paste0('_lpfi',namedata))
matplot(y=lpfi,type='p',lty=1,lwd=1,pch='.',col='black',xlab='A',ylab='F')
dev.off()}
pdff(paste0('_lpff',namedata))
matplot(y=lpff,type='p',lty=1,lwd=1,pch='.',col='black',xlab='A',ylab='F')
dev.off()

pdff(paste0('_numcsamples',namedata))
matplot(x=(0:nnt)/nn,y=t(samplesnu[round(seq(1,nsamples,length.out=min(100,nsamples))),])*nn,type='l',lty=1,lwd=1,pch='.',col='grey',xlab='A',ylab='F')
dev.off()

pdff(paste0('_Fmcsamples',namedata))
matplot(x=(0:nnt)/nn,y=t(samplesF[round(seq(1,nsamples,length.out=min(100,nsamples))),])*nn,type='l',lty=1,lwd=1,pch='.',col='grey',xlab='A',ylab='F')
dev.off()

print('done with plots')
