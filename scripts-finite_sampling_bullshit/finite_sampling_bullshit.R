## Calculations of long-run and sample mutual info ##

pdff <- function(filename){pdf(file=paste0(filename,'.pdf'),paper='a4r',height=11.7,width=16.5)}
#library('ggplot2')
library('RColorBrewer')
#library('cowplot')
#library('png')
#library('plot3D')
library('foreach')
library('LaplacesDemon')
#library('RNetCDF')
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

set.seed(181225)



## freqs[,S] = response freqs for stimulus S
## assumes all stimuli equally probable
mutualinfo <- function(freqs,base=2){##in bits
    stimfreqs <- 1/ncol(freqs)
    jointfreqs <- freqs*stimfreqs
    respfreqs <- rowSums(jointfreqs)

    sum(jointfreqs*log2(jointfreqs/outer(respfreqs,rep(stimfreqs,ncol(freqs)))),na.rm=TRUE)/log2(base)
}

freqsampling <- function(freqs,nsamples=20){
    nresps <- nrow(freqs)
    foreach(s=1:ncol(freqs),.combine=cbind)%do%{
        onesample <- sample(1:nresps,nsamples,replace=TRUE,prob=freqs[,s])
        tabulate(onesample,nbins=nresps)
    }    
}


plotsingle <- function(rfreqs0,filename,tit0,nmcsamples=1000,nbreaks='Sturges',base=2){
    nstim <- ncol(rfreqs0)
    nresp <- nrow(rfreqs0)
    milongrun0 <- mutualinfo(rfreqs0)

    misamples0 <- foreach(i=1:nmcsamples,.combine=c)%do%{
        fsample <- foreach(s=1:nstim,.combine=cbind)%do%{
            c(tabulate(sample(x=1:nresp,size=20,replace=TRUE,prob=rfreqs0[,s]),nbins=nresp)/20)
        }
        mutualinfo(fsample)
    }

    hires <- hist(misamples0,breaks=nbreaks,plot=FALSE)
    maxy <- max(hires$counts)
    minx <- -0.02
    maxx <- 1.02
    ## minx <- min(hires$breaks,milongrun0)-0.02
    ## maxx <- max(hires$breaks,milongrun0)+0.02
    subs <- misamples0[seq(1,length(misamples0),length.out=100)]
    smean <- signif(mean(misamples0),3)
    q1 <- signif(quantile(misamples0,0.16),3)
    q2 <- signif(quantile(misamples0,0.84),3)
    q3 <- signif(quantile(misamples0,0.025),3)
    q4 <- signif(quantile(misamples0,0.975),3)
    smedian <- signif(quantile(misamples0,0.5),3)

    maintext <- paste0('long-run=',signif(milongrun0,3),' bit;  sample: mean=',smean,' bit,  median=',smedian,' bit,  68% in (',q1,', ',q2,') bit,  95% in (',q3,', ',q4,') bit')
    
    pdff(paste0('histo_',filename))
    matplot(x=subs,y=rep(-maxy/20,100),type='p',lty=1,lwd=3,pch=18,cex=1,col=myblue,xlim=c(minx,maxx),ylim=c(-maxy/20,maxy),xlab=paste0('sampled MI/bit'),ylab='',main=maintext)
    hist(misamples0,breaks=nbreaks,xlim=c(minx,maxx),ylim=c(-maxy/20,maxy),xlab='I',ylab='',add=TRUE)
##legend('top',paste0('\nblack: mean, blue: median\nyellow: 16%q'),bty='n')
    matlines(x=rep(milongrun0,2),y=c(-maxy/20,maxy),type='l',lty=1,lwd=3,pch='.',col=myred)
    dev.off()

    pdff(paste0('resp_',filename))
    barplot(t(freqs0),beside=TRUE,xlab='responses',ylab='long-run frequencies',main=tit0,names=1:10,space=c(0.2,1))
    dev.off()
}

freqs0 <- matrix(1/10,10,2)
plotsingle(freqs0,'caseA', 'case A',nmcsamples=10000,nbreaks='Sturges')

pt <- 6
fr <- 100
f1 <- (c((dcoe(0,n=pt))*fr/100, rep(1/(10-pt),(10-pt))*(100-fr)/100))
freqs0 <- cbind(rev(f1),f1)
plotsingle(freqs0,'caseB', 'case B',nmcsamples=10000,nbreaks='Sturges')

pt <- 6
fr <- 90
f1 <- c((dcoe(0,n=pt))*fr/100, rep(1/(10-pt),(10-pt))*(100-fr)/100)
pt <- 6
fr <- 90
f2 <- c(rev(dcoe(0,n=pt))*fr/100, rep(1/(10-pt),(10-pt))*(100-fr)/100)
freqs0 <- cbind(rev(f1),(f2))
plotsingle(freqs0,'caseC', 'case C',nmcsamples=10000,nbreaks='Sturges')

pt <- 5
fr <- 99
f1 <- c(rev(dcoe(1/pt,n=pt))*fr/100, rep(1/(10-pt),(10-pt))*(100-fr)/100)
freqs0 <- cbind(rev(f1),(f1))
plotsingle(freqs0,'caseD', 'case D',nmcsamples=10000,nbreaks='Sturges')



f1 <- c(rep(1/5,5)*90/100,10/100,rep(0,4))
freqs0 <- cbind(f1,rev(f1))
plotsingle(freqs0,'caseB', 'case B',nmcsamples=10000,nbreaks='Sturges')

f1 <- c(rep(1/5,5),rep(0,5))
freqs0 <- cbind(f1,rev(f1))
plotsingle(freqs0,'caseB', 'case B')


dcoe <- function(a,n=10){(1-n*a)*2/(n^2+n)*(1:n)+a}
f1 <- dcoe(0)
freqs0 <- cbind(f1,rev(f1))
plotsingle(freqs0,'caseC', 'case C')


f1 <- dcoe(0)
freqs0 <- cbind(rev(f1),rev(f1))
plotsingle(freqs0,'test', 'case B')

pt <- 10
fr <- 100
f1 <- c(rev(dcoe(1/pt,n=pt))*fr/100, rep(1/(10-pt),(10-pt))*(100-fr)/100)
freqs0 <- cbind(rev(f1),(f1))
plotsingle(freqs0,'caseA', 'case A',nmcsamples=10000,nbreaks='Sturges')


pt <- 7
fr <- 100
f1 <- c(rev(dcoe(0,n=pt))*fr/100, rep(1/(10-pt),(10-pt))*(100-fr)/100)
freqs0 <- cbind(rev(f1),f1)
plotsingle(freqs0,'test', 'case C')








f1 <- c((dcoe(0,n=6)), rep(0,4))
freqs0 <- cbind(f1,rev(f1))
plotsingle(freqs0,'test', 'case C')


plotmulti <- function(afreqs,filename,tit0,nmcsamples=1000,base=2){
    nstim <- 2
    if(is.null(dim(afreqs))){afreqs <- cbind(afreqs,afreqs)}

    misamples <- foreach(i=1:nmcsamples,.combine=rbind)%do%{
        rfreqs0 <- cbind(c(rdirichlet(n=1,alpha=afreqs[,1])),
                         c(rdirichlet(n=1,alpha=afreqs[,2])))
        
        fsample <- foreach(s=1:nstim,.combine=cbind)%do%{
            tabulate(sample(x=1:nresp,size=20,replace=TRUE,prob=rfreqs0[,s]),nbins=nresp)/20
        }
        c(mutualinfo(rfreqs0), mutualinfo(fsample))
    }

    maxmi <- 1 #max(c(misamples))
    pdff(paste0('scatter_',filename))
    par(pty = "s")
    matplot(x=misamples[,1],y=misamples[,2],type='p',lty=1,lwd=3,pch=16,cex=0.2,col=myblue,xlim=c(0,maxmi),ylim=c(0,maxmi),xlab=paste0('long-run MI/bit'),ylab=paste0('sample MI/bit'),main=tit0)
    matlines(x=c(0,maxmi),y=c(0,maxmi),type='l',lty=2,lwd=2,pch='.',col=myyellow)
    dev.off()
}

plotmulti(dcoe(0)*5,filename='test',tit0='',nmcsamples=5000)

plotmulti(rep(10,10),filename='centrepeak',tit0='')

plotmulti(rep(1,10),filename='unif',tit0='')

plotmulti(rep(0.1,10),filename='jeffr',tit0='')

plotmulti(rep(1,10),filename='test',tit0='',nmcsamples=5000)

plotmulti(rep(0.1,10),filename='test',tit0='',nmcsamples=10000)

plotmulti(10*rev(dcoe(0)),filename='test',tit0='',nmcsamples=10000)



f1 <- dcoe(0)*1
plotmulti(cbind(f1,rev(f1)),filename='test',tit0='')

plotmulti(rep(1,10),filename='unif',tit0='')

plotmulti(rep(0.1,10),filename='jeffr',tit0='')



plotmultimean <- function(afreqs,filename,tit0,nmcsamples=10,nsubsamples=500,base=2){
    nstim <- 2
    if(is.null(dim(afreqs))){afreqs <- cbind(afreqs,afreqs)}

    misamples <- foreach(i=1:nmcsamples,.combine=rbind)%do%{
        rfreqs0 <- cbind(c(rdirichlet(n=1,alpha=afreqs[,1])),
                         c(rdirichlet(n=1,alpha=afreqs[,2])))

        meanmi <- mean(foreach(i=1:nsubsamples,.combine=c)%do%{
            fsample <- foreach(s=1:nstim,.combine=cbind)%do%{
                tabulate(sample(x=1:nresp,size=20,replace=TRUE,prob=rfreqs0[,s]),nbins=nresp)/20
            }
            mutualinfo(fsample)
        })
        c(mutualinfo(rfreqs0), meanmi)
    }

    maxmi <- max(c(misamples))
    pdff(paste0('scatter_',filename))
    par(pty = "s")
    matplot(x=misamples[,1],y=misamples[,2],type='p',lty=1,lwd=3,pch='.',cex=4,col=myblue,xlim=c(0,maxmi),ylim=c(0,maxmi),xlab=paste0('long-run MI/bit'),ylab=paste0('sampled MI/bit'),main=tit0)
    matlines(x=c(0,maxmi),y=c(0,maxmi),type='l',lty=2,lwd=2,pch='.',col=myyellow)
    dev.off()
}



plotmulti(rep(10,10),filename='centrepeak',tit0='')

plotmulti(rep(1,10),filename='unif',tit0='')

plotmulti(rep(0.1,10),filename='jeffr',tit0='')

plotmultimean(rep(0.1,10),filename='test',tit0='',nmcsamples=30,nsubsamples=500)


plotmulti(dcoe(0)*5,filename='test',tit0='')

f1 <- dcoe(0)*1
plotmulti(cbind(f1,rev(f1)),filename='test',tit0='')

plotmulti(rep(1,10),filename='unif',tit0='')

plotmulti(rep(0.1,10),filename='jeffr',tit0='')




plotmultih <- function(afreqs0,stre,filename,tit0,nmcsamples=1000,base=2){
    nstim <- 2

    misamples <- foreach(i=1:nmcsamples,.combine=rbind)%do%{
        afreqs <- stre*t(rdirichlet(n=2,alpha=afreqs0))
        
        rfreqs0 <- cbind(c(rdirichlet(n=1,alpha=afreqs[,1])),
                         c(rdirichlet(n=1,alpha=afreqs[,2])))
        
        fsample <- foreach(s=1:nstim,.combine=cbind)%do%{
            tabulate(sample(x=1:nresp,size=20,replace=TRUE,prob=rfreqs0[,s]),nbins=nresp)/20
        }
        c(mutualinfo(rfreqs0), mutualinfo(fsample))
    }

    maxmi <- 1 #max(c(misamples))
    pdff(paste0('scatter_',filename))
    par(pty = "s")
    matplot(x=misamples[,1],y=misamples[,2],type='p',lty=1,lwd=3,pch=16,cex=0.2,col=myblue,xlim=c(0,maxmi),ylim=c(0,maxmi),xlab=paste0('long-run MI/bit'),ylab=paste0('sample MI/bit'),main=tit0)
    matlines(x=c(0,maxmi),y=c(0,maxmi),type='l',lty=2,lwd=2,pch='.',col=myyellow)
    dev.off()
}

plotmultih(rep(1,10),stre=10,filename='test',tit0='',nmcsamples=10000)

plotmultih(rep(1,10),stre=1,filename='test',tit0='',nmcsamples=2000)

plotmultih(rep(1,10),stre=1,filename='test',tit0='',nmcsamples=2000)


plotmulti(dcoe(0)*5,filename='test',tit0='')

f1 <- dcoe(0)*1
plotmulti(cbind(f1,rev(f1)),filename='test',tit0='')



backmultih <- function(afreqs0,stre,samplevalue,nmcsamples=1000,base=2){
    nstim <- 2
    v1 <- samplevalue[1]
    v2 <- samplevalue[2]

    misamples <- foreach(i=1:nmcsamples,.combine=c)%dopar%{
        afreqs <- stre*t(rdirichlet(n=2,alpha=afreqs0))
        
        rfreqs0 <- cbind(c(rdirichlet(n=1,alpha=afreqs[,1])),
                         c(rdirichlet(n=1,alpha=afreqs[,2])))
        
        fsample <- foreach(s=1:nstim,.combine=cbind)%do%{
            tabulate(sample(x=1:nresp,size=20,replace=TRUE,prob=rfreqs0[,s]),nbins=nresp)/20
        }
        smi <- mutualinfo(fsample)
        if(smi<v1 || smi>v2) return(NULL)
        mutualinfo(rfreqs0)
    }
    misamples
}


testsa <- backmultih(rep(1,10),stre=1,c(0.19,0.21),nmcsamples=100000)

saveRDS(testsa,file='example_backmi.rds')



plotmultihlines <- function(afreqs0,stre,filename,vsample,vlr,tit0,nmcsamples=1000,base=2){
    nstim <- 2

    misamples <- foreach(i=1:nmcsamples,.combine=rbind)%do%{
        afreqs <- stre*t(rdirichlet(n=2,alpha=afreqs0))
        
        rfreqs0 <- cbind(c(rdirichlet(n=1,alpha=afreqs[,1])),
                         c(rdirichlet(n=1,alpha=afreqs[,2])))
        
        fsample <- foreach(s=1:nstim,.combine=cbind)%do%{
            tabulate(sample(x=1:nresp,size=20,replace=TRUE,prob=rfreqs0[,s]),nbins=nresp)/20
        }
        c(mutualinfo(rfreqs0), mutualinfo(fsample))
    }

    maxmi <- 1 #max(c(misamples))
    pdff(paste0('scatter_',filename))
    par(pty = "s")
    matplot(x=misamples[,1],y=misamples[,2],type='p',lty=1,lwd=3,pch=16,cex=0.2,col=myblue,xlim=c(0,maxmi),ylim=c(0,maxmi),xlab=paste0('long-run MI/bit'),ylab=paste0('sample MI/bit'),main=tit0)
    matlines(x=c(-0.1,vlr[2]),y=c(vsample,vsample),type='l',lty=2,lwd=2,pch='.',col=mygreen)
    matlines(x=rep(vlr[1],2),y=c(vsample,-0.1),type='l',lty=2,lwd=2,pch='.',col=mygreen)
    matlines(x=rep(vlr[2],2),y=c(vsample,-0.1),type='l',lty=2,lwd=2,pch='.',col=mygreen)
    dev.off()
}

plotmultihlines(rep(1,10),stre=1,vsample=0.2,vlr=quantile(testsa,probs=c(0.16,0.84)),filename='test',tit0='',nmcsamples=10000)
