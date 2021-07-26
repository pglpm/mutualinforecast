## Monte Carlo calculation of posterior for joint frequency distribution
## File with parameters
## *** remember to check the relentropy term in logprob for the "tail"!

## prob. distribution for joint frequencies R and long-run frequencies nu

pdff <- function(filename){pdf(file=paste0(filename,'.pdf'),paper='a4r',height=11.7,width=16.5)}
#library('ggplot2')
#library('RColorBrewer')
#library('cowplot')
#library('png')
#library('plot3D')
#library('Rmpi')
#library('doParallel')
#library('snow')
library('foreach')
library('LaplacesDemon')
#library('MASS')
#library('dplyr')

#cores <- mpi.universe.size() -1
#print(paste0('cores detected: ', cores+1))
## clusters <- 5
## cpus <- 5
## nchains <- 5
nchains <- 7
print(paste0('chains: ', nchains))

filename <- '_mcs_nuR_afssblockcollSnu_adapt2'
scriptname <- 'mcfullnet'
nn <- 2000
nnt <- 75
thinning <- 100
mcadapt <- 2e5
mcburnin <- mcadapt
mciterations <- mcburnin + 1e6
mcstatus <- 1000
L <- 1#5
T <- 417641
dirmult <- 417641*100 #/10
dirbase <- 0 #/10
#prior <- (1-(0:nn)+nn)/((nn+1)*(nn+2))
## 2000: 16:08 - 15:42 

f <- unname(c(unlist(read.csv(paste0('stensola_activity_freqs_bins417641.csv'),header=FALSE,sep=','))))
## f <- unname(c(unlist(read.csv(paste0('stensola_activity_absfreqs_bins417641.csv'),header=FALSE,sep=','))))
nonzero <- which(f>0) ## this must be developed for riehle data

Fme <- unname(c(unlist(read.csv(paste0('stensola_5mom_priorc_N',nn,'.csv'),header=FALSE,sep=','))))
smoothness <- nn/(nnt*8*sum(diff(Fme,differences=4)^2))#round(nn)
Fme <- c(Fme[1:nnt],sum(Fme[(nnt+1):(nn+1)]))
Fme <- Fme/sum(Fme)

ft <- f[f>0]
#ft <- f[1:3]
n <- length(f)-1
nt <- length(ft)-1

## print('creating cluster...')
## if(clusters>1){
## #    cl <- makeCluster(clusters, outfile='/dev/null')
## #     registerDoParallel(cl)
##      registerDoParallel(clusters)
## print(paste0('doparworkers: ', getDoParWorkers()))
## }

## This is the original hyperg.-distribution matrix, truncated for aa>nnt
## hh <- exp(unname(as.matrix(foreach(a=0:nt,.combine=cbind)%:%foreach(aa=0:nnt,.combine=c)%do%{lchoose(n,a)+lchoose(nn-n,aa-a)-lchoose(nn,aa)})))

## This is the new hyperg.-distribution matrix:
## elements are combined and geometrically averaged for aa>=nnt
hh <- rbind(exp(unname(as.matrix(foreach(a=0:nt,.combine=cbind)%:%foreach(aa=0:(nnt-1),.combine=c)%do%{lchoose(n,a)+lchoose(nn-n,aa-a)-lchoose(nn,aa)}))),
            (foreach(a=0:nt,.combine=c)%do%{
                exp(mean(foreach(aa=nnt:(nn-n+a),.combine=c)%do%{
                        lchoose(n,a)+lchoose(nn-n,aa-a)-lchoose(nn,aa)
                    }))*(nn-n+a-nnt+1)
            })/(nn-nnt+1))

## hh2 <- unname(as.matrix(foreach(a=0:nt,.combine=cbind)%:%foreach(aa=0:nnt,.combine=c)%dopar%{choose(n,a)*choose(nn-n,aa-a)/choose(nn,aa)}))

nnt1 <- nnt + 1
nt1 <- nt + 1
nparm <- nnt1 + nt1*(nt1+1)/2 + nt1*(nnt-nt)
print(paste0('number of parameters: ',nparm))
inu <- 1:nnt1
iR <- (nnt1+1):nparm
#i2 <- (nnt+1)+(1:((nt+1)*(nt+2)/2))
#i3 <- (nnt+1)+(nt+1)*(nt+2)/2+(1:((nt+1)*(nnt-nt)))
#ib <- c(i2,i3)
rzeros <- matrix(0,nnt1,nt1)
rpos <- lower.tri(rzeros,diag=T)
nnyR <- n*nn  ##*crossprod(hh,Fme)/ft

rl <- c(nnt1,apply(rpos,2,sum))
##covm <- sapply(rl,function(x){diag(rep(1,x))})

#covm <- diag(rep(nn,nparm))

covm <- list()
B <- list()
sta <- 0
for(i in 1:length(rl)){
    B[[i]] <- (sta+1):(sta+rl[i])
    covm[[i]] <- diag(rep(1,rl[i]))
    sta <- sta+rl[i]
}

stepwidth <- 1 #c(rep(nn/100,length(inu)),rep(nnyR/100,length(iR)))
nsteps <- 1000 #c(rep(nn/100,length(inu)),rep(nnyR/100,length(iR)))
nprev <- 0

## this smoothing-measure matrix leaves the last element independent,
## to be considered as "nnt and more"
smoothm <- sqrt(smoothness)*t(cbind(sapply(1:nnt,function(x){
    diff(diag(nnt)[x,],differences=4)
}),rep(0,nnt-4)))

nunames <- paste0('n',inu-1)
Rnames <- paste0('R',col(rpos)-1,'.',row(rpos)-1)[rpos]

PGF <- function(data){
    ff <- rdirichlet(2,Fme*dirmult+dirbase)
    R <- hh*(ff[2,])
    rR <- colSums(R)/ft
    R <- t(t(R)/rR)
    c(nn*ff[1,],nnyR*R[rpos]) 
}

mydata <- list(y=1, PGF=PGF,
               parm.names=c(nunames,Rnames),
               mon.names=c('lpf')
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

logprob <- function(parm,data){
    parm <- interval(parm,0,Inf)

    nu <- parm[inu]
    rnu <- sum(nu)
    nu <- nu/rnu

    R <- rzeros
    R[rpos] <- parm[iR]
    rR <- colSums(R)/ft
    R <- t(t(R)/rR) # flattening using [rpos] makes it slightly slower

   lpf <- -T*sum(R*log(R/((hh*nu))),na.rm=TRUE)   -sum(crossprod(smoothm, nu)^2) #-sum(crossprod(smoothm, rowSums(R))^2) -L*sum(nu*log(nu),na.rm=TRUE)
##    lpf <- -T*sum(R*log(R/((hh*nu))),na.rm=TRUE)
    LP <- lpf - 5e5*((rnu-nn)^2+sum((rR-nnyR)^2)) ##-jacobian dep. on the rs

    list(LP=LP, Dev=-2*LP, Monitor=lpf, yhat=1, parm=parm)
}

## check logprobs for the maxent distribution

R0 <- hh*Fme
rR0 <- colSums(R0)/ft
R0 <- t(t(R0)/rR0) # flattening using [rpos] makes it slightly slower
parm0 <- c(nn*Fme,n*nn*(R0[rpos]))

lpf0 <- logprob(parm0,mydata)$Monitor
print(paste0('LOGPROB ME DISTR, F SPACE = ',lpf0))
lpx0 <- logprob(parm0,mydata)$LP
print(paste0('LOGPROB ME DISTR, x SPACE = ',lpx0))


fileid <- paste0(filename,'_N',nn,'_t',nnt,'_L',L,'_T',T,'_sm',smoothness,'_s',(mciterations-mcburnin)*nchains/thinning)

print(scriptname)
print(filename)
print(fileid)

print('params read')
