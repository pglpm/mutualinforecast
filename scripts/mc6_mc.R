## Monte Carlo calculation of posterior probability ##

jobid <- as.numeric(commandArgs(trailingOnly=TRUE))
if(length(jobid)==0){jobid <- 1}

pdff <- function(filename){pdf(file=paste0(filename,'.pdf'),paper='a4r',height=11.7,width=16.5)}
#library('ggplot2')
library('RColorBrewer')
#library('cowplot')
#library('png')
#library('plot3D')
library('foreach')
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

set.seed(181225+jobid)

scriptname <- 'mc6'

#################################################################
######################## EQUATIONS SETUP ########################
#################################################################

generalname <- paste0('_f_dirich1')
datafilename <- 'data1'
priorflag <- 0 # 1=NOT prior, 0=prior


nn <- 10L
print(paste0('number of bins/session: ',nn))
bins <- seq(0,1,length.out=nn+1)
binnames <- 1L:nn
bincentres <- (bins[-1] + bins[-(nn+1)])/2
dimnn <- c(nn,nn,nn)
binwidth <- 1L/nn

data0 <- read.csv(paste0(datafilename,'.csv'),header=TRUE,sep=',')
data0 <- data0[,c('Session','Successful.trials')]
dataL <- sapply(1:3,function(i){data0[data0[,'Session']==i,'Successful.trials']})
binL <- sapply(1:3,function(i){as.integer(cut(dataL[,i],breaks=bins,include.lowest = TRUE, labels=binnames))})

nn3 <- nn^3
freqL <- rep(0,nn3)
dim(freqL) <- dimnn
for(i in 1:nrow(binL)){
    freqL[t(binL[i,])]  <- freqL[t(binL[i,])] + 1
}

freqL <- priorflag * freqL

## data0 <- read.csv(paste0('unlearned.csv'),header=TRUE,sep=',')
## data0 <- data0[,c('Session','Successful.trials')]
## dataU <- sapply(1:3,function(i){data0[data0[,'Session']==i,'Successful.trials']})
## binU <- sapply(1:3,function(i){as.integer(cut(dataU[,i],breaks=bins,include.lowest = TRUE, labels=binnames))})
## freqU <- foreach(i=1:3, .combine=cbind)%do%{tabulate(binU[,i], nbins=nn)}

## data0 <- read.csv(paste0('observers.csv'),header=TRUE,sep=',')
## data0 <- data0[,c('Session','Successful.trials')]
## dataX <- sapply(1:3,function(i){data0[data0[,'Session']==i,'Successful.trials']})
## binX <- sapply(1:3,function(i){as.integer(cut(dataX[,i],breaks=bins,include.lowest = TRUE, labels=binnames))})
## freqX <- foreach(i=1:3, .combine=cbind)%do%{tabulate(binX[,i], nbins=nn)}


baseL <- dnorm(bincentres,0.9,0.1)
baseL <- baseL/sum(baseL)
matplot(x=bincentres,y=baseL,type='l',lty=c(1,1,1,1,1),lwd=c(3,3,1,1,1),pch='.',col=c('black',myblue,myyellow,myyellow,myred),xlab='A/N',ylab='N F')
## baseU <- dnorm(bincentres,bincentres[1],0.2)
## baseU <- baseU/sum(baseU)
## ## matplot(x=bincentres,y=baseU,type='l',lty=c(1,1,1,1,1),lwd=c(3,3,1,1,1),pch='.',col=c('black',myblue,myyellow,myyellow,myred),xlab='A/N',ylab='N F')

diff4 <- function(x){for(i in 1:4){x <- x[-1]-x[-length(x)]}
    x}
sc <- 1/10
ac <- 1 #nn^(1/2)

mylegend <- paste0('D^3lnnu, S=',sc,', A=',ac)

alpha <- ac*(baseL %o% baseL %o% baseL)
initdistr <- alpha*nn3


nparm <- nn3 - 1
print(paste0('number of parameters: ',nparm))

pquantum <- (min(baseL)^3)/10
unifd <- rep(1/nn3,nn3)

X2F <- function(x){
    x <- c(x,-sum(x))
    x <- exp(x-max(x))
    x <- x/sum(x)
    dim(x) <- dimnn
    x
}

F2X <- function(probd){
    probd <- c(probd)
    if(min(probd)==0){probd <- (probd+unifd*pquantum)/(1+pquantum)}
    lnu <- log(probd)
    (lnu-mean(lnu))[-nn3]
}

pe2 <- c(2,1,3)
pe3 <- c(3,2,1)

marg1 <- function(x){dim(x) <- dimnn
    rowSums(x)}
marg2 <- function(x){dim(x) <- dimnn
    rowSums(aperm(x,pe2))}
marg3 <- function(x){dim(x) <- dimnn
    rowSums(aperm(x,pe3))}

###################################################################
######################## MONTE CARLO SETUP ########################
###################################################################
## Initial.Values <- as.initial.values(mcoutput)
## mcoutput <- LaplacesDemon(Model, Data=mydata, Initial.Values,
##      Covar=mcoutput$Covar, Iterations=1350, Status=6, Thinning=45,
##      Algorithm="AFSS", Specs=list(A=50, B=NULL, m=mcoutput$Specs$m,
##      n=150, w=mcoutput$CovarDHis[nrow(mcoutput$CovarDHis),]))

adaptphase <- TRUE
thinning <- round(nparm/5)
mciterations <- thinning*50
mcstatus <- thinning*1
singleplots <- FALSE
combinedplots <- TRUE

onlyplots <- FALSE

filesid <- paste0(generalname,'_N',nn)#,'_sm',smc,'_av',ac) #,'_samp',(mciterations-mcburnin)*nchains/thinning)

stagecount <- 0
if(length(list.files(pattern=paste0('^_output[0-9]*_stg',stagecount,filesid,'.rds')))==0){
    ## first run

    if(adaptphase){mcadapt <- Inf} else {mcadapt <- 0}

    listparm <- c(1:nparm,NA)
    dim(listparm) <- dimnn
    condt <- ((1:10) - 1)%/%5
    covm <- list()
    B <- list()
    ll <- 0
    for(i in 0:1){for(j in 0:1){for(k in 0:1){ll <- ll+1
                                    toinsert <- c(listparm[condt==i,condt==j,condt==k])
                                    B[[ll]] <- toinsert[!is.na(toinsert)]
                                    covm[[ll]] <- diag(rep(1,length(B[[ll]])))
                            }}}
    stepwidth <- 1 # for AFSS
    nsteps <- 100
    nprev <- 0
    PGF <- function(data){
        ff <- c(rdirichlet(1,initdistr))
        F2X(ff)
    }
} else {
    ## subsequent runs
    
    stagecount <- 1
    while(length(list.files(pattern=paste0('^_output[0-9]*_stg',stagecount,filesid,'.rds'))) > 0){stagecount <- stagecount+1}

    prevoutputlist <- list.files(pattern=paste0('^_output[0-9]*','_stg',stagecount-1,filesid,'.rds'))
    prevnchains <- length(prevoutputlist)

    prevjobid <- 1 + (jobid-1)%%prevnchains

    prevoutput <- paste0('_output',prevjobid,'_stg',stagecount-1,filesid,'.rds')
    print(paste0('using ',prevoutput))
    mcoutput <- readRDS(prevoutput)

    inivalueid <- nrow(mcoutput$Posterior1) - ((jobid-1)%/%prevnchains)
    print(paste0('using previous sample #',inivalueid))

    if(adaptphase){mcadapt <- Inf} else {mcadapt <- 100}
    B <- mcoutput$Specs$B
    covm <- mcoutput$Covar # for AFSS
    stepwidth <- mcoutput$CovarDHis[nrow(mcoutput$CovarDHis),] # for AFSS
    nsteps <- mcoutput$Specs$m
    nprev <- mcoutput$Iterations + mcoutput$Specs$n
    PGF <- function(data){
        mcoutput$Posterior1[inivalueid,]
    }
}

if(onlyplots){stagecount <- stagecount - 1}
stagename <- paste0('_stg',stagecount)
savefilesid <- paste0(jobid,stagename,filesid)

vnames <- rep(NA,nn3)
dim(vnames) <- dimnn
for(i in 1:nn){for(j in 1:nn){for(k in 1:nn){
                                  vnames[i,j,k] <- paste0(i,'.',j,'.',k)}}}
vnames <- c(vnames)[-1]

mydata <- list(y=1, PGF=PGF,
               parm.names=c(vnames),
               mon.names=c('lpf')
               )

######################## LOG-PROB ########################
newalpha <- alpha + freqL
sc2 <- 1/(2*sc^2)

logprob <- function(parm,data){
    X <- c(parm, -sum(parm))
    X <- exp(X-max(X))
    X <- X/sum(X)
    dim(X) <- dimnn
    LP <- sum(newalpha * log(X)) - sc2*(sum(diff4(rowSums(X))^2) +
                                        sum(diff4(rowSums(colSums(X)))^2) +
                                        sum(diff4(colSums(X,dim=2))^2))
    
    list(LP=LP, Dev=-2*LP, Monitor=1, yhat=1, parm=parm)
}

## probmode <- logprob(F2X(Fmax),mydata)$Monitor
## print(paste0('logprob mode: ', probmode))

print(paste0('test logprob: ', logprob(PGF(mydata),mydata)$LP))

## check logprobs for the maxent distribution

print('...parameters read.')

print(paste0('scriptname: ',scriptname))
print(paste0('generalname: ',generalname))
print(paste0('savefilesid: ',savefilesid))
print(paste0('ADAPTIVE: ',adaptphase))
print(paste0('stage: ',stagecount))
print(paste0('total iterations so far: ',nprev))
print(paste0('iterations to do: ',mciterations))
print(paste0('thinning: ',thinning))
print(paste0('samples: ',round(mciterations/thinning)))


## copy scripts for later checks
if(jobid==1 && !onlyplots){
    file.copy(paste0(scriptname,'_mc.R'),paste0('_scr_mc',savefilesid,'.R'),overwrite=TRUE)
##    file.copy(paste0(scriptname,'_combine.R'),paste0('_scr_combine',stagename,filesid,'.R'),overwrite=TRUE)
    file.copy(paste0(scriptname,'_call.slurm'),paste0('_scr_call',savefilesid,'.slurm'),overwrite=TRUE)
    if(stagecount==0){## very first run
            file.copy(paste0(scriptname,'_mc.R'),paste0('_master_scr_mc',filesid,'.R'),overwrite=TRUE)
##    file.copy(paste0(scriptname,'_combine.R'),paste0('_master_scr_combine',filesid,'.R'),overwrite=TRUE)
    file.copy(paste0(scriptname,'_call.slurm'),paste0('_master_scr_call',filesid,'.slurm'),overwrite=TRUE)
    }
}

#################################################################
######################## MONTE CARLO RUN ########################
#################################################################

if(stagecount==0){
    Initial.Values <- GIV(logprob, mydata, n=1000000, PGF=TRUE)
} else {Initial.Values <- PGF(mydata)}

if(!onlyplots){
print('Running Monte Carlo...')
mcoutput <- LaplacesDemon(logprob, mydata, Initial.Values,
                      ##Covar=NULL,
                      Covar=covm,
                      Thinning=thinning,
                      Iterations=mciterations, Status=mcstatus,
                      Debug=list(DB.chol=TRUE, DB.eigen=TRUE,
                                 DB.MCSE=TRUE, DB.Model=FALSE),
                      LogFile=paste0('_LDlog',savefilesid,'.txt'), #Type="MPI", #Packages=c('MASS'),
                      ##Algorithm="RDMH"#, Specs=list(B=list(1:d,d1:d2,d3:dnp))
                      ##Algorithm="Slice", Specs=list(B=NULL, Bounds=c(0,1), m=Inf, Type="Continuous", w=0.001)
                      ##Algorithm="MALA", Specs=list(A=1e7, alpha.star=0.574, gamma=mcadapt, delta=1, epsilon=c(1e-6,1e-7))
##                      Algorithm="pCN", Specs=list(beta=0.01,n=nprev)
##                      Algorithm="RWM", Specs=list(B=B,n=nprev)
##                      Algorithm="pCN", Specs=list(beta=0.01,n=nprev)
                      Algorithm="AFSS", Specs=list(A=mcadapt, B=B, m=nsteps, n=nprev, w=stepwidth)
                      ##Algorithm="NUTS", Specs=list(A=mcadapt, delta=0.6, epsilon=stepwidth, Lmax=Inf, n=nprev)
                      ##Algorithm="AIES", Specs=list(Nc=4*nparm, Z=NULL, beta=2, CPUs=1, Packages=NULL, Dyn.libs=NULL)
                      
                      ##Algorithm="DRM", Specs=NULL
                      )
    
    ## RWM:no HARM:no pCN:no
    ## mcoutput <- Combine(mcoutputp,mydata)

#    cl <- makeCluster(clusters, outfile='/dev/null')
#     registerDoParallel(cl)

    
    ## save all output for later computation
    saveRDS(mcoutput,file=paste0('_output',savefilesid,'.rds'))
}

###############################################################
######################## PLOTS & DIAGN ########################
###############################################################

## genal function to output plots and diagnostics
plotsanddiagn <- function(mcoutputx,savename,savenc=TRUE){

    qfunction <- function(v){c(mean(v),quantile(v,c(50,16,84)/100),sd(v)*(nsamples+1)/nsamples)}

    samples <- mcoutputx$Posterior1

    lpff <- -mcoutputx$Dev/2

    nsamples <- nrow(samples)

    sstato <- apply(samples,2,qfunction)
    print('means:')
    print(sstato[1,])
    print('sds:')
    print(sstato[5,])
    print('medians:')
    print(sstato[2,])
    write.table(sstato,file=paste0('_statsorig',savename,'.csv'),row.names=FALSE,col.names=FALSE, sep=',')
    
    ## transformation from parameters to frequencies

    samplesF <- t(apply(samples,1,function(x){c(X2F(x))}))
    colnames(samplesF) <- c(paste0('F',vnames), paste0('F',nn,'.',nn,'.',nn))

    samplesM1 <- t(apply(samplesF,1,marg1))
    colnames(samplesM1) <- paste0('M1.',1:nn)

    samplesM2 <- t(apply(samplesF,1,marg2))
    colnames(samplesM2) <- paste0('M2.',1:nn)

    samplesM3 <- t(apply(samplesF,1,marg3))
    colnames(samplesM3) <- paste0('M3.',1:nn)

    tocheck <- samplesF
    
    ## calculate, mean, quantiles, std
    sstatF <- apply(samplesF,2,qfunction)
    write.table(sstatF,file=paste0('_statsF',savename,'.csv'),row.names=FALSE,col.names=FALSE, sep=',')

    sstatM1 <- apply(samplesM1,2,qfunction)
    write.table(sstatM1,file=paste0('_statsM1',savename,'.csv'),row.names=FALSE,col.names=FALSE, sep=',')

    sstatM2 <- apply(samplesM2,2,qfunction)
    write.table(sstatM2,file=paste0('_statsM2',savename,'.csv'),row.names=FALSE,col.names=FALSE, sep=',')

    sstatM3 <- apply(samplesM3,2,qfunction)
    write.table(sstatM3,file=paste0('_statsM3',savename,'.csv'),row.names=FALSE,col.names=FALSE, sep=',')

    if(savenc){print('Saving data as netCDF')

        nc <- create.nc(paste0('_dataF',savename,'.nc'))
        dim.def.nc(nc, 'N3', nn3)
        dim.def.nc(nc, 'N', nn)
        dim.def.nc(nc, 'sample', nsamples)
        dim.def.nc(nc, 'stats', 5)
        
        var.def.nc(nc, 'N3', 'NC_INT', 'N3')
        var.def.nc(nc, 'N', 'NC_INT', 'N')
        var.def.nc(nc, 'sample', 'NC_INT', 'sample')
        var.def.nc(nc, 'stats', 'NC_INT', 'stats')
        
        var.def.nc(nc, 'F', 'NC_DOUBLE', c('sample','N3'))
        var.def.nc(nc, 'Fstats', 'NC_DOUBLE', c('stats','N3'))
        var.def.nc(nc, 'M1', 'NC_DOUBLE', c('sample','N'))
        var.def.nc(nc, 'M1stats', 'NC_DOUBLE', c('stats','N'))
        var.def.nc(nc, 'M2', 'NC_DOUBLE', c('sample','N'))
        var.def.nc(nc, 'M2stats', 'NC_DOUBLE', c('stats','N'))
        var.def.nc(nc, 'M3', 'NC_DOUBLE', c('sample','N'))
        var.def.nc(nc, 'M3stats', 'NC_DOUBLE', c('stats','N'))
        var.def.nc(nc, 'lnP', 'NC_DOUBLE', 'sample')
    
    att.put.nc(nc, "NC_GLOBAL", "title", "NC_CHAR", "MC results")
    att.put.nc(nc, "NC_GLOBAL", "stats", "NC_CHAR", "-1: mean; 16,50,84: quantiles; -2: sd")
    
    var.put.nc(nc, 'N3', 1:nn3)
    var.put.nc(nc, 'N', 1:nn)
    var.put.nc(nc, 'sample', 1:nsamples)
    var.put.nc(nc, 'stats', c(-1,50,16,80,-2))

    var.put.nc(nc, 'F', samplesF)
    var.put.nc(nc, 'Fstats', sstatF)
    var.put.nc(nc, 'M1', samplesM1)
    var.put.nc(nc, 'M1stats', sstatM1)
    var.put.nc(nc, 'M2', samplesM2)
    var.put.nc(nc, 'M2stats', sstatM2)
    var.put.nc(nc, 'M3', samplesM3)
    var.put.nc(nc, 'M3stats', sstatM3)
    var.put.nc(nc, 'lnP', c(lpff))

        close.nc(nc)
        
        saveRDS(samplesF,file=paste0('_Fsamples',savename,'.rds'))

    }
    
    print('Plots for combined data')

    Mstatmax <- max(c(sstatM1,sstatM2,sstatM3))
    ##Msampmax <- max(c(sstatM1,sstatM2,sstatM3))
    
    pdff(paste0('_plotstatsM1',savename))
    matplot(x=bincentres,y=nn*t(sstatM1[-5,]),type='l',lty=c(1,1,1,1,1),lwd=c(3,3,1,1,1),pch='.',col=c('black',myblue,myyellow,myyellow),xlab='X',ylab='N M1')
    ##matlines(x=(0:nt)/n,y=n*c(ft),type='l',lty=3,lwd=5,col=mygreen)
    legend('top',paste0(mylegend,'\nblack: mean, blue: median\nyellow: 16%q'),bty='n')
    dev.off()
 
    pdff(paste0('_plotstatsM2',savename))
    matplot(x=bincentres,y=nn*t(sstatM2[-5,]),type='l',lty=c(1,1,1,1,1),lwd=c(3,3,1,1,1),pch='.',col=c('black',myblue,myyellow,myyellow),xlab='X',ylab='N M2')
    ##matlines(x=(0:nt)/n,y=n*c(ft),type='l',lty=3,lwd=5,col=mygreen)
    legend('top',paste0(mylegend,'\nblack: mean, blue: median\nyellow: 16%q'),bty='n')
    dev.off()

    pdff(paste0('_plotstatsM3',savename))
    matplot(x=bincentres,y=nn*t(sstatM3[-5,]),type='l',lty=c(1,1,1,1,1),lwd=c(3,3,1,1,1),pch='.',col=c('black',myblue,myyellow,myyellow),xlab='X',ylab='N M3')
    ##matlines(x=(0:nt)/n,y=n*c(ft),type='l',lty=3,lwd=5,col=mygreen)
    legend('top',paste0(mylegend,'\nblack: mean, blue: median\nyellow: 16%q'),bty='n')
    dev.off()

    pdff(paste0('_plottraj',savename))
    matplot(x=1:nsamples,y=lpff,type='l',lty=c(1,2),lwd=1,pch='.',col=c(mygreen),xlab='t',ylab='logP')
    for(i in 1:ncol(samplesM1)){
matplot(x=1:nsamples,y=samplesM1[,i],type='l',lty=1,lwd=1,pch='.',col=c('black'),xlab='t',ylab=paste0('M1_',i))
    }
    for(i in 1:ncol(samplesM2)){
matplot(x=1:nsamples,y=samplesM2[,i],type='l',lty=1,lwd=1,pch='.',col=c('black'),xlab='t',ylab=paste0('M2_',i))
    }
    for(i in 1:ncol(samplesM3)){
matplot(x=1:nsamples,y=samplesM3[,i],type='l',lty=1,lwd=1,pch='.',col=c('black'),xlab='t',ylab=paste0('M3_',i))
    }
    dev.off()

    nexamples1 <- min(10,nsamples)
    nexamples2 <- min(100,nsamples)
    nexamples3 <- min(1000,nsamples)
    grows <- 10
    gcols <- 10
    minF <- 0
    samplesindex <- round(seq(1,nsamples,length.out=grows*gcols))
    samplesindex1 <- round(seq(1,nsamples,length.out=nexamples1))
    samplesindex2 <- round(seq(1,nsamples,length.out=nexamples2))
    samplesindex3 <- round(seq(1,nsamples,length.out=nexamples3))
    totake <- unique(c(samplesindex,samplesindex1,samplesindex2,samplesindex3))
    maxM <- max(nn*c(samplesM1[totake,],samplesM2[totake,],samplesM3[totake,]),na.rm=TRUE)
    gapx <- nn+round(nn/5)
    gapy <- maxM+maxM/5
    allmaxM <- (maxM+gapy)*grows

    pdff(paste0('_plotsamplesM1',savename))
    
    matplot(x=-(1:2),y=-(1:2),type='l',lty=1,col='white',xlab='',ylab='',xlim=c(0,gapx*gcols),ylim=c(0,gapy*grows))
    legend('top',mylegend,bty='n')
    k <- 0
    for(i in (1:gcols)-1){for(j in (1:grows)-1){k <- k+1
                             matlines(x=(1:nn)+i*gapx,y=nn*samplesM1[samplesindex[k],]+j*gapy,type='l',lty=1,col='black',xlab='',ylab='',xlim=c(0,gapx*gcols),ylim=c(0,gapy*grows))
##                             matlines(x=(0:nnt)+i*gapx,y=nn*c(Fmax)+j*gapy,type='l',lty=2,col=mygreen,xlab='',ylab='',xlim=c(0,gapx*gcols),ylim=c(0,gapy*grows))
                         }}
    matplot(x=bincentres,y=nn*t(samplesM1[samplesindex1,]),type='l',lty=1,lwd=1,pch='.',xlim=c(0,1),ylim=c(0,maxM),col=mygrey,xlab='N',ylab=paste0(nexamples1,' M1'))
    legend('top',mylegend,bty='n')

    matplot(x=bincentres,y=nn*t(samplesM1[samplesindex2,]),type='l',lty=1,lwd=1,pch='.',xlim=c(0,1),ylim=c(0,maxM),col=mygrey,xlab='N',ylab=paste0(nexamples2,' M1'))
    legend('top',mylegend,bty='n')

    matplot(x=bincentres,y=nn*t(samplesM1[samplesindex3,]),type='l',lty=1,lwd=1,pch='.',xlim=c(0,1),ylim=c(0,maxM),col=mygrey,xlab='N',ylab=paste0(nexamples3,' M1'))
    legend('top',mylegend,bty='n')

    dev.off()

        pdff(paste0('_plotsamplesM2',savename))
    
    matplot(x=-(1:2),y=-(1:2),type='l',lty=1,col='white',xlab='',ylab='',xlim=c(0,gapx*gcols),ylim=c(0,gapy*grows))
    legend('top',mylegend,bty='n')
    k <- 0
    for(i in (1:gcols)-1){for(j in (1:grows)-1){k <- k+1
                             matlines(x=(1:nn)+i*gapx,y=nn*samplesM2[samplesindex[k],]+j*gapy,type='l',lty=1,col='black',xlab='',ylab='',xlim=c(0,gapx*gcols),ylim=c(0,gapy*grows))
##                             matlines(x=(0:nnt)+i*gapx,y=nn*c(Fmax)+j*gapy,type='l',lty=2,col=mygreen,xlab='',ylab='',xlim=c(0,gapx*gcols),ylim=c(0,gapy*grows))
                         }}
    matplot(x=bincentres,y=nn*t(samplesM2[samplesindex1,]),type='l',lty=1,lwd=1,pch='.',xlim=c(0,1),ylim=c(0,maxM),col=mygrey,xlab='N',ylab=paste0(nexamples1,' M2'))
    legend('top',mylegend,bty='n')

    matplot(x=bincentres,y=nn*t(samplesM2[samplesindex2,]),type='l',lty=1,lwd=1,pch='.',xlim=c(0,1),ylim=c(0,maxM),col=mygrey,xlab='N',ylab=paste0(nexamples2,' M2'))
    legend('top',mylegend,bty='n')

    matplot(x=bincentres,y=nn*t(samplesM2[samplesindex3,]),type='l',lty=1,lwd=1,pch='.',xlim=c(0,1),ylim=c(0,maxM),col=mygrey,xlab='N',ylab=paste0(nexamples3,' M2'))
    legend('top',mylegend,bty='n')

    dev.off()

        pdff(paste0('_plotsamplesM3',savename))
    
    matplot(x=-(1:2),y=-(1:2),type='l',lty=1,col='white',xlab='',ylab='',xlim=c(0,gapx*gcols),ylim=c(0,gapy*grows))
    legend('top',mylegend,bty='n')
    k <- 0
    for(i in (1:gcols)-1){for(j in (1:grows)-1){k <- k+1
                             matlines(x=(1:nn)+i*gapx,y=nn*samplesM3[samplesindex[k],]+j*gapy,type='l',lty=1,col='black',xlab='',ylab='',xlim=c(0,gapx*gcols),ylim=c(0,gapy*grows))
##                             matlines(x=(0:nnt)+i*gapx,y=nn*c(Fmax)+j*gapy,type='l',lty=2,col=mygreen,xlab='',ylab='',xlim=c(0,gapx*gcols),ylim=c(0,gapy*grows))
                         }}
    matplot(x=bincentres,y=nn*t(samplesM3[samplesindex1,]),type='l',lty=1,lwd=1,pch='.',xlim=c(0,1),ylim=c(0,maxM),col=mygrey,xlab='N',ylab=paste0(nexamples1,' M3'))
    legend('top',mylegend,bty='n')

    matplot(x=bincentres,y=nn*t(samplesM3[samplesindex2,]),type='l',lty=1,lwd=1,pch='.',xlim=c(0,1),ylim=c(0,maxM),col=mygrey,xlab='N',ylab=paste0(nexamples2,' M3'))
    legend('top',mylegend,bty='n')

    matplot(x=bincentres,y=nn*t(samplesM3[samplesindex3,]),type='l',lty=1,lwd=1,pch='.',xlim=c(0,1),ylim=c(0,maxM),col=mygrey,xlab='N',ylab=paste0(nexamples3,' M3'))
    legend('top',mylegend,bty='n')

    dev.off()


    print('Diagnostics')

    ##pdff(paste0('_diagnostics',savefilesid))
    diagfile <- paste0('_diagnostics',savename,'.txt')
    capture.output(print(Consort(mcoutputx)),file=diagfile,append=T)

    capture.output(print('======================================================================'),file=diagfile,append=T)
    capture.output(print('----- Posterior checks -----'),file=diagfile,append=T)
    capture.output(print(PosteriorChecks(mcoutputx)$Posterior.Summary),file=diagfile,append=T)

capture.output(print('======================================================================'),file=diagfile,append=T)
capture.output(print('----- Hellinger-distance -----'),file=diagfile,append=T)
capture.output(print(t(BMK.Diagnostic(tocheck,batches=2))),file=diagfile,append=T)
capture.output(print('Hellinger-distance > 0.1'),file=diagfile,append=T)
capture.output(print(c(bmsum<-sum(BMK.Diagnostic(tocheck,batches=2)>0.1), bmsum/(nparm+1))),file=diagfile,append=T)
capture.output(print('Hellinger-distance > 0.05'),file=diagfile,append=T)
capture.output(print(c(bmsum<-sum(BMK.Diagnostic(tocheck,batches=2)>0.05), bmsum/(nparm+1))),file=diagfile,append=T)

capture.output(print('======================================================================'),file=diagfile,append=T)
capture.output(print('----- Effective sample size -----'),file=diagfile,append=T)
capture.output(print(ESS(tocheck)),file=diagfile,append=T)
   
}

#######################################################################
######################## 1-CHAIN DATA & PLOTS ########################
#######################################################################

if(singleplots){plotsanddiagn(mcoutputx=mcoutput,savename=savefilesid,savenc=FALSE)}

#######################################################################
######################## COMBINED DATA & PLOTS ########################
#######################################################################

if(jobid==1 && (combinedplots || !adaptphase)){
    
    combinedfilesid <- paste0(stagename,filesid)
    
    ## count how many chains we're running
    nchains <- length(list.files(pattern=paste0('^_LDlog[0-9]*',combinedfilesid,'\\.txt')))
    
    print('Waiting for all chains to finish...')
    while(length(list.files(pattern=paste0('^_output[0-9]*',combinedfilesid))) < nchains){Sys.sleep(10)}
    print('...done waiting.')
    if(!onlyplots){Sys.sleep(30)}
    
    print('Reading data from all chains...')
    ## read Monte Carlo outputs ##
    mcoutputlist <- foreach(jid=1:nchains)%do%{
        readRDS(paste0('_output',jid,combinedfilesid,'.rds'))}

    plotsanddiagn(mcoutputx=Combine(mcoutputlist,mydata),savename=paste0('_ALL',combinedfilesid),savenc=TRUE)

}

print('...Finished.')
