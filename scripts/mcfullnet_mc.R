## Monte Carlo calculation of posterior for joint frequency distribution
## Parameters are in mcfullnet_params.R
## This script is called with
##     sbatch core5ar.slurm mcfullnet_mc.R
## The output must then be combined with
##     Rscript -e "source('mcfullnet_combine.R')" '[inputfile]'
## for example
##     Rscript -e "source('mcfullnet_combine.R')" 'mcs_nuR_afssblock3Snu_adapt2_N5000_t250_L1_T417641_sm5804426401.41682_s7000'
##

## prob. distribution for joint frequencies R and long-run frequencies nu
jobid <- as.numeric(commandArgs(trailingOnly=TRUE))
set.seed(181225+jobid)

#### Read parameters ####
source('mcfullnet_params.R')


## copy script for later checks
if(jobid==1){
    file.copy(paste0(scriptname,'_params.R'),paste0('_sc_params',fileid,'.R'),overwrite=TRUE)
    file.copy(paste0(scriptname,'_mc.R'),paste0('_sc_mc',fileid,'.R'),overwrite=TRUE)
    file.copy(paste0(scriptname,'_combine.R'),paste0('_sc_combine',fileid,'.R'),overwrite=TRUE)
}

#### Sampling ####

#stopCluster(cl)
    ## cl <- makeCluster(clusters, outfile="")
    ## registerDoParallel(cl)


##Initial.Values <- GIV(logprob, mydata, n=10, PGF=TRUE) # fails in parallel for nnt=165 for some reason...
Initial.Values <- PGF(mydata)

print('running Monte Carlo...')
mcoutput <- LaplacesDemon(logprob, mydata, Initial.Values,
                      ##Covar=NULL,
                      Covar=covm,
                      Thinning=thinning,
                      Iterations=mciterations, Status=mcstatus,
                      Debug=list(DB.chol=FALSE, DB.eigen=FALSE,
                                 DB.MCSE=FALSE, DB.Model=FALSE),
                      LogFile=paste0('_LDlog',jobid,fileid), #Type="MPI", #Packages=c('MASS'),
                      ##Algorithm="RDMH"#, Specs=list(B=list(1:d,d1:d2,d3:dnp))
                      ##Algorithm="Slice", Specs=list(B=NULL, Bounds=c(0,1), m=Inf, Type="Continuous", w=0.001)
                      ##Algorithm="MALA", Specs=list(A=1e7, alpha.star=0.574, gamma=mcadapt, delta=1, epsilon=c(1e-6,1e-7))
                      ##Algorithm="pCN", Specs=list(beta=0.001)
                      Algorithm="AFSS", Specs=list(A=mcadapt, B=B, m=nsteps, n=nprev, w=stepwidth)
                      ##Algorithm="AIES", Specs=list(Nc=4*nparm, Z=NULL, beta=2, CPUs=1, Packages=NULL, Dyn.libs=NULL)
                      
                      ##Algorithm="DRM", Specs=NULL
                      )
    
    ## RWM:no HARM:no pCN:no
    ## mcoutput <- Combine(mcoutputp,mydata)

#    cl <- makeCluster(clusters, outfile='/dev/null')
#     registerDoParallel(cl)

    
    ## save all output for later computation
    saveRDS(mcoutput,file=paste0('_mcoutput',jobid,fileid,'.rds'))

print('...MC finished. Consorting')

print(Consort(mcoutput))
