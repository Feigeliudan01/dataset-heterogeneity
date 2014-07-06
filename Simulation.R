## Import: input directory .../prepared/ ONLY contains the prepared 
## and the prepared_rm_ RData files of the following data sets:
## CAL, EXPO, MAINZ, MDA4, MSK, SUPERTAM_HGU133A, SUPERTAM_HGU133PLUS2, TRANSBIG, UNT, VDX 
## XYzranking, seqoverlap and av_vars RData are also in the directory.
## That adds up to 23 files totally in the directory.
## The variable 'seq' should be altered if there are other files.

library(WGCNA)  
library(survival)  
library(pls)
library(CoxBoost)

## the matrix is saved in directory Cmatrix
inputdir="prepared/"
transitdir="samples/" # to store 150 samples from each data sets
outputdir="Cmatrix_simulation/"


## Load from directory .../prepared/ 
load(paste(inputdir, "XYzranking.RData", sep=""))
load(paste(inputdir, "seqoverlap.RData", sep=""))



seq <- 1:10  # which datasets we use
N <- length(seq) # number of datasets
C <- matrix(rep(0,times=N*N), nrow=N, ncol=N) # matrix of C-statistics
Cind <- matrix(rep(0,times=N*N), nrow=N, ncol=N)
P <- length(genesymb)  # number of genes
B <- 2  ## total interactions
fold <- 4 
M <- 10 ## number of boosting steps
penal <- 10 ## penalty in CoxBoost
numofSamples <- 150


## function of generating survival time using CoxBoost method
generate_survTime <- function(par_Y, par_X, par_step, par_penal, par_n){
  cbfit <- CoxBoost(time=par_Y[,1], status=par_Y[,2], x=par_X, stepno=par_step, penalty=par_penal) 
  beta <- coef(cbfit)
  survTime <- numeric(par_n)
  time <- par_Y[, 1]
  print(min(time))
  status <- par_Y[, 2]
  data <- data.frame(time, status)
  SFit <- survfit(Surv(time, status) ~ 1, data = data)
  nelsonaalenEstimate <- cumsum(SFit$n.event / SFit$n.risk)
  for(i in 1:par_n){
    u <- runif(1, min = 0, max = 1)
    par_a <- log(u,base=exp(1))
    par_b <- t(as.matrix(beta))
    par_c <- as.matrix(scale(par_X))[i,]
    max_t <- max(nelsonaalenEstimate[nelsonaalenEstimate <= as.numeric((-par_a) * (exp(-(par_b %*% par_c ))))])
    id <- which(nelsonaalenEstimate == max_t)
    TIME <- max(time[id], min(time[!is.na(time)])) 
    ## ?? Is this right??
    ## ?? How to get the inverse of the Nelsin-Aalen estimator function?
    print(paste(i, TIME))
    survTime[i] <- TIME
  }  
  return(survTime)
}


func_train <- function(par_b, par_i, par_setsID){  
  print(paste("TRAIN: interaction", par_b, "- set", par_setsID[par_i], sep=" "))
  load(paste(transitdir, "interaction_", b, "_", setsID[i], ".RData", sep=""))
  Xi <- scale(newX)
  Yi <- dmfsYL_new

  ## univariate Cox regression on gene level
  alpha <- v <- numeric(P) 
  for(p in 1:P){
    #print(p)
    test <- list(time=Yi[,1], status=Yi[,2], x=Xi[, p])
    funcCSV <- coxph(Yi ~ Xi[, p], test)
    ## testing proportional hazard assumption
    plot(cox.zph(funcCSV))
    alpha[p] <- funcCSV$coefficients
  
    ## Calculate vp
    if(alpha[p] > 0) v[p] <- 1
    if(alpha[p] == 0) v[p] <- 0
    if(alpha[p] < 0) v[p] <- -1
    v[p] <- v[p] / sqrt(P)
    rm(newX)
  }
  return(v)
}


func_test <- function(par_b, par_j, par_setsID, par_v){
  print(paste("TEST: interaction", par_b, "- set", par_setsID[par_j], sep=" "))
  load(paste(transitdir, "interaction_", b, "_", par_setsID[par_j], ".RData", sep=""))
  Xj <- scale(newX)
  Yj <- dmfsYL_new
  Sj <- Xj %*% par_v 
  rm(newX)
  return(Sj)
}


calcu_cstat <- function(par_numTest, par_Yj, par_Sj){
  count <- 0 
  count1 <- 0 
  for(a in 1:par_numTest){
    for(d in 1:par_numTest){
      if((!is.na(par_Yj[a,1])) & (!is.na(par_Yj[d,1]))) {
        if(par_Yj[a,1] < par_Yj[d,1]) count <- count + 1
        if(par_Yj[a,1] < par_Yj[d,1] & par_Sj[a] > par_Sj[d]) count1 <- count1 + 1
      }
    }
  }
  if(count > 0 & count1 > 0) result <- count1 / count
  else result <- 0  
  return(result)
}


## Simulation
for(b in 1:B){  
  print(paste("interaction = ", b, sep=""))
  ## Draw from data sets with replacement
  print("drawing sets")
  prob <- rep((1/(N-2)),times=N-2)
  setsID <- sample(seq[-c(2, 4)], prob=prob, replace=TRUE)
  print(setsID)

  
  ## Draw samples from set with replacement
  ## Generating survival time
  print("drawing samples")
  uniqueID <- unique(setsID)
  uniqueN <- length(uniqueID)
  for(i in 1:uniqueN){
    print(uniqueID[i])
    load(paste(inputdir, "prepared_rm_", names[uniqueID[i]], ".RData", sep=""))
    n <- nrow(newX)
    sampleind <- sample(1:n, numofSamples, replace=TRUE)
    newX <- newX[sampleind, ]
    dmfsYL_new <- dmfsYL[[uniqueID[i]]][sampleind, ]
    ## fill in missing value
    dmfsYL_new[which(is.na(dmfsYL_new[,1])), 1] =  mean(dmfsYL_new[,1])
    dmfsYL_new[which(is.na(dmfsYL_new[,2])), 2] =  0
    survTime <- generate_survTime(dmfsYL_new, newX, M, penal, numofSamples)   
    dmfsYL_new <- Surv(survTime, dmfsYL_new[,2])        
    filepath <- paste(transitdir, "interaction_", b, "_", uniqueID[i], ".RData", sep="")
    save(newX, dmfsYL_new, file=filepath)
    rm(newX)
  }
  
  
  ## Validation process
  for(i in 1:(N-2)){
    for(j in 1:(N-2)){
      if(setsID[i] != setsID[j]){ # train on setsID[i], validate on setsID[j]
        v <- func_train(par_b=b, par_i=i, par_setsID=setsID)
        Sj <- func_test(par_b=b, par_j=j, par_setsID=setsID, par_v=v)
        load(paste(transitdir, "interaction_", b, "_", setsID[j], ".RData", sep=""))
        Yj <- dmfsYL_new
        numTest <- numofSamples  ### number of samples in j
        result <- calcu_cstat(par_numTest=numTest, par_Yj=Yj, par_Sj=Sj)
        C[setsID[i], setsID[j]] <- C[setsID[i], setsID[j]] + result                              
        Cind[setsID[i], setsID[j]] <- Cind[setsID[i], setsID[j]] + 1
      }      
      
      else if(setsID[i] == setsID[j]){
        print(paste("CV: interaction:", b, "set:", setsID[i], sep=" "))
        load(paste(transitdir, "interaction_", b, "_", setsID[i], ".RData", sep=""))
        sizeofTest <- round(numofSamples / fold) # size of the testing set in data set ii
        
        ##  Randomly separate the samples into four groups
        ind <- 1:numofSamples
        setind1 <- sample(ind, sizeofTest, replace=FALSE) 
        setind2 <- sample(ind[-setind1], sizeofTest, replace=FALSE)
        setind3 <- sample(ind[-c(setind1,setind2)], sizeofTest, replace=FALSE)
        setind4 <- ind[-c(setind1,setind2,setind3)]
        dividedSets <- list(Set1=(scale(newX))[setind1, ], Set2=(scale(newX))[setind2, ],
                            Set3=(scale(newX))[setind3, ], Set4=(scale(newX))[setind4, ])
        dividedResponse <- list(Res1=dmfsYL_new[setind1, ], Res2=dmfsYL_new[setind2, ],
                                Res3=dmfsYL_new[setind3, ], Res4=dmfsYL_new[setind4, ])
        
        cstat <- numeric(fold)
        for(k in 1:fold){
          print(paste("fold =", k, sep=" "))          
          testX <- dividedSets[[k]]
          testY <- dividedResponse[[k]]
          trainX <- rbind(dividedSets[-k][[1]], dividedSets[-k][[2]], dividedSets[-k][[3]])
          trainY <- rbind(dividedResponse[-k][[1]], dividedResponse[-k][[2]], dividedResponse[-k][[3]])
                    
          alpha <- v <- numeric(P)   
          for(p in 1:P){
            #print(p)
            test <- list(time=trainY[,1], status=trainY[,2], x=trainX[, p])
            try(funcCV <- coxph(Surv(trainY[,1],trainY[,2]) ~ trainX[, p], test))
            ## testing proportional hazard assumption
            plot(cox.zph(funcCV))
            alpha[p] <- funcCV$coefficients
            ## Calculate vp
            if(alpha[p] > 0) v[p] <- 1
            if(alpha[p] == 0) v[p] <- 0
            if(alpha[p] < 0) v[p] <- -1
            v[p] <- v[p] / sqrt(P)
          }
          scoreY <- testX %*% v
          numTest <- nrow(testX)
          cstat[k] <- calcu_cstat(par_numTest=numTest, par_Yj=testY, par_Sj=scoreY)
        }       
        C[setsID[i], setsID[i]] <- C[setsID[i], setsID[i]] + sum(cstat) / fold
        Cind[setsID[i], setsID[i]] <- Cind[setsID[i], setsID[i]] + 1
      }
    }
  }
}

C <- C / Cind


print(C)
## save C matrix for later use in the output directory
save(C, file=paste(outputdir, "C_matrix.RData", sep=""))
