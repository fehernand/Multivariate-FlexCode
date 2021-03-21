library(glmnet)
library(SAM)
library(fields)
library(SDMTools)
library(randomForest)
#source("../../CMU/Astro/specSeriesInference_4/R/auxFunctions.r")
source("seriesRegression.R")


calculateBasis=function(z,nZ)
{
  # z is a n X dim vector
  # nZ is the number of basis of each coordinate
  nEval=dim(z)[1]
  nZOld=nZ
  if(nZ<3)
  {
    nZ=3
  }
  
  zAux=z[,1,drop=T]
  sinBasisZ=apply(as.matrix(1:round((nZ)/2)),1,function(xx) sqrt(2)*sin(2*xx*pi*zAux))
  cosBasisZ=apply(as.matrix(1:round((nZ)/2)),1,function(xx) sqrt(2)*cos(2*xx*pi*zAux))
  basisZ1=matrix(NA,length(zAux),2*round((nZ)/2))
  basisZ1[,seq(from=1,length.out=dim(sinBasisZ)[2],by=2)]=sinBasisZ
  basisZ1[,seq(from=2,length.out=dim(cosBasisZ)[2],by=2)]=cosBasisZ
  basisZ1=cbind(rep(1,length(zAux)),basisZ1)
  basisZ1=basisZ1[,1:nZOld,drop=F]
  
  zAux=z[,2,drop=T]
  sinBasisZ=apply(as.matrix(1:round((nZ)/2)),1,function(xx) sqrt(2)*sin(2*xx*pi*zAux))
  cosBasisZ=apply(as.matrix(1:round((nZ)/2)),1,function(xx) sqrt(2)*cos(2*xx*pi*zAux))
  basisZ2=matrix(NA,length(zAux),2*round((nZ)/2))
  basisZ2[,seq(from=1,length.out=dim(sinBasisZ)[2],by=2)]=sinBasisZ
  basisZ2[,seq(from=2,length.out=dim(cosBasisZ)[2],by=2)]=cosBasisZ
  basisZ2=cbind(rep(1,length(zAux)),basisZ2)
  basisZ2=basisZ2[,1:nZOld,drop=F]
  
  basisZ=matrix(t(apply(basisZ1,2,function(x)as.vector(t(x*basisZ2)))),nrow(basisZ1),ncol(basisZ1)*ncol(basisZ2),byrow = T)
  
  return(basisZ)
  
}



fitCDE=function(xTrain,zTrain,nIMax=length(zTrain),regressionFunction,regressionFunction.extra=NULL)
{
  objectCDE=NULL
  objectCDE$zMax=apply(zTrain,2,max)
  objectCDE$zMin=apply(zTrain,2,min)
  zTrain[,1]=(zTrain[,1]-objectCDE$zMin[1])/(objectCDE$zMax[1]-objectCDE$zMin[1])
  zTrain[,2]=(zTrain[,2]-objectCDE$zMin[2])/(objectCDE$zMax[2]-objectCDE$zMin[2])
  
  responseFourier=calculateBasis(zTrain,nIMax)
  regressionObject=regressionFunction(x=xTrain,responses=responseFourier,extra=regressionFunction.extra)
  objectCDE$nIMax=nIMax
  objectCDE$zTrain=zTrain
  objectCDE$xTrain=xTrain
  objectCDE$regressionObject=regressionObject
  rm(regressionObject,xTrain,zTrain,responseFourier)
  gc(verbose = FALSE)
  return(objectCDE)
}


fitCDE.tune=function(objectCDE,xValidation,zValidation)
{
  zValidation[,1]=(zValidation[,1]-objectCDE$zMin[1])/(objectCDE$zMax[1]-objectCDE$zMin[1])
  zValidation[,2]=(zValidation[,2]-objectCDE$zMin[2])/(objectCDE$zMax[2]-objectCDE$zMin[2])
  
  basisZValidation=calculateBasis(zValidation,objectCDE$nIMax) # returns matrix length(z)xnIMax with the basis for z
  
  coefficientsXValidation=predict(objectCDE$regressionObject,xValidation)
  
  term1=1/2*colMeans(coefficientsXValidation^2)
  
  term1SubSums=t(apply(as.matrix(1:objectCDE$nIMax),1,function(x)
  {
    sequenceX=outer(1:x ,seq(0,objectCDE$nIMax^2-2,by=objectCDE$nIMax),"+")
    subSums=apply(sequenceX,2,function(y)sum(term1[y]))
    
    cumsum(subSums)
  }))
  
  term2=colMeans(coefficientsXValidation*basisZValidation[,1:ncol(coefficientsXValidation),drop=F])
  term2SubSums=t(apply(as.matrix(1:objectCDE$nIMax),1,function(x)
  {
    sequenceX=outer(1:x ,seq(0,objectCDE$nIMax^2-2,by=objectCDE$nIMax),"+")
    subSums=apply(sequenceX,2,function(y)sum(term2[y]))
    
    cumsum(subSums)
  }))
  aux=term1SubSums-term2SubSums
  bestI=which(aux==min(aux),arr.ind = T)[1,]
  
  objectCDE$errors=aux
  objectCDE$bestI=bestI # first variable, second variable
  objectCDE$bestError=min(objectCDE$errors)
  return(objectCDE)
  
}

predictCDE=function(objectCDE,xNew,B=100)
{
  zGrid=seq(from=0,to=1,length.out=B)
  zGrid=expand.grid(zGrid,zGrid)
  
  binSize=((1)/(B+1))^2
  delta=ifelse(!is.null(objectCDE$bestDelta),objectCDE$bestDelta,0)
  
  
  coeff=predict(objectCDE$regressionObject,xNew,maxTerms=objectCDE$nIMax^2)
  basisZNew=calculateBasis(zGrid,objectCDE$nIMax) # returns matrix length(z)xnIMax with the basis for z
  sequenceX=outer(1:objectCDE$bestI[1] ,seq(0,objectCDE$nIMax*objectCDE$bestI[2]-1,by=objectCDE$nIMax),"+")
  estimates=coeff[,as.vector(sequenceX)]%*%t(basisZNew[,as.vector(sequenceX)])
  
  
  #estimates=t(apply(estimates,1,function(xx).normalizeDensity(xx,delta)))
  
  estimates=estimates/(sum(estimates)*binSize)
  
  estimates=estimates/((objectCDE$zMax[1]-objectCDE$zMin[1])*(objectCDE$zMax[2]-objectCDE$zMin[2]))
  
  returnValue=NULL
  returnValue$CDE=estimates
  returnValue$z1=seq(from=objectCDE$zMin[1],to=objectCDE$zMax[1],length.out=B)
  returnValue$z2=seq(from=objectCDE$zMin[2],to=objectCDE$zMax[2],length.out=B)
  return(returnValue)
}

.normalizeDensity=function(estimates,delta=0)
{
  # internal function of renormalization of density
  estimatesThresh=estimates
  estimatesThresh[estimatesThresh<1e-8]=0
  
  m=matrix(estimatesThresh,sqrt(length(estimatesThresh)),sqrt(length(estimatesThresh)))
  connectedComp=ConnCompLabel(m>0)
  
  areas=tapply(as.vector(m), as.vector(connectedComp), sum)
  areas=areas/sum(areas)
  
  whichRemove=as.numeric(names(which(areas<delta)))
  m[connectedComp%in%whichRemove]=0
  
  return(as.vector(m))
  
}

chooseDelta = function(objectCDE, xValidation,zValidation,deltaGrid=seq(0,0.4,0.05)) 
{
  error=rep(NA,length(deltaGrid))
  cat("\n Progress Bar:\n")
  for(ii in 1:length(deltaGrid))
  {
    cat(paste(c(rep("|",ii),rep(" ",length(deltaGrid)-ii),"|\n"),collapse=""))
    objectCDE$bestDelta=deltaGrid[ii]
    estimateErrors=estimateErrorCDE(objectCDE=objectCDE,xTest=xValidation,zTest=zValidation,se=FALSE)
    error[ii]=estimateErrors
  }
  whichMin=(1:length(error))[error==min(error)]
  bestDelta=deltaGrid[max(whichMin)]
  return(bestDelta)
}



estimateErrorCDE=function(objectCDE=objectCDE,xTest,zTest,se=FALSE)
{
  zGrid1=seq(objectCDE$zMin[1],objectCDE$zMax[1],length.out=50)
  zGrid2=seq(objectCDE$zMin[2],objectCDE$zMax[2],length.out=50)
  
  predictedComplete=predictCDE(objectCDE,xNew = xTest,B=length(zGrid1))
  predictedComplete=predictedComplete$CDE*(objectCDE$zMax[1]-objectCDE$zMin[1])*(objectCDE$zMax[2]-objectCDE$zMin[2])
  
  colmeansComplete=colMeans(predictedComplete^2)
  sSquare=mean(colmeansComplete)
  
  zGrid=expand.grid(zGrid1,zGrid2)
  n=nrow(zTest)
  predictedObserved=apply(as.matrix(1:n),1,function(xx) { index=which.min(abs(rdist(zTest[xx,],zGrid)))
  return(predictedComplete[xx,index])
  })
  likeli=mean(predictedObserved)
  gc()
  
  return(1/2*sSquare-likeli)
  
  
}

regressionFunction.SAM=function(x,responses,extra=NULL)
{
  # Both x and responses are matrices
  n=dim(x)[1]
  random=sample(1:n)
  nTrain=round(0.7*n)
  xTrain=x[random[1:nTrain],]
  responsesTrain=responses[random[1:nTrain],]
  xValidation=x[random[-c(1:nTrain)],,drop=FALSE]
  responsesValidation=responses[random[-c(1:nTrain)],,drop=FALSE]
  
  
  sVec=extra$sVec
  if(is.null(sVec))
    sVec=round(seq(1,14,length.out = 6))
  
  
  
  fittedReg=apply(as.matrix(1:ncol(responsesTrain)),1,function(ii){
    
    bestCol=bestError=rep(NA,length(sVec))
    for(s in 1:length(bestCol))
    {
      fit=try(samQL(xTrain,responsesTrain[,ii,drop=FALSE],p = sVec[s]),silent = TRUE)
      if(class(fit)=="try-error")
      {
        bestError[s]=Inf
        next;
      }
      out = predict(fit,xValidation)    
      out$values[is.na(out$values)]=0
      errorPerLambda=colMeans((out[[1]]-matrix(responsesValidation[,ii,drop=FALSE],nrow(out[[1]]),ncol(out[[1]])))^2)
      bestCol[s]=which.min(errorPerLambda)
      bestError[s]=min(errorPerLambda)
    }
    bestS=sVec[which.min(bestError)]
    fit=samQL(x,responses[,ii,drop=FALSE],p = bestS)
    object=NULL
    object$fit=fit
    object$bestS=bestS
    object$bestCol=bestCol[which.min(bestError)]
    object$whichCovariates=fit$func_norm[,object$bestCol]>1e-23
    return(object)
  })
  
  regressionObject=NULL
  regressionObject$fittedReg=fittedReg
  class(regressionObject)="SAM"
  rm(fittedReg,responsesTrain,xTrain,xValidation,responsesValidation)
  gc(verbose=FALSE)
  return(regressionObject)
}

predict.SAM=function(object,xNew,maxTerms=NULL)
{
  if(class(object)!="SAM")
    stop("Object has wrong class, should be SAM")
  
  if(!is.null(maxTerms))
  {
    maxTerms=min(maxTerms,length(object$fittedReg))
  } else {
    maxTerms=length(object$fittedReg)
  }
  
  predictedValidation=apply(as.matrix(1:maxTerms),1,function(xx)
  {
    out.tst = try(predict(object$fittedReg[[xx]]$fit,xNew),silent = TRUE)
    if(class(out.tst)=="try-error")
    {
      return(NA)
    }
    predicted=out.tst[[1]][,object$fittedReg[[xx]]$bestCol]
    return(predicted)
  })
  
  whichNa=sapply(predictedValidation,function(x)any(is.na(x)))
  if(all(!whichNa))
  {
    return(predictedValidation)
  }
  whichSelect=rep(F,length(predictedValidation))
  whichSelect[1:(which.max(whichNa)-1)]=T
  predictedValidation=subset(predictedValidation,subset = whichSelect)
  predictedValidation=sapply(predictedValidation,function(x)x)
  return(predictedValidation)
  
}

#### Functions

require(fields)
require(ggplot2)
require(scales)
require(grid)


radialKernel=function(distances,extraKernel=list("eps.val"=1))
{
  # Given the distances and the bandwidth eps.val, computes the matrix of radial kernel
  return(exp(-distances^2/(4*extraKernel$eps.val)))
}


seriesRegression=function(xTrain=NULL,xValidation=NULL,kernelTrainTrain=NULL,kernelValidationTrain=NULL,
                          zTrain,zValidation,nXMax=100,methodEigen=ifelse(length(zTrain)>500,"SVDRandom","SVD"),
                          kernelFunction=ifelse(is.null(kernelTrainTrain),radialKernel,NA),extraKernel=NULL,basis="KPCA",verbose=TRUE)
{
  # main function
  if(!is.matrix(zTrain))
    zTrain=matrix(zTrain,length(zTrain),1)
  
  if(!is.matrix(zValidation))
    zValidation=matrix(zValidation,length(zValidation),1)
  
  #zTrain=as.numeric(zTrain)
  #zValidation=as.numeric(zValidation)
  
  if(!is.null(kernelTrainTrain)|!is.null(kernelValidationTrain))
  {
    if(!is.null(xTrain)|!is.null(xTrain))
      stop("Please provide wither xTrain-type of arguments or kernelTrainTrain-type of arguments")
  }
  
  if(is.null(xTrain)&!is.null(xValidation)) stop("Please provide xTrain if you provide xValidation")
  
  if(!is.null(xTrain)&is.null(xValidation)) stop("Please provide xValidation if you provide xTrain")
  
  if(!is.null(xTrain)&!is.null(kernelTrainTrain)) stop("Please provide either xTrain of kernelTrainTrain")
  
  if(!is.null(xValidation)&!is.null(kernelValidationTrain)) stop("Please provide either xValidation of kernelValidationTrain")
  
  
  if(is.null(xTrain)&(is.null(kernelTrainTrain)|is.null(kernelTrainTrain))) stop("You have to provide either kernelTrainTrain and kernelValidationTrain or xTrain and xValidation")
  
  if(!is.null(xValidation))
  {
    if(!is.matrix(xValidation))
      stop("xValidation has to be a Matrix")
    
    if(!is.matrix(xTrain))
      stop("xTrain has to be a Matrix")
    
    
    if(dim(xValidation)[2]!=dim(xTrain)[2]) stop("Dimensions of x don't match!")
    
    if(dim(xTrain)[1]!=dim(zTrain)[1]) stop("Dimensions of x and z train don't match!")
    
    if(dim(xValidation)[1]!=dim(zValidation)[1]) stop("Dimensions of x and z validation don't match!")
  }
  
  if(basis!="KPCA" & basis!="Diffusion")
    stop("basis not implemented. Has to be KPCA or Diffusion")
  
  if(verbose) cat("\n Fit model")  
  
  if(!is.null(xTrain))
  {
    object=.seriesRegressionFit(x=xTrain,z=zTrain,nXMax=nXMax,methodEigen=methodEigen,kernelFunction=kernelFunction,extraKernel=extraKernel,basis=basis,userCall=FALSE,verbose=verbose)
  } else {
    if(!is.na(kernelFunction))
      warning("Argument's kernelFunction and extraKernel is being ignored: inputs given by the user are already Gram matrices")
    object=.seriesRegressionFit(kernelMatrix=kernelTrainTrain,z=zTrain,nXMax=nXMax,methodEigen=methodEigen,basis=basis,userCall=FALSE,verbose=verbose)
  }
  
  
  if(verbose) cat("\n Tune I")  
  
  if(!is.null(xTrain))
  {
    object=.seriesRegressionTuning(object=object,xValidation=xValidation,zValidation=zValidation,userCall=F)
  } else {
    object=.seriesRegressionTuning(object=object,kernelNewOld=kernelValidationTrain,zValidation=zValidation,userCall=F)
  }
  
  
  if(verbose) cat("\n Compute Validation set error")
  object$bestError=min(object$errors)
  
  
  if(abs(object$nXBest-object$nX)<5) warning("nXBest is close to nXMax, you may consider increasing nXMax")
  
  return(object)  
}


.seriesRegressionFit=function(x=NULL,kernelMatrix=NULL,z,nXMax=NULL,methodEigen,kernelFunction=radialKernel,extraKernel=NULL,basis=NULL,userCall=TRUE,verbose=TRUE)
{
  # x is n by d
  # z is n by 1, and its elements are assumed to be between 0 and 1
  # estimates f(z|X)
  # returns all coefficients up to nZ, the maximum number of components
  # for Z
  # zTrain between 0 and 1
  if(userCall)
    cat("\n A user will typically use function 'seriesRegression' instead, are you sure about what you are doing?")
  
  
  object=list()
  class(object)="seriesRegression"
  
  
  if(is.null(kernelMatrix)) # if kernel matrix is not provided
  {
    distancesX=rdist(x,x) # of all training samples to the selected ones
    object$xSelectedTrain=x
    rm(x)
    object$kernelFunction=kernelFunction
    
    if(identical(kernelFunction,radialKernel)&is.null(extraKernel)) # default choice of bandwidth
    {
      extraKernel=list("eps.val"=median(distancesX^2)/8)
    }
    kernelMatrix=kernelFunction(distancesX,extraKernel)
    object$extraKernel=extraKernel
    rm(distancesX)
  }
  
  if(any(is.na(kernelMatrix))) stop("Kernel with NA")
  
  nAll = length(kernelMatrix[,1]) # sample size of all training  
  
  
  if(verbose) cat("\n Computing EigenVectors for basis for X")  
  
  if(basis=="KPCA")
  {
    
    if(methodEigen=="SVD")
    {    
      results=eigen(kernelMatrix,symmetric=T)
      nX=dim(results$vectors)[2]
      
      if(!is.null(nXMax)) 
      {
        nXMax=min(nX,nXMax)
      } else {
        nXMax=nX
      }
      
      basisX=sqrt(nAll)*results$vectors[,1:nXMax,drop=FALSE]
      eigenValues=results$values[1:nXMax]/nAll
      
      rm(results)
      
    } else if(methodEigen=="SVDRandom") {
      if(!is.null(nXMax)) 
      {
        nXMax=min(length(z)-15,nXMax)
      } else {
        nXMax=length(z)-15
      }
      
      p=10
      
      Omega=matrix(rnorm(nAll*(nXMax+p),0,1),nAll,nXMax+p)
      Z=kernelMatrix%*%Omega
      Y=kernelMatrix%*%Z
      Q=qr(x=Y)
      Q=qr.Q(Q)
      B=t(Q)%*%Z%*%solve(t(Q)%*%Omega)
      eigenB=eigen(B)
      lambda=eigenB$values
      U=Q%*%eigenB$vectors
      
      basisX=Re(sqrt(nAll)*U[,1:nXMax])
      eigenValues=Re((lambda/nAll)[1:nXMax])
      
      
      rm(Omega,Z,Y,Q,B,U,eigenB)
      gc()
      
    } else   {
      stop("Method not implemented")
    }
    
    
    if(verbose) cat("\n Computing coefficients")  
    coefficients=1/nAll*t(z)%*%basisX
    
  } else { # diffusion basis
    
    stationary=rowSums(kernelMatrix)/sum(sum(kernelMatrix))
    kernelMatrixN=kernelMatrix/(sqrt(rowSums(kernelMatrix))%*%t(sqrt(colSums(kernelMatrix))))
    rm(kernelMatrix)
    gc()
    
    if(methodEigen=="SVD")
    {    
      eigenDec=eigen(kernelMatrixN,symmetric=T)
      nX=dim(eigenDec$vectors)[2]
      
      if(!is.null(nXMax)) 
      {
        nXMax=min(nX,nXMax)
      } else {
        nXMax=nX
      }
      
      resultsVectors=Re(sqrt(nAll)*eigenDec$vectors[,1:nXMax])
      resultsValues=Re((eigenDec$values/nAll)[1:nXMax])
      rm(eigenDec)
      gc()
      
    } else if(methodEigen=="SVDRandom") {
      if(!is.null(nXMax)) 
      {
        nXMax=min(length(z)-15,nXMax)
      } else {
        nXMax=length(z)-15
      }
      
      p=10
      Omega=matrix(rnorm(nAll*(nXMax+p),0,1),nAll,nXMax+p)
      Z=kernelMatrixN%*%Omega
      Y=kernelMatrixN%*%Z
      #Y=kernelMatrixN%*%Y
      #Y=kernelMatrixN%*%Y
      Q=qr(x=Y)
      Q=qr.Q(Q)
      B=t(Q)%*%Z%*%solve(t(Q)%*%Omega)
      eigenB=eigen(B)
      lambda=Re(eigenB$values)
      U=Q%*%Re(eigenB$vectors)
      resultsVectors=Re(sqrt(nAll)*Re(U[,1:nXMax]))
      resultsValues=Re((lambda/nAll)[1:nXMax])
      rm(Omega,Z,Y,Q,B,U,eigenB)
      gc()
      
    }
    
    stationary=stationary*nAll
    stationaryMatrix=matrix(rep(stationary,nXMax),nAll,nXMax,byrow=F)
    basisX=resultsVectors*1/sqrt(stationaryMatrix)
    
    matrixEigenDensity=basisX*stationaryMatrix
    
    if(verbose) cat("\n Computing coefficients")  
    coefficients=1/nAll*t(z)%*%matrixEigenDensity
    
    eigenValues=resultsValues
    
    rm(stationaryMatrix,resultsValues)
    gc()
  }
  
  gc(verbose=F)
  object$coefficients=coefficients
  object$methodEigen=methodEigen
  object$nX=nXMax
  object$Z=z
  object$eigenX=basisX
  object$eigenValuesX=eigenValues
  object$basis=basis
  
  return(object)  
  
}  


plot.seriesRegression=function(object,component1=1,component2=2)
{
  data=NULL
  data$x=object$eigenX[,component1]
  data$y=object$eigenX[,component2]
  data$z=object$Z
  data=as.data.frame(data)
  
  
  
  theme_set(theme_gray(base_size = 35))
  
  
  g2<- ggplot(data, aes(x = y, y = x)) + 
    geom_point(size=5, aes(color=z)) + 
    theme(panel.background = element_rect(fill='white', colour='grey1'))+
    theme(legend.position = "top" ) + xlab(bquote(psi[.(component1)]))+ylab(bquote(psi[.(component2)]))+
    #scale_color_gradient(name="Z",low="navyblue",high="red")+
    scale_colour_gradientn(colours=rainbow(5))+
    theme(legend.text = element_text(size = 20))+
    theme(legend.title = element_text(size = 24))+
    theme(legend.key.height=unit(1.5,"line")) +# Change 3 to X
    theme(legend.key.width=unit(2.2,"line"))+ # Change 3 to X
    theme(legend.justification=c(1,1), legend.position=c(1,1))+
    theme(legend.background = element_rect(fill=alpha('blue', 0)))+
    guides(color = guide_colorbar(ticks = FALSE))
  plot(g2)
}


.seriesRegressionTuning=function(object,xValidation=NULL,kernelNewOld=NULL,zValidation,userCall=TRUE)
{
  # returns the best choice of I and J, the number of components to use in each
  # direction
  # matrixErrors: should the matrix with errors be returned? default is false
  # zValidation expected to be between 0 and 1
  if(userCall)
    cat("\n A user will typically use function 'seriesRegression' instead, are you sure about what you are doing?")
  
  
  if(class(object)!="seriesRegression") stop("Object should be of class seriesRegression")
  
  if(is.null(kernelNewOld))
  {
    kernelNewOld=object$kernelFunction(rdist(xValidation,object$xSelectedTrain),object$extraKernel)
  }
  
  
  if(any(is.na(kernelNewOld))) stop("Kernel with NA")
  
  nX=object$nX
  
  m=dim(kernelNewOld)[1] # New
  n=dim(kernelNewOld)[2] # Old
  
  if(object$basis=="KPCA")
  {
    basisX=kernelNewOld %*% object$eigenX
    basisX=1/n*basisX*matrix(rep(1/object$eigenValuesX,m),m,nX,byrow=T)
  } else { # Diffusion basis
    sumsKernel=rowSums(kernelNewOld) #  for each new feature, sum of kernel wrt old features
    finalResult=kernelNewOld%*%object$eigenX
    finalResult=finalResult/sumsKernel
    basisX=finalResult*matrix(rep(1/(n*object$eigenValuesX),m),m,nX,byrow=T)
  }
  
  
  errors=apply(as.matrix(1:nX),1,function(kk)
  {
    betaHat=t(object$coefficients[1:kk,drop=FALSE])
    predictedY=basisX[,1:kk,drop=F]%*%t(betaHat)
    return(mean((predictedY-zValidation)^2))
  })
  object$nXBest=(1:nX)[which.min(errors)]  
  object$bestError=min(errors)
  gc(verbose=F)
  object$errors=errors
  
  
  
  return(object)
}




predictRegression=function(object,xTest=NULL,kernelTestTrain=NULL)
{
  
  if(class(object)!="seriesRegression") stop("Object should be of class seriesRegression")
  
  if(is.null(kernelTestTrain)&is.null(object$xSelectedTrain)) 
    stop("It seems you trained the regression using the kernelTrainTrain argument, please use kernelTestTrain instead of xTest here")
  
  
  if(is.null(kernelTestTrain))
  {
    kernelNewOld=object$kernelFunction(rdist(xTest,object$xSelectedTrain),object$extraKernel)
  } else {
    kernelNewOld=kernelTestTrain
    rm(kernelTestTrain)
  }
  
  
  nXBest=object$nX
  
  if(!is.null(object$nXBest)) # if CV was done
  {
    nXBest=object$nXBest
  }  
  
  m=dim(kernelNewOld)[1] # New
  n=dim(kernelNewOld)[2] # Old
  
  if(object$basis=="KPCA")
  {
    basisX=kernelNewOld %*% object$eigenX[,1:nXBest,drop=FALSE]
    basisX=1/n*basisX*matrix(rep(1/object$eigenValuesX[1:nXBest],m),m,nXBest,byrow=T)
  } else { # Diffusion
    sumsKernel=rowSums(kernelNewOld) #  for each new feature, sum of kernel wrt old features
    finalResult=kernelNewOld%*%object$eigenX[,1:nXBest,drop=F]
    finalResult=finalResult/sumsKernel
    basisX=finalResult*matrix(rep(1/(n*object$eigenValuesX[1:nXBest]),m),m,nXBest,byrow=T)
    
  }
  nTextX=dim(kernelNewOld)[1] # how many test points for x
  
  predictedY=basisX%*%object$coefficients[1:object$nXBest]
  
  return(predictedY)  
}
