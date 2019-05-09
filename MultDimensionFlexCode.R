require(slam)

#' Fit MultDimension FlexCode
#'
#' @description
#' This function fit MultDimension FlexCode
#'
#' @section Warning:
#' For  some unknown reason the R check won't run the examples, but they work as
#' expected.
#'
#' @section Maintainers:
#' FelipeEstat
#'
#' @author FelipeEstat
#'
#' @param xTrain - Object of class dataframe. Covariate Database train.
#' @param zTrain - Object of class dataframe. Response Database train.
#' @param xValidation - Object of class dataframe. Covariate Database Validation.
#' @param zValidation - Object of class dataframe. Response Database Validation.
#' @param nIMax - Maximum possible number of components of the series expansion (that is, the function will find the best I<nIMax).
#'  Default is 30.
#' @param regressionFunction - a function indicating which regression method will be used to estimate the expansion coefficients.
#'  Currently can be one of regressionFunction.NN, regressionFunction.NW, regressionFunction.SpAM, regressionFunction.Series,
#'  regressionFunction.Lasso, regressionFunction.Forest or regressionFunction.XGBoost. Type ?regressionFunction.XX to find out
#'  more about method XX.
#' @param copulaFunction - Copula function to be used. Can be 'gumbel', 'clayton', 'ig'.
#' @param weigth - A vector of weigths. 
#'   
#' 
#' @examples
#' \dontrun{
#' 
#'  library(FlexCoDE)
#'  reg <- regressionFunction.XGBoost
#' 
#'  model1 <- MultDimensionFlexCoDE(xTrain = xTrain,
#'                                zTrain = zTrain,
#'                                xValidation = xValidation,
#'                                zValidation = zValidation,      
#'                                nIMax = 30,         
#'                                regressionFunction = reg,
#'                                copulaFunction = 'gumbel', 
#'                                weigth = FALSE) 
#'                                   
#'                                                                                                                                
#' }
#'
#' @return Object of class MultFlexCode
#'
#'
MultDimensionFlexCoDE <- function(xTrain,                     
                                zTrain,                     
                                xValidation,                
                                zValidation,  
                                regressionFunction,
                                copulaFunction = 'gumbel',  
                                weigth = FALSE){
  
  conditionalDensities <- multMarginalConditionalDensityFlexCode__(xTrain=xTrain,
                                                                   zTrain = zTrain,
                                                                   xValidation = xValidation,
                                                                   zValidation = zValidation,
                                                                   regressionFunction = regressionFunction)
  
  resultWeigths <- if(weigth == TRUE) randomForestWeigth__(xTrain, zTrain, xValidation)
             else NULL
  
  copulaAjust <- copulaParamAjust__(conditionalDensities,
                                    copulaFunction,
                                    weigths = resultWeigths$weigths)
  
  out <- structure(list(op = copulaAjust$op,
                        weigth = weigth,
                        fitWeigth = resultWeigths$fitWeigths,
                        conditionalFit = conditionalDensities$fit,
                        copulaFunction = copulaFunction), class = 'multFlexCode')
  
  return(out)
}  

#' (Internal) Density Function Estimated by FlexCode
#'
#' @description
#' This function calculate probability density and (1 - culmulative density). Values to be used in copula function.
#'
#' @section Warning:
#' For  some unknown reason the R check won't run the examples, but they work as
#' expected.
#'
#' @section Maintainers:
#' FelipeEstat
#'
#' @author FelipeEstat
#'
#' @param fit - List of fitted objects
#' 
#' @examples
#' \dontrun{
#' 
#'                                   
#'                                                                                                                                
#' }
#'
#' @return list with terms f - probability density. S -  1- cumulative density.
#' 
#' @keywords internal
#'
#'
densityFlexCode__ <- function(pred, z){
  
  fDensity <- rep(NA,length(z)) 
  fCumulative <- rep(NA,length(z)) 
  
  for(i in 1:length(z)){
    index=FNN::knnx.index(data=pred$z,query = z[i],k=1)
    
    fDensity[i] <- pred$CDE[i,index]
    
    fCumulative[i] <- sum((pred$z[2]-pred$z[1])*pred$CDE[i,1:index])
  }
  
  fDensity[fDensity == 0] <- 0.0000001
  fCumulative[fCumulative >= 1] <- 0.9999999
  fCumulative[fCumulative == 1] <- 0.0000001
  
  return(list("f" = fDensity, "S" = 1 - fCumulative))
}


#' (Internal) Multiple Marginal Conditional Density Estimated by FlexCode
#'
#' @description
#' This function performes multiples conditional density 
#'
#' @section Warning:
#' For  some unknown reason the R check won't run the examples, but they work as
#' expected.
#'
#' @section Maintainers:
#' FelipeEstat
#'
#' @author FelipeEstat
#'
#' @param fit - List of fitted objects
#' 
#' @examples
#' \dontrun{
#' 
#'                   
#'                                                                                        
#' }
#'
#' @return list
#' 
#' @keywords internal
#'
#'
multMarginalConditionalDensityFlexCode__ <- function(xTrain,
                                                    zTrain,
                                                    xValidation,
                                                    zValidation,
                                                    regressionFunction){
  
  result <- list()
  for(i in 1:ncol(zTrain)){
    
    fit <- fitFlexCoDE(xTrain, zTrain[,i], xValidation, zValidation[,i],
                        nIMax = 30,regressionFunction = regressionFunction) 
   
    densityValidation <- densityFlexCode__(pred = predict(fit ,xValidation),
                               z = zValidation[,i])
    
    result[["validation"]][["s"]][[i]] <- densityValidation$S
    result[["validation"]][["f"]][[i]] <- densityValidation$f
  
    result[["fit"]][[i]] <- fit
  }
  
  return(result)
}


#' (Internal) Select Copula Function
#'
#' @description
#' This function select copula function.
#'
#' @section Warning:
#' For  some unknown reason the R check won't run the examples, but they work as
#' expected.
#'
#' @section Maintainers:
#' FelipeEstat
#'
#' @author FelipeEstat
#'
#' @param copulaFunction - can be 'gumbel', 'clayton, 'ig'
#' 
#' @examples
#' \dontrun{
#' 
#'                                   
#'                                                                                                                                
#' }
#'
#' @return Class function
#' 
#' @keywords internal
#'
#'
selectllCopulaFunction__ <- function(copulaFunction){
  #require(CouplaFunctions)
  source("CopulaFunctions.R")
  
  if(copulaFunction == 'gumbel'){llfn = llGumbel}
  if(copulaFunction == 'clayton'){llfn = llClayton}
  if(copulaFunction == 'ig'){llfn = llIG}
  
  return(llfn)
}

#' (Internal) Select Copula Function
#'
#' @description
#' This function select copula function.
#'
#' @section Warning:
#' For  some unknown reason the R check won't run the examples, but they work as
#' expected.
#'
#' @section Maintainers:
#' FelipeEstat
#'
#' @author FelipeEstat
#'
#' @param copulaFunction - can be 'gumbel', 'clayton, 'ig'
#' 
#' @examples
#' \dontrun{
#' 
#'                                   
#'                                                                                                                                
#' }
#'
#' @return Class function
#' 
#' @keywords internal
#'
#'
copulaParamAjust__ <- function(conditionalDensities,
                               copulaFunction,
                               weigths){
  
  llCopulaFunction <- selectllCopulaFunction__(copulaFunction)
  
  if(is.null(weigths) == FALSE){
    
    op <- matrix(rep(NA,nrow(weigths)*3),nrow=nrow(weigths))
    
    for(i in 1:nrow(weigths)){
      op_temp <- optim(par=0.5,fn=llCopulaFunction,w=weigths[i,],Data=conditionalDensities$validation
                       ,method="Brent",control=list(fnscale=-1),lower=0.0001,upper=0.9999)
      op[i,] <- c(op_temp$value,op_temp$par,op_temp$convergence)
    }
  }
  
  if(is.null(weigths) == TRUE){
    op <- optim(par = 0.5, fn = llCopulaFunction, Data=conditionalDensities$validation, method="Brent",
                control=list(fnscale=-1),lower=0.0001,upper=0.9999)
    
  }
  return(list(op = op))
}

#' (Internal) Random Forest weigth
#'
#' @description
#' This function performes random forest technique to calculate weigths.
#'
#' @section Warning:
#' For  some unknown reason the R check won't run the examples, but they work as
#' expected.
#'
#' @section Maintainers:
#' FelipeEstat
#'
#' @author FelipeEstat
#'
#' @param copulaFunction - can be 'gumbel', 'clayton, 'ig'
#' 
#' @examples
#' \dontrun{
#' 
#'                                   
#'                                                                                                                                
#' }
#'
#' @return Class function
#' 
#' @keywords internal
#'
#'
randomForestWeigth__ <- function(xTrain,
                                 zTrain,
                                 xValidation){
  
  colnames(xTrain) = NULL
  wCoordenate <- list()
  
  for(i in 1:ncol(zTrain)){
    fitForest <- randomForest::randomForest(x=xTrain,y=zTrain[,i])
    wCoordenate[[i]] <- wForest__(xValidation,xTest,fitForest)
    wCoordenate[[i]] <- wCoordenate[[i]]/rowSums(wCoordenate[[i]])
  }
  
  weights <- Reduce('+', wCoordenate)/length(wCoordenate)
  
  return(list(weigths = weights, fitWeigth = fitForest))
}

#' (Internal) Random Forest auxiliar proximate weigth
#'
#' @description
#' This function performes random forest technique to calculate weigths.
#'
#' @section Warning:
#' For  some unknown reason the R check won't run the examples, but they work as
#' expected.
#'
#' @section Maintainers:
#' FelipeEstat
#'
#' @author FelipeEstat
#'
#' @param copulaFunction - can be 'gumbel', 'clayton, 'ig'
#' 
#' @examples
#' \dontrun{
#' 
#'                                   
#'                                                                                                                                
#' }
#'
#' @return Class function
#' 
#' @keywords internal
#'
#'
#'
wForest__ = function(xValidation,xNew,fit_forest_train){
  aux=predict(fit_forest_train,rbind(xNew,xValidation),proximity = TRUE)
  
  proximidade=aux$proximity[1:nrow(xNew),-c(1:nrow(xNew))] 
  
  return(proximidade)
}

#' (Internal) Predict Mult Dimensional FlexCode
#'
#' @description
#' This function performes random forest technique to calculate weigths.
#'
#' @section Warning:
#' For  some unknown reason the R check won't run the examples, but they work as
#' expected.
#'
#' @section Maintainers:
#' FelipeEstat
#'
#' @author FelipeEstat
#'
#' @param copulaFunction - can be 'gumbel', 'clayton, 'ig'
#' 
#' @examples
#' \dontrun{
#' 
#'                                   
#'                                                                                                                                
#' }
#'
#' @return Class function
#' 
#' @keywords internal
#'
#'

predict.multFlexCode = function(model, newX, weigths = NULL){
  
  cPred <- conditionalPred(fit = model$conditionalFit, newX)
  
  #### Calcula a densidade no ponto e 1 - a acumulada para todo grid pred$z
  Data <- predMultMarginalConditionalDensityFlexCode__(cPred)
  
  # Select like log copula function
  llCopula <- selectllCopulaFunction__(model$copulaFunction)
  
  
  CDE <- list()
  for(k in 1:nrow(newX)){ # para cada Z
    
    d <- length(Data$s)
    
    distS <- as.data.frame(sapply(Data$s, function(x) x[k,]))
    
    dimDist <- which(distS != 0, arr.ind = T)
    
    distf <- as.data.frame(sapply(Data$f, function(x) x[k,]))
    
    grid <- list()
    
    for(i in 1:d){ 
      
      grid[[i]] <- dimDist[dimDist[,2] == i,1]
      
      
    }
    
    expGrid <- expand.grid(grid)
    
    dataTemp <- data.frame("s" = rep(NA,d), "f" = rep(NA,d))
    values <- rep(NA, nrow(expGrid))
    for(j in 1:nrow(expGrid)){
      for(dTemp in 1:d){
          dataTemp[dTemp,1] = distS[expGrid[j,dTemp],dTemp]
          dataTemp[dTemp,2] = distf[expGrid[j,dTemp],dTemp]
        }
        
        values[j] <- exp(llCopula(model$op$'par', w =1, Data = dataTemp, op = 0))
     }
    
    CDE[[k]] <- simple_sparse_array(as.matrix(expGrid), values, dim = rep(1000, d))
  }
  
  z <- list()
  
  for(i in 1:d){
    z[[i]] <- cPred[[i]]$z
  }
  
  return(list('CDE' = CDE, 'z' = z, 'newX' = newX, 'copulaFunction' = model$copulaFunction, 'conditionalFit' = model$conditionalFit, 'par' = model$op$par))
}

#' (Internal) Conditional Predict Mult Dimensional FlexCode
#'
#' @description
#' This function performes random forest technique to calculate weigths.
#'
#' @section Warning:
#' For  some unknown reason the R check won't run the examples, but they work as
#' expected.
#'
#' @section Maintainers:
#' FelipeEstat
#'
#' @author FelipeEstat
#'
#' @param copulaFunction - can be 'gumbel', 'clayton, 'ig'
#' 
#' @examples
#' \dontrun{
#' 
#'                                   
#'                                                                                                                                
#' }
#'
#' @return Class function
#' 
#' @keywords internal
#'
#'
#'
conditionalPred <- function(fit = model$fit, newX){
  pred <- list()
  for(i in 1:length(fit)){
    pred[[i]]  <- predict(fit[[i]],newX)
  }
  return(pred)
}

#' (Internal) Predict Multiple Marginal Conditional Density Estimated by FlexCode
#'
#' @description
#' This function performes multiples conditional density 
#'
#' @section Warning:
#' For  some unknown reason the R check won't run the examples, but they work as
#' expected.
#'
#' @section Maintainers:
#' FelipeEstat
#'
#' @author FelipeEstat
#'
#' @param fit - List of fitted objects
#' 
#' @examples
#' \dontrun{
#' 
#'                   
#'                                                                                        
#' }
#'
#' @return list
#' 
#' @keywords internal
#'
#'
predMultMarginalConditionalDensityFlexCode__ <- function(cPred){
  
  result <- list()
  for(i in 1:length((cPred))){
    predConditionalDensity <- conditionalDensityResult__(cPred[[i]])
    
    result[["s"]][[i]] <- predConditionalDensity$S
    result[["f"]][[i]] <- predConditionalDensity$f
  }
  
  return(result)
}

conditionalDensityResult__ <- function(pred){
  # Z1 Validação density 
  
  pred$CDE <- pred$CDE/mean(rowSums(pred$CDE)*(pred$z[2]-pred$z[1]))
  
  dimValuesCDE <- which(pred$CDE !=0, arr.ind = T)
  
  dimValuesCDE <- dimValuesCDE[order(dimValuesCDE[,1]),]
  
  distF <- matrix(0,nrow = nrow(pred$CDE), ncol = ncol(pred$CDE))
  for(i in 1:nrow(pred$CDE)){
    
    dimTemp <- dimValuesCDE[dimValuesCDE[,1] == i,2]
    
    for(j in 1:length(dimTemp)){
      
      distF[i,j] <- sum((pred$z[2]-pred$z[1])*pred$CDE[i,1:dimTemp[j]])
    }
  }
  
  distF <- pred$CDE
  
  distF[distF>=1] <- 1
  
  distS <- 1-distF
  
  distS[distS == 1] = 0
  
  return(list("f" = distF, "S" = distS))
}






calculateRisk <- function(pred, zTest){
  
  result <- list()
  for(i in 1:ncol(zTest)){
    densityValidation <- densityFlexCode__(pred = predict(pred$conditionalFit[[i]] , pred$newX),
                                           z = zTest[,i])
    
    result[["s"]][[i]] <- densityValidation$S
    result[["f"]][[i]] <- densityValidation$f
  }
  
  llCopulaFunction <- selectllCopulaFunction__(pred$copulaFunction)
  
  risco <- llCopulaFunction(pred$par,result)/length(result$f[[1]])
  sd <- sd(llCopulaFunction(pred$par,result,op=2))/sqrt(length(result$f[[1]]))
  
  # Weight
  #risco <- llCopulaFunction(op[,2],conditionalDensities)/length(conditionalDensities$f[[1]])
  #sd <- sd(llCopulaFunction(op[,2],conditionalDensities,op=2))/sqrt(length(conditionalDensities$f[[1]]))
  
  return(list(risco = risco, sd = sd))
}
  

