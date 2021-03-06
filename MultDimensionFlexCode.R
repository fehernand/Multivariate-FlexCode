require(slam)
require(FlexCoDE)
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
#' @param weight - A vector of weights. 
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
#'                                weight = FALSE) 
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
                                  weight = FALSE){
  
  conditionalDensities <- multMarginalConditionalDensityFlexCode__(xTrain=xTrain,
                                                                   zTrain = zTrain,
                                                                   xValidation = xValidation,
                                                                   zValidation = zValidation,
                                                                   regressionFunction = regressionFunction)
  
  if(weight == TRUE){
    
  
      resultParamsweight <- fitRandomForestweight__(xTrain = xTrain,
                                                  zTrain = zTrain,
                                                  xValidation = xValidation,
                                                  zValidation = zValidation,
                                                  conditionalDensities = conditionalDensities,
                                                  copulaFunction = copulaFunction)
    
    fitweight <- resultParamsweight$fitForest
    
    out <- structure(list(weight = weight,
                          fitweight = resultParamsweight$fitForest,
                          hweight = resultParamsweight$h,
                          conditionalFit = conditionalDensities$fit,
                          conditionalDensities = conditionalDensities,
                          copulaFunction = copulaFunction,
                          xValidation = xValidation), class = 'multFlexCode')
    
  }
  if(weight == "euclidean"){
    
    
    
    resultParamsweight <- crossValidationWeight__(conditionalDensities = conditionalDensities,
                            copulaFunction =  copulaFunction,
                            xValidation = xValidation,
                            type = "euclidean",
                            kfold = 3,
                            h_limits = c(0.05, 2),
                            h_by = 0.25)
    
    out <- structure(list(weight = weight,
                          hweight = resultParamsweight$h_best,
                          conditionalFit = conditionalDensities$fit,
                          conditionalDensities = conditionalDensities,
                          copulaFunction = copulaFunction,
                          xValidation = xValidation), class = 'multFlexCode')
    
  }
  
  if(weight == FALSE){ 
    
    copulaAjust <- copulaParamAjust__(conditionalDensities,
                                      copulaFunction)
    
    out <- structure(list(op = copulaAjust$op,
                          weight = weight,
                          conditionalFit = conditionalDensities$fit,
                          copulaFunction = copulaFunction), class = 'multFlexCode')
  }
  
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
  fCumulative[fCumulative == 0] <- 0.0000001
  
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
                       nIMax = 40,regressionFunction = regressionFunction, chooseSharpen = T) 
    
    densityValidation <- densityFlexCode__(pred = predict(fit ,xValidation),
                                           z = zValidation[,i])
    
    densityTrain <- densityFlexCode__(pred = predict(fit ,xTrain),
                                           z = zTrain[,i])
    
    result[["validation"]][["s"]][[i]] <- densityValidation$S
    result[["validation"]][["f"]][[i]] <- densityValidation$f
    
    result[["train"]][["s"]][[i]] <- densityTrain$S
    result[["train"]][["f"]][[i]] <- densityTrain$f
    
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
                               weights = NULL){
  
  llCopulaFunction <- selectllCopulaFunction__(copulaFunction)
  
  if(copulaFunction == 'gumbel'){
    lowerLimit = 0
    upperLimit = 1
  }
  if(copulaFunction == 'clayton'){
    lowerLimit = -1
    upperLimit = 100
  }
  if(copulaFunction == 'ig'){
    lowerLimit = 0
    upperLimit = 100
  }
  
  
  if(is.null(weights) == FALSE){
    
    op <- matrix(rep(NA,nrow(weights)*3),nrow=nrow(weights))
    
    for(i in 1:nrow(weights)){
      op_temp <- optim(par=0.5,fn=llCopulaFunction,w=weights[i,],Data=conditionalDensities$validation
                       ,method="Brent",control=list(fnscale=-1),lower=lowerLimit,upper=upperLimit)
      op[i,] <- c(op_temp$value,op_temp$par,op_temp$convergence)
    }
  }
  
  if(is.null(weights) == TRUE){
    op <- optim(par = 0.5, fn = llCopulaFunction, Data=conditionalDensities$validation, method="Brent",
                control=list(fnscale=-1),lower=lowerLimit,upper=upperLimit)
    
  }
  return(list(op = op))
}

#' (Internal) Random Forest weight
#'
#' @description
#' This function performes random forest technique to calculate weights.
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
fitRandomForestweight__ <- function(xTrain,
                                    zTrain,
                                    xValidation,
                                    zValidation,
                                    conditionalDensities,
                                    copulaFunction
                                    ){
  
  colnames(xTrain) = NULL
  
  fitForest <- list()
  
  for(i in 1:ncol(zTrain)){
    fitForest[[i]] <- randomForest::randomForest(x=xTrain,y=zTrain[,i])
  }
  
  
  result_cv <- crossValidationWeight__(conditionalDensities = conditionalDensities,
                                             copulaFunction = copulaFunction,
                                             fitForest = fitForest,
                                             xValidation = xValidation,
                                             kfold = 2,
                                       h_limits = c(0.05, 7), h_by = 1)
  
  return(list("fitForest" = fitForest, "h" = result_cv$h_best))
}

crossValidationWeight__ <- function(conditionalDensities, copulaFunction , xValidation, type = "randomforest", fitForest = NULL, kfold = 3, h_limits = c(0.05, 5), h_by = 0.5 ){
  
  indiceKfold <- rep(1:kfold, length.out = nrow(xValidation))[sample(1:nrow(xValidation))]
  
  h_grid = seq(from = h_limits[1], to = h_limits[2], by = h_by)
  risk_cross <- c()
  for(h in 1:length(h_grid)){
    
    risk_kfold <- c()
    for(i in 1:kfold){
    
      if(type == "randomforest"){
        weights <- randomForestweight__(fitweight = fitForest,
                                        newX = xValidation[indiceKfold == i,],
                                        xValidation = xValidation[indiceKfold != i,],
                                        h = h_grid[h])
      }
      if(type == "euclidean"){
        weights <- euclideanDistance(xValidation = xValidation[indiceKfold != i,],
                                     xTest = xValidation[indiceKfold == i,],
                                     columns = 1:ncol(xValidation),
                                     h = h_grid[h])
      }
      
      
      # quebrar conditional densities em teste e treino
      
      dataTemp_treino <- list()
      dataTemp_teste <- list()
      
      dataTemp_treino$validation$s <- lapply(conditionalDensities$validation$s, function(x) x[indiceKfold != i])
      dataTemp_treino$validation$f <- lapply(conditionalDensities$validation$f, function(x) x[indiceKfold != i])
      
      dataTemp_teste$validation$s <- lapply(conditionalDensities$validation$s, function(x) x[indiceKfold == i])
      dataTemp_teste$validation$f <- lapply(conditionalDensities$validation$f, function(x) x[indiceKfold == i])
        
      par <- copulaParamAjust__(conditionalDensities = dataTemp_treino,
                                      copulaFunction,
                                      weights = weights)$op[,2]
      
      llCopula <- selectllCopulaFunction__(copulaFunction)
      
      values_temp <- c()
      
      for(w in 1:length(par)){
        individual_data <- data.frame("s" = c(NA,NA), "f" = c(NA,NA))
        individual_data[,1] <- as.data.frame(sapply(dataTemp_teste$validation$s, function(x) x[w]))
        individual_data[,2] <- as.data.frame(sapply(dataTemp_teste$validation$f, function(x) x[w]))
        
        
        values_temp[w] <- llCopula(par[w], w = 1, Data = individual_data, op = 0) 
      }
      
      risk_kfold[i] <- mean(values_temp)
    }
    risk_cross[h] <- mean(risk_kfold)
  }
  
  return(list("risk" = risk_cross,"h" = h_grid, "h_best" = h_grid[risk_cross == max(risk_cross)] ))
}

randomForestweight__ <- function(fitweight , newX, xValidation, h){
  
  
  wCoordenate <- list()
  
  for(i in 1:length(fitweight)){
    wCoordenate[[i]] <- wForest__(xValidation = xValidation, newX = newX, fit_forest_train = fitweight[[i]])
    wCoordenate[[i]] <- exp(-(1/(wCoordenate[[i]]*h)))
    wCoordenate[[i]] <- wCoordenate[[i]]/rowSums(wCoordenate[[i]])
  }
  
  
  # wCoordenate <- wForest__(xValidation = xValidation, newX = newX, fit_forest_train = fitweight)
  # wCoordenate <- exp(-(1/(wCoordenate*h)))
  # wCoordenate <- wCoordenate/rowSums(wCoordenate)
  
  #weights <- wCoordenate
  weights <- Reduce('+', wCoordenate)/length(wCoordenate)
  
  return(weights)
}


#### Euclidean Distance

euclideanDistance <- function(xValidation, xTest, h, columns = 1:5){
  
  weights <- matrix(NA,ncol = nrow(xValidation), nrow = nrow(xTest))
  
  names(xTest) = names(xValidation) = as.character(1:ncol(xValidation))
  
  tempScale <- scale(rbind(xValidation, xTest))
  
  xValidation <- tempScale[1:nrow(xValidation),]
  
  xTest <- tempScale[-c(1:nrow(xValidation)),]
  
  
  for(i in 1:nrow(xTest)){
    
    distanciaUnitaria <- as.matrix(xValidation[,columns]) - t(matrix(rep(as.numeric(xTest[i,columns]),nrow(xValidation)), ncol = nrow(xValidation)))
    
    distanciaUnitaria <- (distanciaUnitaria)^2
    
    distanciaUnitaria[is.na(distanciaUnitaria)] = 0
    
    distanceTemp <- sqrt(rowSums(distanciaUnitaria) + 0.0001)
    
    weights[i,] <- exp(-(distanceTemp/(h)))  
    
  }
  
  return(weights)
  
}



#' (Internal) Random Forest auxiliar proximate weight
#'
#' @description
#' This function performes random forest technique to calculate weights.
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
wForest__ <- function(xValidation, newX, fit_forest_train){
  newX <- as.matrix(newX)
  xValidation <- as.matrix(xValidation)
  
  xPred <- rbind(newX,xValidation)
  colnames(xPred) = NULL
  
  aux=predict(fit_forest_train,xPred,proximity = TRUE)
  proximidade=aux$proximity[1:nrow(newX),-c(1:nrow(newX))] 
  return(proximidade)
}





#' (Internal) Predict Mult Dimensional FlexCode
#'
#' @description
#' This function performes random forest technique to calculate weights.
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

predict.multFlexCode = function(model, newX, h = NULL){
  
  cPred <- conditionalPred(fit = model$conditionalFit, newX)
  
  #### Calcula a densidade no ponto e 1 - a acumulada para todo grid pred$z
  Data <- predMultMarginalConditionalDensityFlexCode__(cPred)
  
  # Select like log copula function
  llCopula <- selectllCopulaFunction__(model$copulaFunction)
  
  # If weighted
  if(model$weight == TRUE){
    
    if(is.null(h)){
      h <- model$hweight
    }
    
    
    F1 = 1 - model$conditionalDensities$validation$s[[1]]
    
    weights <- randomForestweight__(fitweight = model$fitweight,
                                    newX = newX,
                                    xValidation = model$xValidation,
                                    h = h)
    
    par <- copulaParamAjust__(conditionalDensities = model$conditionalDensities,
                              copulaFunction = model$copulaFunction,
                              weights = weights)$op[,2]
  }
  if(model$weight == "euclidean"){
    
    if(is.null(h)){
      h <- model$hweight
    }
    
    weights <- euclideanDistance(xValidation = model$xValidation,
                                 xTest = newX,h = h, columns = 1:ncol(model$xValidation))
    
    par <- copulaParamAjust__(conditionalDensities = model$conditionalDensities,
                              copulaFunction = model$copulaFunction,
                              weights = weights)$op[,2]
    
  }
  
  if(model$weight == FALSE){
    
    par <- rep(model$op$'par', nrow(newX))
  }
  
  
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
      
      values[j] <- exp(llCopula(par[k], w = 1, Data = dataTemp, op = 0))
      
    }
    
    CDE[[k]] <- simple_sparse_array(as.matrix(expGrid), values, dim = rep(1000, d))
  }
  
  z <- list()
  
  for(i in 1:d){
    z[[i]] <- cPred[[i]]$z
  }
  
  return(list('CDE' = CDE, 'z' = z, 'newX' = newX, 'copulaFunction' = model$copulaFunction, 'conditionalFit' = model$conditionalFit, 'par' = par))
}

#' (Internal) Conditional Predict Mult Dimensional FlexCode
#'
#' @description
#' This function performes random forest technique to calculate weights.
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
    predConditionalDensity <- conditionalDensityResult__(pred = cPred[[i]])
    
    result[["s"]][[i]] <- predConditionalDensity$S
    result[["f"]][[i]] <- predConditionalDensity$f
  }
  
  return(result)
}

conditionalDensityResult__ <- function(pred){
  # Z1 Valida??o density 
  
  pred$CDE <- pred$CDE/mean(rowSums(pred$CDE)*(pred$z[2]-pred$z[1]))
  
  dimValuesCDE <- which(pred$CDE !=0, arr.ind = T)
  
  dimValuesCDE <- dimValuesCDE[order(dimValuesCDE[,1]),]
  
  distF <- matrix(0,nrow = nrow(pred$CDE), ncol = ncol(pred$CDE))
  for(i in 1:nrow(pred$CDE)){
    
    dimTemp <- dimValuesCDE[dimValuesCDE[,1] == i,2]
    
    for(j in 1:length(dimTemp)){
      
      distF[i,dimTemp[j]] <- sum((pred$z[2]-pred$z[1])*pred$CDE[i,1:dimTemp[j]])
    }
  }
  
  distF[distF>=1] <- 1
  
  distS <- 1-distF
  
  distS[distS == 1] = 0
  
  return(list("f" = pred$CDE, "S" = distS))
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


