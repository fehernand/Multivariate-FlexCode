

gridDados <- function(caso, xTest, z){
  
  n = nrow(xTest)
  
  x1 <- xTest[,1]
  x2 <- xTest[,2]
  x3 <- xTest[,3]
  x4 <- xTest[,4]
  x5 <- xTest[,5]
  
  mu1   <- 3 + 0.5*x1 + 5*x2 - 1*x3 + 3*x4 - 1.1*x5
  mu2   <- 5 + 0.1*x1 + 2*x2 - 3*x3 + 3*x4 + 1.2*x5
  mu    <- cbind(mu1,mu2)
  Sigma <- matrix(c(1,0.5,0.5,1),ncol=2,byrow=TRUE)
  
  
  grid <- expand.grid(z[[1]], z[[2]])
  CDE <- list()
  if(caso == 1){
    for(i in 1:n){ CDE[[i]] <- matrix(mvtnorm::dmvnorm(x = grid,mean = as.numeric(mu[i,]),sigma = Sigma), ncol = 1000, nrow = 1000) }
  }
  
  if(caso == 2){
    
    for(i in 1:n){
      Sigma <- matrix(c(x1[i]+0.5, 0,  0, x1[i]+0.5), ncol = 2, nrow = 2, byrow = T)
      
      CDE[[i]] <- matrix(mvtnorm::dmvnorm(x = grid,mean = as.numeric(mu[i,]),sigma = Sigma), ncol = 1000, nrow = 1000)
    }
  }
  
  if(caso == 3){
    
    for(i in 1:n){
      Sigma <- matrix(c(x1[i]+0.5, x1[i]/2,  x1[i]/2, x1[i]+0.5), ncol = 2, nrow = 2, byrow = T)
      
      CDE[[i]] <- matrix(mvtnorm::dmvnorm(x = grid,mean = as.numeric(mu[i,]),sigma = Sigma), ncol = 1000, nrow = 1000)
    }
  }
  return(CDE)
}



gerarDados <- function(caso = 1, n.Test = 250, n.TV = 5000){
  
  n <- n.Test+n.TV
  
  x1 <- runif(n, 0,1)
  x2 <- rgamma(n,shape=3,rate=4)
  x3 <- rnorm(n,3,4)
  x4 <- rpois(n,10)
  x5 <- rbeta(n,1/2,1/4)
  
  # Cov. Irrelevantes
  
  aux.ci <- matrix(rnorm(n*95),ncol=95) # 95 cov. irrelevante
  
  
  mu1   <- 3 + 0.5*x1 + 5*x2 - 1*x3 + 3*x4 - 1.1*x5
  mu2   <- 5 + 0.1*x1 + 2*x2 - 3*x3 + 3*x4 + 1.2*x5
  mu    <- cbind(mu1,mu2)
  Sigma <- matrix(c(1,0.5,0.5,1),ncol=2,byrow=TRUE)
  
  
  if(caso == 1){
    data <- matrix(rep(NA,n*2),ncol=2)
    for(i in 1:n){ data[i,] <- mvtnorm::rmvnorm(1,as.numeric(mu[i,]),Sigma) }
  }
  
  if(caso == 2){
    
    data <- matrix(rep(NA,n*2),ncol=2)
    for(i in 1:n){
      Sigma <- matrix(c(x1[i]+0.5, 0,  0, x1[i]+0.5), ncol = 2, nrow = 2, byrow = T)
      
      data[i,] <- mvtnorm::rmvnorm(1,as.numeric(mu[i,]),Sigma)
    }
  }
  
  if(caso == 3){
    
    data <- matrix(rep(NA,n*2),ncol=2)
    for(i in 1:n){
      Sigma <- matrix(c(x1[i]+0.5, x1[i]/2,  x1[i]/2, x1[i]+0.5), ncol = 2, nrow = 2, byrow = T)
      
      data[i,] <- mvtnorm::rmvnorm(1,as.numeric(mu[i,]),Sigma)
    }
  }
  
  if(caso == 4){
    aux.1 <- rep(NA,n)
    aux.2 <- rep(NA,n)
    for(i in 1:n){
      aux.1[i] <- rlogis(1,mu1[i],5) 
      aux.2[i] <- rnorm(1,-mu2[i]*aux.1[i],2)
    }
    data <- cbind(aux.1,aux.2)
  }
  
  out <- as.data.frame(cbind(data,x1,x2,x3,x4,x5,aux.ci))
  
  names(out) = c("z1", "z2",paste("X", 1:100, sep = ""))
  
  return(out)
}


splitData <- function(data, n.Test = 25){
  
  n <- nrow(data)
  
  ro <- sample(1:n)
  
  xTest       <- data[ro[1:n.Test],-c(1:2)] # nrow(xTest)=1000 (fixo)
  xTrain      <- data[ro[(n.Test+1):(n.Test+(n-n.Test)*0.7)],-c(1:2)] 
  xValidation <- data[ro[((n.Test+(n-n.Test)*0.7)+1):n],-c(1:2)]
  
  colnames(xTest)=colnames(xTrain)=colnames(xValidation)=NULL
  
  zTest       <- data[ro[1:n.Test],c(1:2)] # nrow(xTest)=1000 (fixo)
  zTrain      <- data[ro[(n.Test+1):(n.Test+(n-n.Test)*0.7)],c(1:2)] 
  zValidation <- data[ro[((n.Test+(n-n.Test)*0.7)+1):n],c(1:2)]
  
  return(list(xTest = xTest,
              xTrain = xTrain,
              xValidation = xValidation,
              zTest = zTest,
              zTrain = zTrain,
              zValidation = zValidation))
}


library(np)

library(randomForest)

wForest <- function(xValidation, newX, fit_forest_train){
  newX <- as.matrix(newX)
  xValidation <- as.matrix(xValidation)
  aux=predict(fit_forest_train,rbind(newX,xValidation),proximity = TRUE)
  proximidade=aux$proximity[1:nrow(newX),-c(1:nrow(newX))] 
  return(proximidade)
}


## Função Kernel
riscoKernelFunction <- function(xTrain, zTrain, xValidation, zValidation, xTest, zTest){
  
  fit1 <- fitFlexCoDE(xTrain,zTrain[,1],xValidation,zValidation[,1],
                      nIMax = 30,regressionFunction = regressionFunction.XGBoost) #XGBoost
  
  # CD p/ Z2
  fit2 <- fitFlexCoDE(xTrain,zTrain[,2],xValidation,zValidation[,2],
                      nIMax = 30,regressionFunction = regressionFunction.XGBoost) #XGBoost
  
  
  pred1  <- predict(fit1,xValidation) 
  pred2 <- predict(fit2,xValidation) 
  
  # lat
  f1 <- rep(NA,length(xValidation[,1]))
  F1 <- rep(NA,length(xValidation[,1]))
  for(i in 1:length(xValidation[,1])){
    index=FNN::knnx.index(data=pred1$z,query = zValidation[i,1],k=1)
    f1[i] <- pred1$CDE[i,index]
    F1[i] <- sum((pred1$z[2]-pred1$z[1])*pred1$CDE[i,1:index])
  }
  
  for(i in 1:length(F1)){
    if(f1[i]==0) f1[i] <- 0.0000001
    if(F1[i]>=1) F1[i] <- 0.9999999
    if(F1[i]==0) F1[i] <- 0.0000001
  }
  S1 <- 1-F1
  
  # lon
  f2 <- rep(NA,length(xValidation[,1]))
  F2 <- rep(NA,length(xValidation[,1]))
  for(i in 1:length(xValidation[,1])){
    index=FNN::knnx.index(data=pred2$z,query = zValidation[i,2],k=1)
    f2[i] <- pred2$CDE[i,index]
    F2[i] <- sum((pred2$z[2]-pred2$z[1])*pred2$CDE[i,1:index])
  }
  
  for(i in 1:length(F2)){
    if(f2[i]==0) f2[i] <- 0.0000001
    if(F2[i]>=1) F2[i] <- 0.9999999
    if(F2[i]==0) F2[i] <- 0.0000001
  }
  S2 <- 1-F2
  
  pred11 <- predict(fit1,xTest)
  pred22 <- predict(fit2,xTest)
  
  # lat
  f11 <- rep(NA,length(xTest[,1]))
  F11 <- rep(NA,length(xTest[,1]))
  F11_CDE <- matrix(NA,ncol = 1000, nrow = length(xTest[,1]))
  for(i in 1:length(xTest[,1])){
    index=FNN::knnx.index(data=pred11$z,query = zTest[i,1],k=1)
    f11[i] <- pred11$CDE[i,index]
    F11[i] <- sum((pred11$z[2]-pred11$z[1])*pred11$CDE[i,1:index])
    for(j in 1:1000){
      F11_CDE[i,j] <- sum((pred11$z[2]-pred11$z[1])*pred11$CDE[i,1:j])
    }
  }
  
  for(i in 1:length(F11)){
    if(f11[i]==0) f11[i] <- 0.0000001
    if(F11[i]>=1) F11[i] <- 0.9999999
    if(F11[i]==0) F11[i] <- 0.0000001
  }
  
  F11_CDE[F11_CDE == 0] <- 0.0000001
  F11_CDE[F11_CDE >= 1] <- 0.9999999
  
  
  # lon
  f22 <- rep(NA,length(xTest[,1]))
  F22 <- rep(NA,length(xTest[,1]))
  F22_CDE <- matrix(NA,ncol = 1000, nrow = length(xTest[,1]))
  for(i in 1:length(xTest[,1])){
    index=FNN::knnx.index(data=pred22$z,query = zTest[i,2],k=1)
    f22[i] <- pred22$CDE[i,index]
    F22[i] <- sum((pred22$z[2]-pred22$z[1])*pred22$CDE[i,1:index])
    for(j in 1:1000){
      F22_CDE[i,j] <- sum((pred22$z[2]-pred22$z[1])*pred22$CDE[i,1:j])
    }
  }
  
  for(i in 1:length(F22)){
    if(f22[i]==0) f22[i] <- 0.0000001
    if(F22[i]>=1) F22[i] <- 0.9999999
    if(F22[i]==0) F22[i] <- 0.0000001
  }
  
  F22_CDE[F22_CDE == 0] <- 0.0000001
  F22_CDE[F22_CDE >= 1] <- 0.9999999
  
  colnames(xTrain) =NULL
  fit.forest.train1=randomForest::randomForest(x=xTrain,y=zTrain[,1])
  w.coordenada.1=wForest(xValidation,as.matrix(xTest),fit.forest.train1)
  w.coordenada.1=w.coordenada.1/rowSums(w.coordenada.1)
  
  fit.forest.train2=randomForest::randomForest(x=xTrain,y=zTrain[,2])
  w.coordenada.2=wForest(xValidation,as.matrix(xTest),fit.forest.train2)
  w.coordenada.2=w.coordenada.2/rowSums(w.coordenada.2)
  
  wForest_value =(w.coordenada.1+w.coordenada.2)/2
  
  u_validation <- cbind(F1,F2)
  bandwidth <- ks::Hpi(u_validation)
  
  a1 <- rep(NA,length(zTest[,1]))
  CDE <- list()
  for(i in 1:length(zTest[,1])){
    grid.points <- expand.grid(F11_CDE[i,], F22_CDE[i,])
    
    pred_temp <- ks::kde(u_validation,
                         h = bandwidth, 
                         eval.points = grid.points,
                         compute.cont = FALSE,
                         unit.interval = FALSE, # fazer no transformado!
                         w = wForest_value[i,])
    
    a1[i] <- ks::kde(u_validation,
                     h = bandwidth, 
                     compute.cont = FALSE,
                     unit.interval = FALSE,
                     eval.points = cbind(F11[i],F22[i]), # fazer no transformado!
                     w = wForest_value[i,])$estimate
    
    CDE[[i]] <- matrix(pred_temp$estimate, ncol = 1000, nrow = 1000, byrow = F)
    
  }
  
  RKernel    <- sum(log(f11*f22*a1))/length(a1)
  SD_RKernel <- sd(log(f11*f22*a1))
  
  return(list("risk" = RKernel, "sd" = SD_RKernel,"CDE" = CDE, "grid" = grid.points))
}


