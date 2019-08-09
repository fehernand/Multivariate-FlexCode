require(RFCDE)
setwd("~/MultFlexCode_1")
forest <- list()
predForest <- list()
timeF <- c()
load("resultG.rda")

require(RFCDE)
for(i in 1:4){
  
  load(file = paste("dados",i,".rda", sep=""))
  
  
  n_trees <- 500
  mtry <- 10
  node_size <- 30
  n_basis <- 15
  
  start.time <- Sys.time()
  forest[[i]] <- RFCDE(data.split$xTrain, data.split$zTrain, n_trees = n_trees, mtry = mtry,
                  node_size = node_size, n_basis = n_basis)
  end.time <- Sys.time()
  timeF[i] <- end.time - start.time
  
  
  n_grid <- 30
  z_grid <- expand.grid(z1 = predG[[i]]$z[[1]],
                        z2 = predG[[i]]$z[[2]])
  
  # Verificar
  CDE <- predict(forest[[i]], as.matrix(data.split$xTest), "CDE", z_grid, bandwidth = matrix(c(10,0,0,10), ncol = 2, nrow = 2, byrow = T))
  
  risk <- predict(forest[[i]], as.matrix(data.split$xTest), "CDE", data.split$zTest, bandwidth = "cv")

  predForest[[i]] <- list(CDE = CDE, z.eval = risk)
}

save(predForest, file = "resultForest.rda")

estimate_loss(forest[[i]], x_test = data.split$xTest, z_test = data.split$zTest)
