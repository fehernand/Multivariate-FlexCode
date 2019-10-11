require(np)
load("resultG.rda")
load("dados2.rda")

# Estimando modelo (band)
bw1 <- npcdensbw(xdat=data.split$xTrain, ydat=data.split$zTrain, ftol = 0.01, tol=0.01,nmult=1,bwmethod="normal-reference")

# 
z1_grid <- predG[[2]]$z[[1]]
z2_grid <- predG[[2]]$z[[2]]

size <- 300

z1_grid_temp <- seq(min(z1_grid),max(z1_grid), length.out = 300)
z2_grid_temp <- seq(min(z2_grid),max(z2_grid), length.out = 300)
risk <- c()
CDE <- matrix(NA, nrow = nrow(data.split$zTest), ncol = size^2)

for(i in 1:nrow(data.split$zTest)){
  z1_test <- data.split$zTest[i,1]
  z2_test <- data.split$zTest[i,2]
  
  # #size valor mais proximos do z observado
  #z1_index <- sort(FNN::knnx.index(data=z1_grid,query =z1_test,k=size))
  #z2_index <- sort(FNN::knnx.index(data=z2_grid,query =z2_test,k=size))
  
  # grid de valores de z1 e z2 calculados
  #z1_grid_temp[i,] <- z1_grid[z1_index]
  #z2_grid_temp[i,] <- z2_grid[z2_index]
  
  # mudando formato para estimacao da densidade condicional
  #Z.eval <- expand.grid(z1 = z1_grid[z1_index], z2 = z2_grid[z2_index])
  
  Z.eval <- expand.grid(z1 = z1_grid_temp, z2 = z2_grid_temp)
  
  gridX <- data.split$xTest[rep(i, size^2),]
  
  # Calculando CDE
  CDE[i,] <- fitted(npcdens(bw1, exdat = gridX, eydat = Z.eval  ,ftol = 0.01, tol=0.01))
  
  # Calculando densidade estimada nos pontos observados zTest
  risk[i] <- fitted(npcdens(bw1, exdat = data.split$xTest[i,], eydat = data.split$zTest[i,] ,ftol = 0.01, tol=0.01))
  
}

predKernel2 <- list(CDE = CDE, z.eval = risk, z1_grid = z1_grid_temp, z2_grid_temp = z2_grid_temp)

save(predKernel2, file = "resultKernel2.rda")

