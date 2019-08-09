set.seed(801)

source("MultDimensionFlexCode.R")
source("ExemplosFlexCodeMulti.R")

library(FlexCoDE)

kernel <- list()

for(i in 1:3){
  load(file = paste("dados",i,".rda", sep=""))
  
  kernel[[i]] <- riscoKernelFunction(xTrain = data.split$xTrain,
                                     zTrain = data.split$zTrain,
                                     xValidation = data.split$xValidation,
                                     zValidation = data.split$zValidation,
                                     xTest = data.split$xTest,
                                     zTest = data.split$zTest)
  
}

save(kernel, file = "riscoKernel.rda")
