#' (Internal) Gumbel Copula log like Function
#'
#' @description
#' Gumbel copula log like Function
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
#' @param alpha - 
#' @param Data - 
#' @param w - weigths
#' @param op - If 1 return sum.
#' @param log - Bolean operator. If true return log.
#' 
#' @examples
#' \dontrun{
#' 
#'                                   
#'                                                                                                                                
#' }
#'
#' @return Values 
#' 
#' @keywords internal
#'
#'
llGumbel <- function(alpha,Data,w=1,op=1){	
  llike <- rep(-10^(-10),length(Data$f[[1]]))
  if ((length(alpha[alpha>0])==length(alpha)) && (length(alpha[alpha<1])==length(alpha))) { 
    
    if(op == 0){
      f <- as.data.frame(matrix(Data$f, ncol = length(Data$f)), col.names = as.character(1:length(Data$f)))
      S <- as.data.frame(matrix(Data$s, ncol = length(Data$s)), col.names = as.character(1:length(Data$s)))
    }else{
      f <- as.data.frame(Data$f, col.names = as.character(1:length(Data$f)))
      S <- as.data.frame(Data$s, col.names = as.character(1:length(Data$s)))
    }
    
    A <- -log(S)
    sumlogA <- apply(log(A), 1, sum)
    sumAlphaA <- apply(A^(1/alpha), 1, sum) 
    sumlogS <- apply(log(S), 1, sum)
    sumlogf <- apply(log(f), 1, sum)
    
    
    llike <- w*((1/alpha-1)*sumlogA - sumAlphaA^alpha + (alpha - 2)*log(sumAlphaA) + log(sumAlphaA^alpha + (1-alpha)/alpha) -
         sumlogS + sumlogf)
  }
  ifelse(op==1,return(sum(llike)),return(llike))
}

#' (Internal) Clayton Copula log like Function
#'
#' @description
#' Clayton copula log like function.
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
#' @param eta - 
#' @param Data - 
#' @param w - weigths
#' @param op - If 1 return sum.
#' @param log - Bolean operator. If true return log.
#' 
#' @examples
#' \dontrun{
#' 
#'                                   
#'                                                                                                                                
#' }
#'
#' @return Values 
#' 
#' @keywords internal
#'
#'

llClayton <- function(eta,Data,w=1,op=1){
  llike <- -10^(-10)
  if ((length(eta[eta>-1])==length(eta))){
    if(op == 0){
      f <- as.data.frame(matrix(Data$f, ncol = length(Data$f)), col.names = as.character(1:length(Data$f)))
      S <- as.data.frame(matrix(Data$s, ncol = length(Data$s)), col.names = as.character(1:length(Data$s)))
    }else{
    
      f <- as.data.frame(Data$f, col.names = as.character(1:length(Data$f)))
      S <- as.data.frame(Data$s, col.names = as.character(1:length(Data$s)))
    }
    A <- apply(S^(-eta), 1, sum) - 1
    logS <- apply(log(S), 1, sum)
    logf <- apply(log(f), 1, sum)
    
    llike <- w*((-eta-1)*logS + (-1/eta-2)*log(A) + log(1+eta) + logf)
  }
  ifelse(op==1,return(sum(llike)),return(llike))
}

#' (Internal) Ig Copula log like Function
#'
#' @description
#' Ig copula log like Function
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
#' @param eta - 
#' @param Data - 
#' @param w - weigths
#' @param op - If 1 return sum.
#' @param log - Bolean operator. If true return log.
#' 
#' @examples
#' \dontrun{
#' 
#'                                   
#'                                                                                                                                
#' }
#'
#' @return Values 
#' 
#' @keywords internal
#'
#'

llIG <- function(eta,Data,w=1,op=1){
  llike <- -10^(-10)
  if((length(eta[eta>0])==length(eta))){
    if(op == 0){
      f <- as.data.frame(matrix(Data$f, ncol = length(Data$f)), col.names = as.character(1:length(Data$f)))
      S <- as.data.frame(matrix(Data$s, ncol = length(Data$s)), col.names = as.character(1:length(Data$s)))
      A <- -log(S)
    }else{
    f <- as.data.frame(Data$f, col.names = as.character(1:length(Data$f)))
    S <- as.data.frame(Data$s, col.names = as.character(1:length(Data$s)))
    }
    alpha <- 0.5
    
    logS <- apply(log(S), 1, sum)
    logf <- apply(log(f), 1, sum)
    A <- eta^alpha - alpha*eta^(alpha-1)*log(S)
    logA <- apply(log(A), 1,sum)
    Aeta <- apply(A^(1/alpha), 1, sum) - eta
    cS <- -(1/alpha)*(eta^(1-alpha)*Aeta^alpha-eta)
    
    llike <- w*((1/alpha-1)*logA + cS +(alpha-2)*log(Aeta) + log(Aeta^alpha + (1-alpha)*eta^(alpha-1)) - logS + logf)
  
  }
  ifelse(op==1,return(sum(llike)),return(llike))
}


