
#seed = as.numeric(commandArgs(trailingOnly = TRUE))
seed = 1

set.seed(seed)

n <- 1000 # sample size

pmax = 30

U = runif(n)
X = matrix(runif(n*pmax),ncol=pmax)


source("Models.R")

for(model in names(models)){
  name = models[[model]][["name"]]
  p = models[[model]][["p"]] ; q = models[[model]][["q"]]
  B_0 = models[[model]][["B_0"]] ; B_1 = models[[model]][["B_1"]]
  xi = models[[model]][["xi"]] ; ell = models[[model]][["models"]]
  nameFile = paste(
    "Data-",name,
    "-p_",p,
    "-q_",q,
    "-n_",n,
    ".RData",
    sep=""
  )
  

  
  y <- U^(-apply(t(B_0)%*%t(X),2,xi))*apply(cbind(U,t(B_1)%*%t(X)),2,ell)
  
  

}