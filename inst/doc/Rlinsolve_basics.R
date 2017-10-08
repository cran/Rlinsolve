## ---- message=F, warning=F-----------------------------------------------
library(microbenchmark)
library(pcg)
library(optR)
library(Rlinsolve)

## ------------------------------------------------------------------------
A = matrix(rnorm(10000*50),nrow=10000)
x = matrix(rnorm(50))
b = A%*%x

## ------------------------------------------------------------------------
# for SOLVE in R base, it needs to be transformed into normal equation form.
Anormal = t(A)%*%A
bnormal = t(A)%*%b
microbenchmark(solve(Anormal,bnormal))

## ------------------------------------------------------------------------
microbenchmark(pcg(Anormal,bnormal),
               optR(A,b,method="gauss"),
               times=10)

## ------------------------------------------------------------------------
microbenchmark(lsolve.sor(A,b,verbose=FALSE),
               lsolve.bicgstab(A,b,verbose=FALSE),
               times=10)

## ------------------------------------------------------------------------
Psparse = aux.fisch(10,sparse=TRUE) # sparse matrix
Pdense = aux.fisch(10,sparse=FALSE) # dense  matrix
x = matrix(rnorm(100))
b = Pdense %*% x

## ------------------------------------------------------------------------
microbenchmark(solve(Pdense,b),
               optR(Pdense,b,method="gauss"),
               pcg(Pdense,b),
               times=10)
microbenchmark(lsolve.bicg(Psparse,b,verbose=FALSE),
               lsolve.cgs(Psparse,b,verbose=FALSE),
               lsolve.gs(Psparse,b,verbose=FALSE),
               lsolve.sor(Psparse,b,verbose=FALSE),
               times=10)

