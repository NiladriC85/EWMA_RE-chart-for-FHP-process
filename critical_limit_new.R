rm(list=ls())

######### Required packages ########
library("parallel",quietly = T)
library("iterators",quietly = T)
library("foreach",quietly = T)
library("doParallel",quietly = T)
require("clue", quietly=T)
require("energy", quietly = T)
require("randtoolbox", quietly = T)
require("pracma", quietly = T)
require("kernlab", quietly = T)
require("crossmatch", quietly = T)
require("HHG", quietly = T)
require("gTests", quietly = T)
library("MASS",quietly = T)
library("class",quietly = T)
library("robustbase",quietly = T)
library("sfsmisc",quietly = T)
library("geometry",quietly = T)
library("ddalpha",quietly = T)
library("mvtnorm",quietly = T)
library("writexl",quietly = T)
library("readxl",quietly = T)
library("purrr",quietly=T)
library("GA",quietly=T)

find_L <- function(m, n, p, lambda, rho, B, sim, FAP0, aux_sims = 5000, tol = 1e-7, max.iter = 100) {
  mu0 <<- rep(0, p)
  sigma <<- rho * outer(seq_len(p), seq_len(p), function(i, j) pmin(i, j) / pmax(i, j))
  lambda <<- lambda
  B <<- B
  sim <<- sim
  
  compute_RE=function(x,y,m=nrow(x),n=nrow(y),dim=ncol(x),gridch=torus(m+n,dim))
  {
    comdata=rbind(x,y)
    distmat=matrix(0,nrow=m+n,ncol=m+n)
    for(i in 1:(m+n))
      distmat[i,]=apply((comdata[i,]-t(gridch)),2,Norm)^2
    assignmentFUN=solve_LSAP(distmat)
    assignmentSOL=cbind(seq_along(assignmentFUN),assignmentFUN)
    randenergySTAT=eqdist.etest(gridch[assignmentSOL[,2],],sizes = c(m,n), R=1)
    return(randenergySTAT$statistic)
  }
  
  t_stat <- function() {
    s1 <- mvrnorm(m, mu0, sigma)
    s2 <- mvrnorm(n, mu0, sigma)
    grid <- torus(m + n, p)
    as.numeric(compute_RE(s1, s2))
  }
  
  RE_auxiliary <- replicate(aux_sims, t_stat())
  mean_RE <<- mean(RE_auxiliary)
  sd_RE <<- sd(RE_auxiliary)
  
  maximum_EWMA <- function(B, s1, grid) {
    RE <- numeric(B)
    EWMA <- numeric(B)
    s2 <- mvrnorm(n, mu0, sigma)
    RE[1] <- as.numeric(compute_RE(s1, s2, grid))
    EWMA[1] <- lambda * RE[1] + (1 - lambda) * mean_RE
    for (k in 2:B) {
      s2 <- mvrnorm(n, mu0, sigma)
      RE[k] <- as.numeric(compute_RE(s1, s2))
      EWMA[k] <- lambda * RE[k] + (1 - lambda) * EWMA[k - 1]
    }
    max(EWMA)
  }
  
  conditional_FAP <<- function(L, m, n, p, rho) {
    ucl <- mean_RE + L * sqrt(lambda / (2 - lambda)) * sd_RE
    s1 <- mvrnorm(m, mu0, sigma)
    grid <- torus(m + n, p)
    cm <- replicate(sim, maximum_EWMA(B, s1, grid))
    mean(cm > ucl)
  }
  
  objective_function1 <- function(Lval) {
    as.numeric(quantile(replicate(sim, conditional_FAP(Lval, m, n, p, rho)), 0.95) - FAP0)
  }
  
  bisection <- function(f, a, b, tol = tol, max.iter = max.iter) {
    fa <- f(a); fb <- f(b)
    while (fa * fb >= 0) {
      if (fb < 0) a <- a - 1 else b <- b + 1
      fa <- f(a); fb <- f(b)
    }
    for (i in 1:max.iter) {
      c <- (a + b) / 2
      fc <- f(c)
      if (abs(fc) < tol || (b - a) / 2 < tol) return(c)
      if (fc * fa > 0) { a <- c; fa <- fc } else { b <- c }
    }
    stop("No convergence")
  }
  
  bisection(objective_function1, 1, 3)
}

m <- 50; n <- 5; p <- 3; lambda <- 0.5; rho <- 3; B <- 10; sim <- 10000; FAP0 <- 0.1
L_star <- find_L(m, n, p, lambda, rho, B, sim, FAP0)
L_star

lambda <- 0.5
quantile(replicate(sim, conditional_FAP(3.084, m, n, p, rho)), 0.95)



