library(ggplot2)

#1D emulation
###########################
franke2d <- function(xx)
{
  ##########################################################################
  #
  # FRANKE'S FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUT:
  #
  # xx = c(x1, x2)
  #
  ##########################################################################
  
  x1 <- xx[1]
  x2 <- xx[2]
  
  term1 <- 0.75 * exp(-(9*x1-2)^2/4 - (9*x2-2)^2/4)
  term2 <- 0.75 * exp(-(9*x1+1)^2/49 - (9*x2+1)/10)
  term3 <- 0.5 * exp(-(9*x1-7)^2/4 - (9*x2-3)^2/4)
  term4 <- -0.2 * exp(-(9*x1-4)^2 - (9*x2-7)^2)
  
  y <- term1 + term2 + term3 + term4
  return(y)
}

franke2d_x1 <- function(xx)
{
  x1 <- xx
  x2 <- 0.3
  
  term1 <- 0.75 * exp(-(9*x1-2)^2/4 - (9*x2-2)^2/4)
  term2 <- 0.75 * exp(-(9*x1+1)^2/49 - (9*x2+1)/10)
  term3 <- 0.5 * exp(-(9*x1-7)^2/4 - (9*x2-3)^2/4)
  term4 <- -0.2 * exp(-(9*x1-4)^2 - (9*x2-7)^2)
  
  y <- term1 + term2 + term3 + term4
  return(y)
}

#A brief exploration of data points for interesting emulation

x <- seq(0,1, by = 0.001)

plot(y = franke2d_x1(x), x = x)

simple_BL_emulator_v1 <- function(x,              # the emulator prediction point
                                  xD,             # the run input locations xD
                                  D,              # the run outputs D = (f(x^1),...,f(x^n))
                                  theta = 1,      # the correlation lengths
                                  sigma = 1,      # the prior SD sigma sqrt(Var[f(x)])
                                  E_f = 0         # prior expectation of f: E(f(x)) = 0 
){
  
  # store length of runs D  
  n <- length(D)
  
  ### Define Covariance structure of f(x): Cov[f(x),f(xdash)] ###
  Cov_fx_fxdash <- function(x,xdash) sigma^2 * exp(-(x-xdash)^2/theta^2) 
  
  
  ### Define 5 objects needed for BL adjustment ###
  # Create E[D] vector
  E_D <- rep(E_f,n)
  
  # Create Var_D matrix:
  Var_D <- matrix(0,nrow=n,ncol=n)
  for(i in 1:n) for(j in 1:n) Var_D[i,j] <- Cov_fx_fxdash(xD[i],xD[j])
  
  # Create E[f(x)]
  E_fx <- E_f
  
  # Create Var_f(x) 
  Var_fx <- sigma^2
  
  # Create Cov_fx_D row vector
  Cov_fx_D <- matrix(0,nrow=1,ncol=n)
  for(j in 1:n) Cov_fx_D[1,j] <- Cov_fx_fxdash(x,xD[j])
  
  
  ### Perform Bayes Linear adjustment to find Adjusted Expectation and Variance of f(x) ###
  ED_fx   <-  E_fx + Cov_fx_D %*% solve(Var_D) %*% (D - E_D)   
  VarD_fx <-  Var_fx - Cov_fx_D %*% solve(Var_D) %*% t(Cov_fx_D)  
  
  ### return emulator expectation and variance ###
  return(c("ExpD_f(x)"=ED_fx,"VarD_f(x)"=VarD_fx))  
  
}

plot_BL_emulator_V1 <- function(
    em_out,         # a nP x 2 matrix of emulator outputs
    xP,             # vector of nP inputs where emulator evaluated
    xD,             # the run input locations
    D,              # the run outputs D = (f(x^1),...,f(x^n))
    maintitle=NULL  # title of the plot, blank if NULL
){
  
  ### plot emulator output ###
  plot(xP,em_out[,"ExpD_f(x)"],ylim=c(0.2,1.4),ty="l",col="blue",lwd=2.5,
       xlab="Input parameter x",ylab="Output f(x)",main=maintitle) # main argument: plot title
  lines(xP,em_out[,"ExpD_f(x)"]+3*sqrt(em_out[,"VarD_f(x)"]),col="red",lwd=2.5)
  lines(xP,em_out[,"ExpD_f(x)"]-3*sqrt(em_out[,"VarD_f(x)"]),col="red",lwd=2.5)
  
  ### plot true function: we would not normally be able to do this! ###
  lines(xP,franke2d_x1(xP),lwd=2,lty=1)
  
  ### Plot the runs ###
  points(xD,D,pch=21,col=1,bg="green",cex=1.5)
  legend('topright',legend=c("Emulator Expectation","Emulator Prediction Interval",
                             "True function f(x)","Model Evaluations"),
         lty=c(1,1,1,NA),pch=c(NA,NA,NA,16),col=c("blue","red",1,"green"),lwd=2.5,pt.cex=1.3)
}

#Now we produce emulation of the function at an interesting value

par(mfrow = c(1,2))

mean(franke2d_x1(seq(0,1, len = 1000)))


xP <- seq(-0.001, 1.001, len = 1000)
nP <- length(xP)
xD <- c(0, 0.25, 0.5, 0.75, 1)
D <- franke2d_x1(xD)

em_out <- t(sapply(xP, simple_BL_emulator_v1,xD=xD, D=D, theta = 0.3, sigma = 0.5, E_f = 0))

plot_BL_emulator_V1(em_out = em_out, xP = xP, xD = xD, D = D, maintitle = 'Emulator Output: 5 runs')

#7 runs

xD <- c(0, 0.125, 0.25 ,0.5, 0.75, 0.875, 1)
D <- franke2d_x1(xD)

em_out <- t(sapply(xP, simple_BL_emulator_v1,xD=xD, D=D, theta = 0.3, sigma = 0.5, E_f = 0))

plot_BL_emulator_V1(em_out = em_out, xP = xP, xD = xD, D = D, maintitle = 'Emulator Output: 7 runs')

#Now we want to showcase a variety of theta, sigma, and E_f runs

theta <- c(0, 0.25, 0.5, 0.75, 1)
sigma <- c(0, 0.5, 1)
E_f <- c()

em_out <- t(sapply(xP, simple_BL_emulator_v1,xD=xD, D=D, theta = 0.01, sigma = 0.5, E_f = 0))

plot_BL_emulator_V1(em_out = em_out, xP = xP, xD = xD, D = D, maintitle = 'Emulator Output: 7 runs')

#It is worth saying we should also do further analysis of the emulator at a later date.
#These involve altering of theta, sigma and E_f
par(mfrow  =c(3,2))
em_out <- t(sapply(xP, simple_BL_emulator_v1,xD=xD, D=D, theta = 0.01, sigma = 0.5, E_f = 0))
plot_BL_emulator_V1(em_out = em_out, xP = xP, xD = xD, D = D, maintitle = 'Emulator Output 1: 5 runs')
em_out <- t(sapply(xP, simple_BL_emulator_v1,xD=xD, D=D, theta = 0.25, sigma = 0.5, E_f = 0))
plot_BL_emulator_V1(em_out = em_out, xP = xP, xD = xD, D = D, maintitle = 'Emulator Output 2: 5 runs')
em_out <- t(sapply(xP, simple_BL_emulator_v1,xD=xD, D=D, theta = 0.5, sigma = 0.5, E_f = 0))
plot_BL_emulator_V1(em_out = em_out, xP = xP, xD = xD, D = D, maintitle = 'Emulator Output 3: 5 runs')
em_out <- t(sapply(xP, simple_BL_emulator_v1,xD=xD, D=D, theta = 1, sigma = 0.5, E_f = 0))
plot_BL_emulator_V1(em_out = em_out, xP = xP, xD = xD, D = D, maintitle = 'Emulator Output 4: 5 runs')
em_out <- t(sapply(xP, simple_BL_emulator_v1,xD=xD, D=D, theta = 0.25, sigma = 0.25, E_f = 0))
plot_BL_emulator_V1(em_out = em_out, xP = xP, xD = xD, D = D, maintitle = 'Emulator Output 5: 5 runs')
em_out <- t(sapply(xP, simple_BL_emulator_v1,xD=xD, D=D, theta = 0.25, sigma = 0.5, E_f = 0))
plot_BL_emulator_V1(em_out = em_out, xP = xP, xD = xD, D = D, maintitle = 'Emulator Output 6: 5 runs')

par(mfrow = c(1,1))

###########################
#2D emulation
###########################

franke2d_x2<- function(xx)
{
  x1 <- xx[,1]
  x2 <- xx[,2]
  
  term1 <- 0.75 * exp(-(9*x1-2)^2/4 - (9*x2-2)^2/4)
  term2 <- 0.75 * exp(-(9*x1+1)^2/49 - (9*x2+1)/10)
  term3 <- 0.5 * exp(-(9*x1-7)^2/4 - (9*x2-3)^2/4)
  term4 <- -0.2 * exp(-(9*x1-4)^2 - (9*x2-7)^2)
  
  y <- term1 + term2 + term3 + term4
  return(y)
}

simple_BL_emulator_v2 <- function(x,              # the emulator prediction point
                                  xD,             # the run input locations xD
                                  D,              # the run outputs D = (f(x^1),...,f(x^n))
                                  theta = 1,      # the correlation lengths
                                  sigma = 1,      # the prior SD sigma sqrt(Var[f(x)])
                                  E_f = 0         # prior expectation of f: E(f(x)) = 0 
){
  
  # store length of runs D  
  n <- length(D)
  
  ### Define Covariance structure of f(x): Cov[f(x),f(xdash)] ###
  # Cov_fx_fxdash <- function(x,xdash) sigma^2 * exp(-(x-xdash)^2/theta^2)    # XXX Old 1D version
  Cov_fx_fxdash <- function(x,xdash) sigma^2 * exp(-sum((x-xdash)^2)/theta^2) # XXX New 2D version
  
  
  ### Define 5 objects needed for BL adjustment ###
  # Create E[D] vector
  E_D <- rep(E_f,n)
  
  # Create Var_D matrix:
  Var_D <- matrix(0,nrow=n,ncol=n)
  # for(i in 1:n) for(j in 1:n) Var_D[i,j] <- Cov_fx_fxdash(xD[i],xD[j])  # XXX Old 1D version
  for(i in 1:n) for(j in 1:n) Var_D[i,j] <- Cov_fx_fxdash(xD[i,],xD[j,])  # XXX New 2D version
  
  # Create E[f(x)]
  E_fx <- E_f
  
  # Create Var_f(x) 
  Var_fx <- sigma^2
  
  # Create Cov_fx_D row vector
  Cov_fx_D <- matrix(0,nrow=1,ncol=n)
  # for(j in 1:n) Cov_fx_D[1,j] <- Cov_fx_fxdash(x,xD[j])    # XXX Old 1D version
  for(j in 1:n) Cov_fx_D[1,j] <- Cov_fx_fxdash(x,xD[j,])    # XXX New 2D version
  
  
  ### Perform Bayes Linear adjustment to find Adjusted Expectation and Variance of f(x) ###
  ED_fx   <-  E_fx + Cov_fx_D %*% solve(Var_D) %*% (D - E_D)   
  VarD_fx <-  Var_fx - Cov_fx_D %*% solve(Var_D) %*% t(Cov_fx_D)  
  
  ### return emulator expectation and variance ###
  return(c("ExpD_f(x)"=ED_fx,"VarD_f(x)"=VarD_fx))  
  
}

#Now we construct data to simulate over

D_grid <- seq(0.05,0.95,len = 4)
xD <- as.matrix(expand.grid('x1' = D_grid, 'x2' = D_grid))

D <- franke2d_x2(xD)

x_grid <- seq(-0.001, 1.001, len = 50)
xP <- as.matrix(expand.grid('x1' = x_grid, 'x2' = x_grid))
dim(xD)

em_out <- t(apply(xP, 1, simple_BL_emulator_v2, xD = xD, D = D, theta = 0.35, sigma = 1, E_f = 0))

E_D_franke_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
Var_D_franke_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid))

filled.contour(x_grid, x_grid, E_D_franke_mat, xlab = 'x1', ylab = 'x2')

filled.contour(x_grid, x_grid, E_D_franke_mat, xlab = 'x1', ylab = 'x2',
               plot.axes = {
                 axis(1);axis(2)
                 points(xD, pch = 21, col = 1, bg = 'green', cex = 1.5)
                 }
               )

emul_fill_cont <- function(
    cont_mat,            # matrix of values we want contour plot of 
    cont_levs=NULL,      # contour levels (NULL: automatic selection)
    nlev=20,             # approx no. of contour levels for auto select  
    plot_xD=TRUE,        # plot the design runs TRUE or FALSE
    xD=NULL,             # the design points if needed
    xD_col="green",      # colour of design runs
    x_grid,              # grid edge locations that define xP
    ...                  # extra arguments passed to filled.contour
){
  
  ### Define contour levels if necessary ###
  if(is.null(cont_levs)) cont_levs <- pretty(cont_mat,n=nlev)     
  
  ### create the filled contour plot ###
  filled.contour(x_grid,x_grid,cont_mat,levels=cont_levs,xlab="x1",ylab="x2",...,  
                 plot.axes={axis(1);axis(2)                 # sets up plotting in contour box
                   contour(x_grid,x_grid,cont_mat,add=TRUE,levels=cont_levs,lwd=0.8)   # plot contour lines
                   if(plot_xD) points(xD,pch=21,col=1,bg=xD_col,cex=1.5)})  # plot design points
}

library(viridisLite)

exp_cols <- magma
diag_cols <- turbo

par(mfrow = c(1,2))

emul_fill_cont(E_D_franke_mat, seq(-2, 2, 0.05), xD = xD, x_grid = x_grid,
               color.palette = exp_cols,
               main = 'Franke Emulator Adjusted Expectation')

emul_fill_cont(Var_D_franke_mat, NULL, xD = xD, x_grid = x_grid,
               main = 'Franke Emulator Adjusted Variance')

franke_mat <- matrix(franke2d_x2(xP), nrow = length(x_grid), ncol = length(x_grid))

emul_fill_cont(franke_mat, seq(-2, 2, 0.05), xD = xD, x_grid = x_grid,
               color.palette = exp_cols,
               main = 'Franke True Computer Model')

S_diag_mat <- (E_D_franke_mat - franke_mat) / sqrt(Var_D_franke_mat)

emul_fill_cont(S_diag_mat, cont_levels = seq(-3, 3, 0.25), xD = xD, x_grid = x_grid,
               xD_col = 'purple',
               color.palette = diag_cols,
               main = 'Franke Emulator Diagnostics')

#Now we want to evaluate 2D emulator for different values of certain statistics

theta_seq <- c(0.05,0.1,0.15,0.2,0.35,0.4,0.45)

sigma_seq <- 

### loop over vector of theta values ###
for(i in 1:length(theta_seq)){
  
  ### Evaluate emulator over 201 prediction points xP and store in matrices ###
  em_out <- t(apply(xP,1,simple_BL_emulator_v2,xD=xD,D=D,theta=theta_seq[i],sigma=1,E_f=0))   
  E_D_franke_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  
  ### Plot filled contour plot of emulator expectation ###
  emul_fill_cont(cont_mat=E_D_franke_mat,cont_levs=seq(-2,2,0.05),xD=xD,x_grid=x_grid,
                 color.palette=magma,main=paste("Emul. Adjusted Expectation E_D[f(x)], theta =",theta_seq[i]))
}

for(i in 1:length(theta_seq)){
  
  ### Evaluate emulator over 201 prediction points xP and store in matrices ###
  em_out <- t(apply(xP,1,simple_BL_emulator_v2,xD=xD,D=D,theta=theta_seq[i],sigma=1,E_f=0))   
  Var_D_franke_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  
  ### Plot filled contour plot of emulator expectation ###
  emul_fill_cont(cont_mat=Var_D_franke_mat,nlev=12,xD=xD,x_grid=x_grid,
                 main=paste("Emul. Adjusted Variance Var_D[f(x)], theta =",theta_seq[i]))
}

#Now start fresh with LHD

lhd_maximin <- function(nl=16){                    # nl = number of points in LHD 
  
  x_lhd <- cbind("x1"=sample(0:(nl-1)),"x2"=sample(0:(nl-1))) / nl  +  0.5/nl  # create LHD
  
  ### Maximin loop: performs swaps on 1st of two closest points with another random point
  for(i in 1:1000){
    mat <- as.matrix(dist(x_lhd)) + diag(10,nl) # creates matrix of distances between points
    # note the inflated diagonal 
    closest_runs <- which(mat==min(mat),arr.ind=TRUE)   # finds pairs of closest runs
    ind <- closest_runs[sample(nrow(closest_runs),1),1] # chooses one of close runs at random
    swap_ind <- sample(setdiff(1:nl,ind),1)       # randomly selects another run to swap with
    x_lhd2 <- x_lhd                               # creates second version of LHD
    x_lhd2[ind[1],1]   <- x_lhd[swap_ind,1] # swaps x_1 vals between 1st close run & other run
    x_lhd2[swap_ind,1] <- x_lhd[ind[1],1]   # swaps x_1 vals between 1st close run & other run
    if(min(dist(x_lhd2)) >= min(dist(x_lhd))-0.00001) {  # if min distance between points is same or better
      x_lhd <- x_lhd2                                    # we replace LHD with new LHD with the swap
      # cat("min dist =",min(dist(x_lhd)),"Iteration = ",i,"\n") # write out min dist 
    }
  }
  
  ### plot maximin LHD ###
  plot(x_lhd,xlim=c(0,1),ylim=c(0,1),pch=16,xaxs="i",yaxs="i",col="blue",
       xlab="x1",ylab="x2",cex=1.4)
  abline(h=(0:nl)/nl,col="grey60")
  abline(v=(0:nl)/nl,col="grey60")
  return(x_lhd)
}

set.seed(1)
x_lhd <- lhd_maximin(16)

xD <- x_lhd

D <- franke2d_x2(xD)

em_out <- t(apply(xP, 1, simple_BL_emulator_v2, xD = xD, D = D, theta = 0.35, sigma = 1, E_f = 0))
E_D_franke_mat <- matrix(em_out[,'ExpD_f(x)'], nrow = length(x_grid), ncol = length(x_grid))
Var_D_franke_mat <- matrix(em_out[, 'VarD_f(x)'], nrow = length(x_grid), ncol = length(x_grid))

franke_mat <- matrix(franke2d_x2(xP), nrow = length(x_grid), ncol = length(x_grid))

S_diag_mat <- (E_D_franke_mat - franke_mat) / sqrt(Var_D_franke_mat)

emul_fill_cont(E_D_franke_mat, seq(-2, 2, 0.05), xD = xD, x_grid = x_grid,
               color.palette = exp_cols, main = 'Franke Emulator Adjusted Expectation')

emul_fill_cont(franke_mat, seq(-2, 2, 0.05), xD = xD, x_grid = x_grid,
               color.palette = exp_cols, main = 'Franke True Computer Model')

emul_fill_cont(Var_D_franke_mat, NULL, xD = xD, x_grid = x_grid,
               main = 'Franke Emulator Adjuster Variance')

emul_fill_cont(S_diag_mat, cont_levels = seq(-3, 3, 0.25), xD = xD, x_grid = x_grid,
               xD_col = 'purple', color.palette = diag_cols, main = 'Franke Emulator Diagnostics')


theta_seq <- c(0.05,0.1,0.15,0.25,0.35,0.45)

### loop over vector of theta values ###
for(i in 1:length(theta_seq)){
  
  ### Evaluate emulator over 201 prediction points xP and store in matrices ###
  em_out <- t(apply(xP,1,simple_BL_emulator_v2,xD=xD,D=D,theta=theta_seq[i],sigma=1,E_f=0))   
  E_D_franke_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  
  ### Plot filled contour plot of emulator expectation ###
  emul_fill_cont(cont_mat=E_D_franke_mat,cont_levs=seq(-2,2,0.05),xD=xD,x_grid=x_grid,
                 color.palette=magma,main=paste("Franke Adjusted Expectation, theta =",theta_seq[i]))
}

for(i in 1:length(theta_seq)){
  
  ### Evaluate emulator over 201 prediction points xP and store in matrices ###
  em_out <- t(apply(xP,1,simple_BL_emulator_v2,xD=xD,D=D,theta=theta_seq[i],sigma=1,E_f=0))   
  Var_D_franke_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  
  ### Plot filled contour plot of emulator expectation ###
  emul_fill_cont(cont_mat=Var_D_franke_mat,nlev=12,xD=xD,x_grid=x_grid,
                 main=paste("Franke Adjusted Variance, theta =",theta_seq[i]))
}

###########################
#History Matching
###########################
emul_fill_cont_V2 <- function(
    cont_mat,            # matrix of values we want contour plot of 
    cont_levs=NULL,      # contour levels (NULL: automatic selection)
    cont_levs_lines=NULL,   # contour levels for lines (NULL: automatic selection)
    nlev=20,             # approx no. of contour levels for auto select  
    plot_xD=TRUE,        # plot the design runs TRUE or FALSE
    xD=NULL,             # the design points if needed
    xD_col="green",      # colour of design runs
    x_grid,              # grid edge locations that define xP
    ...                  # extra arguments passed to filled.contour
){
  
  ### Define contour levels if necessary ###
  if(is.null(cont_levs)) cont_levs <- pretty(cont_mat,n=nlev)
  
  ### create the filled contour plot ###
  filled.contour(x_grid,x_grid,cont_mat,levels=cont_levs,xlab="x1",ylab="x2",...,  
                 plot.axes={axis(1);axis(2)                 # sets up plotting in contour box
                   if(is.null(cont_levs_lines)) contour(x_grid,x_grid,cont_mat,add=TRUE,levels=cont_levs,lwd=0.8) # plot usual contour lines 
                   if(!is.null(cont_levs_lines)) {
                     contour(x_grid,x_grid,cont_mat,add=TRUE,levels=cont_levs,lwd=0.4,labels="")   # plot thin contour lines 
                     contour(x_grid,x_grid,cont_mat,add=TRUE,levels=cont_levs_lines,lwd=2)   # plot thick contour lines
                   }
                   if(plot_xD) points(xD,pch=21,col=1,bg=xD_col,cex=1.5)})  # plot design points
}


x_grid <- seq(-0.001, 1.001, len = 50)
xP <- as.matrix(expand.grid('x1' = x_grid, 'x2' = x_grid))

z <- 0.8
sigma_e <- 0.15
sigma_epsilon <- 0.15

set.seed(1)
xD_w1 <- lhd_maximin(16)

D <- franke2d_x2(xD_w1)

em_out <- t(apply(xP, 1, simple_BL_emulator_v2, xD = xD_w1, D = D, theta = 0.4, sigma = 1, E_f = 0))
E_D_franke_mat <- matrix(em_out[,'ExpD_f(x)'], nrow = length(x_grid), ncol = length(x_grid))
Var_D_franke_mat <- matrix(em_out[,'VarD_f(x)'], nrow = length(x_grid), ncol = length(x_grid))

Imp_mat <- sqrt((E_D_franke_mat - z)^2 / (Var_D_franke_mat + sigma_e^2 + sigma_epsilon^2))

imp_cols <- function(n) turbo(n,begin=0.15,end=1)
imp_levs <- c(0,seq(1,2.75,0.25),seq(3,18,2),20,30,41)

emul_fill_cont_V2(cont_mat=Imp_mat,cont_levs=imp_levs,cont_levs_lines=3,xD=xD_w1,x_grid=x_grid,
                  xD_col="purple",color.palette=imp_cols,main="Implausibility I(x): Wave 1")

#Now we want to implement the following points to fully explore the area

xD_w2 <- matrix(
  c(0.78, 0.25,
    0.95, 0.05,
    0.5, 0.2, 
    0.05, 0.05,
    0.37, 0.05,
    0.05, 0.55,
    0.05, 0.95, 
    0.5, 0.95,
    0.95, 0.95,
    0.4, 0.3
    ), ncol = 2, byrow = TRUE
)

emul_fill_cont_V2(cont_mat=Imp_mat,cont_levs=imp_levs,cont_levs_lines=3,xD=rbind(xD,xD_w2),
                  x_grid=x_grid,xD_col=rep(c("purple","pink"),c(nrow(xD),nrow(xD_w2))),
                  color.palette=imp_cols,main="Implausibility I(x): Wave 1")

xD_w3 <- rbind(xD, xD_w2)

for(k in 0:10){                          # k=0: wave 1, k>0 add k wave 2 runs sequentially
  
  xD <- xD_w1                           # the 14 point wave 1 design
  if(k>0) xD <- rbind(xD,xD_w2[1:k,])   # k=0: wave 1, k>0 add wave 2 runs sequentially
  
  ### Perform 14 + k runs of model and store as D (would take days for realistic example!) ###
  D <- franke2d_x2(xD)
  
  ### Evaluate emulator over 50x50=2500 prediction points xP and store as matrices ###
  em_out <- t(apply(xP,1,simple_BL_emulator_v2,xD=xD,D=D,theta=0.4,sigma=1,E_f=0))   
  E_D_franke_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  Var_D_franke_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  
  ### Evaluate true function and store in matrix for diagnostic comparisons ###
  frankexP_mat <- matrix(franke2d_x2(xP),nrow=length(x_grid),ncol=length(x_grid)) 
  
  ### Calculate Implausibility Measure Over All 50x50 = 2500 input points in xP ###
  Imp_mat <- sqrt( (E_D_franke_mat - z)^2 / (Var_D_franke_mat + sigma_e^2 + sigma_epsilon^2) ) 
  
  ### Calculate Imp Measure for True f(x) Over All 50x50 = 2500 input points in xP ###
  Imp_true_mat <- sqrt( (frankexP_mat - z)^2 / (sigma_e^2 + sigma_epsilon^2) ) 
  
  ### Define colours and levels for implausibility plots ###
  imp_cols <- function(n) turbo(n,begin=0.15,end=1)
  imp_levs <- c(0,seq(1,2.75,0.25),seq(3,18,2),20,30,41)
  
  ### if k=0 plot wave 1 runs only ###
  if(k==0) emul_fill_cont_V2(cont_mat=Imp_mat,cont_levs=imp_levs,cont_levs_lines=3,xD=xD,
                             x_grid=x_grid,xD_col="purple",color.palette=imp_cols,
                             main="Implausibility I(x): Wave 1")
  ### plot current runs in purple and remaining unevaluated wave 2 runs in pink ###
  emul_fill_cont_V2(cont_mat=Imp_mat,cont_levs=imp_levs,cont_levs_lines=3,xD=rbind(xD,xD_w2),
                    x_grid=x_grid,xD_col=rep(c("purple","pink"),c(nrow(xD)+k,nrow(xD_w2)-k)),  # cover unevaluated w2 points in pink
                    color.palette=imp_cols,main="Implausibility I(x): Wave 2")
  ### once last run done so k=8, plot implausibility for true function f(x) to compare ###
  if(k==nrow(xD_w2)) emul_fill_cont_V2(cont_mat=Imp_true_mat,cont_levs=imp_levs,
                                       cont_levs_lines=3,xD=xD,x_grid=x_grid,xD_col="purple",plot_xD=FALSE,
                                       color.palette=imp_cols,main="Implausibility I(x) using True Model")
}

#Hence we have found our implausiblity region for the observed value

###########################
#Efficiency and Optimisation
###########################

library(pdist)

simple_BL_emulator_v3 <- function(xP,             # the set of emulator prediction points
                                  xD,             # the run input locations xD
                                  D,              # the run outputs D = (f(x^1),...,f(x^n))
                                  theta = 1,      # the correlation lengths (can be a vector)
                                  sigma = 1,      # the prior SD sigma sqrt(Var[f(x)])
                                  E_f = 0,         # prior expectation of f: E(f(x)) = 0 
                                  using_pdist = 1  # if you have installed pdist package
){
  
  # store length of runs D and number of prediction points xP
  n <- length(D)
  nP <- nrow(xP)       # XXX New V3
  
  # # Rescale each input by theta. Works for different theta for each input and for same theta
  xP <- t(t(xP)/theta)     # XXX New V3: solution to Exercise 8.3, CP 3.
  xD <- t(t(xD)/theta)     # XXX New V3: solution to Exercise 8.3, CP 3.
  
  ### Define Cov structure of f(x): Cov[f(x),f(xdash)], now to act on matrix of distances ###
  # Cov_fx_fxdash <- function(x,xdash) sig^2 * exp(-sum((x-xdash)^2)/theta^2) # XXX Old V2
  Cov_fx_fxdash <- function(dist_matrix) sigma^2 * exp(-(dist_matrix)^2) # XXX New dist V3
  
  
  ### Define 5 objects needed for BL adjustment ###
  # Create E[D] vector
  E_D <- rep(E_f,n)
  
  # Create Var_D matrix:
  # Var_D <- matrix(0,nrow=n,ncol=n)                                        # XXX Old V2
  # for(i in 1:n) for(j in 1:n) Var_D[i,j] <- Cov_fx_fxdash(xD[i,],xD[j,])  # XXX Old V2
  Var_D <- Cov_fx_fxdash( as.matrix(dist(xD)) )                       # XXX New dist V3
  
  # Create E[f(x)]
  E_fx <- rep(E_f,nP)
  
  # Create Var_f(x) 
  Var_fx <- rep(sigma^2,nP)
  
  # Create Cov_fx_D row vector now using pdist() function if available, if not use dist()
  # Cov_fx_D <- matrix(0,nrow=1,ncol=n)                       # XXX Old V2
  # for(j in 1:n) Cov_fx_D[1,j] <- Cov_fx_fxdash(x,xD[j,])    # XXX Old V2
  if(using_pdist)  Cov_fx_D <- Cov_fx_fxdash( as.matrix(pdist(xP,xD)) )   # XXX New V3
  if(!using_pdist) 
    Cov_fx_D <- Cov_fx_fxdash( as.matrix(dist(rbind(xP,xD)))[1:nP,(nP+1):(nP+n)] )# XXX NewV3
  
  # find inverse of Var_D using Cholesky decomposition (check Wikipedia if interested!)
  Var_D_inv <- chol2inv(chol(Var_D))     # more efficient way to calculate inverse of Cov mat
  
  ### Perform Bayes Linear adjustment to find Adjusted Expectation and Variance of f(x) ###
  cov_fx_D_Var_D_inv <- Cov_fx_D %*% Var_D_inv  # Need this twice so pre-calculate here
  ED_fx   <-  E_fx + cov_fx_D_Var_D_inv %*% (D - E_D) # adj. expectation of ALL xP pts at once
  # VarD_fx <-  Var_fx - cov_fx_D_Var_D_inv %*% t(Cov_fx_D)       # XXX Old V2
  VarD_fx   <-  Var_fx - apply(cov_fx_D_Var_D_inv * Cov_fx_D,1,sum) # fast way to get diagonals 
  # and hence all nP variances (Note: does not do full np x np covariance matrix)
  
  ### return emulator expectation and variance ###
  return(cbind("ExpD_f(x)"=c(ED_fx),"VarD_f(x)"=VarD_fx))  
  
}

set.seed(1)
xD_w4 <- lhd_maximin(16)

x_grid <- seq(-0.001,1.001,len=50)
xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
xD <- xD_w4
xD <- rbind(xD, c(0.05, 0.05), c(0.2, 0.03))

D <- franke2d_x2(xD)

em_out <- simple_BL_emulator_v3(xP = xP, xD = xD, D = D, theta = c(0.45, 0.45), sigma = 1, E_f = 0)
E_D_franke_mat <- matrix(em_out[,'ExpD_f(x)'], nrow = length(x_grid), ncol = length(x_grid))
Var_D_franke_mat <- matrix(em_out[,'VarD_f(x)'], nrow = length(x_grid), ncol = length(x_grid))

franke_plus <- max(D)
sigma_epsilon <- 0.05

frankexP_mat <- matrix(franke2d_x2(xP), nrow = length(x_grid), ncol = length(x_grid))
franke_plus_true <- max(frankexP_mat)

emul_fill_cont(cont_mat = E_D_franke_mat, cont_levs = NULL, xD = xD, x_grid = x_grid,
              color.palette = exp_cols, main = 'Franke Emulator Adjusted Expectation')

emul_fill_cont(cont_mat=frankexP_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
               color.palette=exp_cols,plot_xD=FALSE,main="True Computer Franke Model")

#Implausibility for optimisation

### Calculate new Implausibility Measure Over All 50x50 = 2500 input points in xP ###
Imp_mat <-  (franke_plus-E_D_franke_mat) / sqrt(Var_D_franke_mat + sigma_epsilon^2) 
Imp_mat[Imp_mat<0] <- 0.0001   # replace negatives with small value for plotting purposes

### Calculate true Implausibility Measure Over All 201 input points in xP ###
Imp_true_mat <- (franke_plus_true-frankexP_mat) / sqrt(sigma_epsilon^2 ) 
Imp_true_mat[Imp_true_mat>40] <- 40  # limit high implausibilities for plotting purposes

imp_cols <- function(n) turbo(n,begin=0.15,end=1)
imp_levs <- c(0,seq(1,2.75,0.25),seq(3,18,2),20,30,41) 

### Plot New Optimisation Implausibility ###
emul_fill_cont_V2(cont_mat=Imp_mat,cont_levs=imp_levs,cont_levs_lines=3,
                  xD=xD,x_grid=x_grid,xD_col="purple",color.palette=imp_cols,
                  main="Max Implausibility I_M(x) Wave 2")

emul_fill_cont_V2(cont_mat=Imp_true_mat,cont_levs=imp_levs,cont_levs_lines=3,xD=xD,x_grid=x_grid,
                  xD_col="purple",plot_xD=FALSE,color.palette=imp_cols,
                  main="Implausibility I(x) using Franke")

