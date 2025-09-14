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
library("purrr",quietly=T)
library(ggplot2)
library(VIM)
library(glmnet)
library(lattice)
library(caret)
library(ggrepel)
library(dplyr)
library(tidyr)
library(reshape2)
library(gridExtra)
library(rebmix)

##############################################################################
####### the random vectors should be row vector ##################################
######### Computing Rank Energy Statistic #########
computestatistic=function(x,y,m=nrow(x),n=nrow(y),dim=ncol(x),gridch=torus(m+n,dim))
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
####################################################################################
#################################################################
data(sensorlessdrive)
data= sensorlessdrive


data_IC = data[data$Class==1,]
data_OC = data[data$Class==9,]

##-------------------------------Correlation plot for the data ----------------------------------------------
data = cbind(data_IC$V5,data_IC$V7,data_IC$V11)

# Convert the matrix to a data frame
data <- as.data.frame(data)

library(GGally)
ggpairs(data)+theme_classic()

library(PerformanceAnalytics)
chart.Correlation(data,histogram = T,pch = 19)

##--------------------------------Marginal density plot for the IC data------------------------------------------
# Add a row index as a column
data$sample <- 1:nrow(data)

# Convert the data frame to long format using melt()
long_data <- melt(data, id.vars = "sample", variable.name = "variable", value.name = "value")

# Create a distribution plot for each variable
plots <- lapply(unique(long_data$variable), function(var) {
  ggplot(long_data[long_data$variable == var, ], aes(x = value)) +
    geom_density(fill = "lightblue", color = "darkblue") +
    labs(x = var, y = "Density") +
    theme_classic()+
    theme(
      panel.background = element_rect(fill = "seashell", color = NA),
      legend.position = "none"  # Remove legend for simplicity
    )
})

plots[[1]]
plots[[2]]
plots[[3]]

# Arrange and display the plots
grid.arrange(grobs = plots, nrow = 2)
grid.arrange(
  plots[[1]], plots[[2]],          # First row with two plots
  plots[[3]],                      # Second row with the third plot in the center
  layout_matrix = rbind(c(1, NA,2),    # Define layout: first row with two plots side by side
                        c(NA, 3, NA))   # Second row with the third plot centered
)

##-------------------------------------------------------------------------------------------------------

dim(data_IC)
dim(data_OC)

ref_data = matrix(unlist(data_IC[1:50,-4]),ncol = 3)
test_data = matrix(unlist(data_OC[1:50,-4]),ncol = 3)

as.numeric(unlist(computestatistic(ref_data,test_data)))



##--------------- Normality test -------------------------------------------------------
shapiro.test(ref_data[,1])
shapiro.test(ref_data[,2])
shapiro.test(ref_data[,3])
##--------------------------------------------------------------------------------------
##########################################################################################
###############################################################################################################

set.seed(100) 

####################################################################################
m = 50  ## Reference sample size
n = 5  ## Test sample size
p = 3   ## dimension of the multivariate data
lambda = 0.1  ## smoothing constant

#L = 10       ## The constant for L*sigma limit 
rho = 3    ## correlation coefficient for the correlation matrix
B = 10       ## time length of the finite horizon process 
sim = 100    ## number of repetitions of the conditional FAP 
FAP0 = 0.1   ## standard FAP
#######################################################################################
###### Parameters for the multivariate normal distribution under H0 ###################

# set mean vector and covariance matrix
mu0 = rep(0,p)
sigma = matrix(NA,nrow = p,ncol = p)
for(i in 1:p){
  for(j in 1:p){
    sigma[i,j] = rho*(min(i,j)/max(i,j))
  }
}

###############################################################################
######################## Calculate the mean of the RE statistic #################################################

t.stat = function(cnt){
  # set mean vector and covariance matrix
  mu0 = rep(0,p)
  sigma = matrix(NA,nrow = p,ncol = p)
  for(i in 1:p){
    for(j in 1:p){
      sigma[i,j] = rho*(min(i,j)/max(i,j))
    }
  }
  # generate random sample
  s1 = mvrnorm(m, mu = mu0, Sigma = sigma)
  s2 = mvrnorm(n, mu = mu0, Sigma = sigma)
  as.numeric(unlist(computestatistic(s1,s2)))
}

## RE test is distribution-free
## So we can estimate the mean of the auxiliary RE statistic for given m,n,p
RE_auxiliary = map_dbl(c(1:10000),~t.stat(.x)) 
quantile(RE_auxiliary,0.99)
mean_RE = mean(RE_auxiliary)
sd_RE = sd(RE_auxiliary)
#######################################################################
RE = numeric(10)
EWMA = numeric(10)
s1 = as.matrix(ref_data)
s2 = as.matrix(test_data[1:5,])  

RE[1] = as.numeric(unlist(computestatistic(s1,s2)))
EWMA[1] = lambda*RE[1]+(1-lambda)*mean_RE

for(i in 2:10){
  s2 = as.matrix(test_data[(5*i-4):(5*i),])  
  RE[i] = as.numeric(unlist(computestatistic(s1,s2)))
  EWMA[i] = lambda*RE[i]+(1-lambda)*EWMA[i-1]
}

EWMA1=EWMA
ucl1 = mean_RE + 1.592*sqrt(lambda/(2-lambda))*sd_RE
ucl1

#######################################################################################################################
#------------------Runtime metrics
#######################################################################################################################
RE   <- numeric(10)
EWMA <- numeric(10)
runtime <- numeric(10)   # to store runtimes for each RE[i]

s1 <- as.matrix(ref_data)
s2 <- as.matrix(test_data[1:5,])  

# First value
t <- system.time({
  RE[1] <- as.numeric(unlist(computestatistic(s1, s2)))
})
runtime[1] <- t["elapsed"]
EWMA[1] <- lambda * RE[1] + (1 - lambda) * mean_RE

# Loop for i = 2,...,10
for (i in 2:10) {
  s2 <- as.matrix(test_data[(5 * i - 4):(5 * i),])  
  t <- system.time({
    RE[i] <- as.numeric(unlist(computestatistic(s1, s2)))
  })
  runtime[i] <- t["elapsed"]
  EWMA[i] <- lambda * RE[i] + (1 - lambda) * EWMA[i - 1]
}

EWMA1  <- EWMA
runtime  # this vector gives runtime for each RE[i]

total_runtime = sum(runtime)
total_runtime

########################################################################################################################

lambda = 0.5  ## smoothing constant
#######################################################################
RE = numeric(10)
EWMA = numeric(10)
s1 = as.matrix(ref_data)
s2 = as.matrix(test_data[1:5,])  

RE[1] = as.numeric(unlist(computestatistic(s1,s2)))
EWMA[1] = lambda*RE[1]+(1-lambda)*mean_RE

for(i in 2:10){
  s2 = as.matrix(test_data[(5*i-4):(5*i),])  
  RE[i] = as.numeric(unlist(computestatistic(s1,s2)))
  EWMA[i] = lambda*RE[i]+(1-lambda)*EWMA[i-1]
}

EWMA2=EWMA
ucl2 = mean_RE + 3.084*sqrt(lambda/(2-lambda))*sd_RE
ucl2
########################################################################################################################
########################################################################################################################

EWMA1
EWMA2

######################################################################################################

# Build plotting data (unchanged)
data <- data.frame(
  Index = rep(seq_along(EWMA1), 2),
  EWMA_Value = c(EWMA1, EWMA2),
  Lambda = factor(rep(c("0.1", "0.5"), each = length(EWMA1)))
)
lambda_labels <- c(expression(lambda == 0.1), expression(lambda == 0.5))

# Thresholds
thr_lo <- 0.770
thr_hi <- 1.171

p <- ggplot(data, aes(x = Index, y = EWMA_Value, color = Lambda, linetype = Lambda)) +
  geom_line(linewidth = 1.2, lineend = "round") +
  geom_hline(yintercept = thr_lo, color = "blue", linetype = "solid", linewidth = 0.7) +
  geom_hline(yintercept = thr_hi, color = "orangered4", linetype = "dashed", linewidth = 0.7) +
  annotate("label", x = Inf, y = thr_lo, label = "",
           hjust = 1.02, vjust = -0.3, size = 3.5, fill = "white", label.size = 0) +
  annotate("label", x = Inf, y = thr_hi, label = "",
           hjust = 1.02, vjust = -0.3, size = 3.5, fill = "white", label.size = 0) +
  labs(x = "Index", y = "EWMA value", color = NULL, linetype = NULL) +
  scale_color_manual(values = c("blue", "orangered4"), labels = lambda_labels) +
  scale_linetype_manual(values = c("solid", "dashed"), labels = lambda_labels) +
  scale_x_continuous(breaks = 1:length(EWMA1), minor_breaks = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    legend.position = "top",
    legend.text = element_text(size = 14),   # << bigger λ labels
    legend.key.width = unit(22, "pt"),
    legend.margin = margin(0, 0, 6, 0)
  )

p <- p +
  guides(
    color    = guide_legend(label.theme = element_text(margin = margin(r = 10))),
    linetype = guide_legend(label.theme = element_text(margin = margin(r = 10)))
  )

# Draw it
p

# High-resolution export (PNG and PDF)
ggsave("EWMA_plot.png", p, width = 7.5, height = 4.2, dpi = 320)
ggsave("EWMA_plot.pdf", p, width = 7.5, height = 4.2, device = cairo_pdf)

########################################################################################################
#---------------- Delayed detection---------------------------------------------------------------------
########################################################################################################

test_data1 = matrix(unlist(rbind(data_IC[51:75,-4],data_OC[1:25,-4])),ncol = 3)

#######################################################################
lambda = 0.1
#######################################################################
RE = numeric(10)
EWMA = numeric(10)
s1 = as.matrix(ref_data)
s2 = as.matrix(test_data1[1:5,])  

RE[1] = as.numeric(unlist(computestatistic(s1,s2)))
EWMA[1] = lambda*RE[1]+(1-lambda)*mean_RE

for(i in 2:10){
  s2 = as.matrix(test_data1[(5*i-4):(5*i),])  
  RE[i] = as.numeric(unlist(computestatistic(s1,s2)))
  EWMA[i] = lambda*RE[i]+(1-lambda)*EWMA[i-1]
}

EWMA1=EWMA
ucl1 = mean_RE + 1.592*sqrt(lambda/(2-lambda))*sd_RE
ucl1
########################################################################################################################

lambda = 0.5  ## smoothing constant
#######################################################################
RE = numeric(10)
EWMA = numeric(10)
s1 = as.matrix(ref_data)
s2 = as.matrix(test_data1[1:5,])  

RE[1] = as.numeric(unlist(computestatistic(s1,s2)))
EWMA[1] = lambda*RE[1]+(1-lambda)*mean_RE

for(i in 2:10){
  s2 = as.matrix(test_data1[(5*i-4):(5*i),])  
  RE[i] = as.numeric(unlist(computestatistic(s1,s2)))
  EWMA[i] = lambda*RE[i]+(1-lambda)*EWMA[i-1]
}

EWMA2=EWMA
ucl2 = mean_RE + 3.084*sqrt(lambda/(2-lambda))*sd_RE
ucl2
########################################################################################################################
########################################################################################################################

EWMA1
EWMA2
#########################################################################################################################

# Build plotting data (unchanged)
data <- data.frame(
  Index = rep(seq_along(EWMA1), 2),
  EWMA_Value = c(EWMA1, EWMA2),
  Lambda = factor(rep(c("0.1", "0.5"), each = length(EWMA1)))
)
lambda_labels <- c(expression(lambda == 0.1), expression(lambda == 0.5))

# Thresholds
thr_lo <- 0.770
thr_hi <- 1.171

p <- ggplot(data, aes(x = Index, y = EWMA_Value, color = Lambda, linetype = Lambda)) +
  geom_line(linewidth = 1.2, lineend = "round") +
  geom_hline(yintercept = thr_lo, color = "blue", linetype = "solid", linewidth = 0.7) +
  geom_hline(yintercept = thr_hi, color = "orangered4", linetype = "dashed", linewidth = 0.7) +
  annotate("label", x = Inf, y = thr_lo, label = "",
           hjust = 1.02, vjust = -0.3, size = 3.5, fill = "white", label.size = 0) +
  annotate("label", x = Inf, y = thr_hi, label = "",
           hjust = 1.02, vjust = -0.3, size = 3.5, fill = "white", label.size = 0) +
  labs(x = "Index", y = "EWMA value", color = NULL, linetype = NULL) +
  scale_color_manual(values = c("blue", "orangered4"), labels = lambda_labels) +
  scale_linetype_manual(values = c("solid", "dashed"), labels = lambda_labels) +
  scale_x_continuous(breaks = 1:length(EWMA1), minor_breaks = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    legend.position = "top",
    legend.text = element_text(size = 14),   # << bigger λ labels
    legend.key.width = unit(22, "pt"),
    legend.margin = margin(0, 0, 6, 0)
  )

p <- p +
  guides(
    color    = guide_legend(label.theme = element_text(margin = margin(r = 10))),
    linetype = guide_legend(label.theme = element_text(margin = margin(r = 10)))
  )

# Draw it
p

# High-resolution export (PNG and PDF)
ggsave("EWMA_plot.png", p, width = 7.5, height = 4.2, dpi = 320)
ggsave("EWMA_plot.pdf", p, width = 7.5, height = 4.2, device = cairo_pdf)


