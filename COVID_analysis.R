######################################1#########################################

#function to calculate k
calculateK <- function(data) {
  n=length(data)
  kSum=0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      kSum=kSum + sum(data[[i]] != data[[j]])
    }
  }
  k_ = kSum / (n * (n - 1) / 2)
  return(k_)
}

#function to calculate W
calculateW <- function(data) {
  n=length(data)
  S=length(data[[1]])
  a1=sum(1 / (1:(n-1)))
  
  w_=S / a1
  return(w_)
}

#function to calculate D
calculateD <- function(data) {
  n=length(data)
  S=length(data[[1]])
  
  k=calculateK(data)
  w=calculateW(data)
  
  a1=sum(1 / (1:n))
  a2=sum(1 / (1:n)^2)
  b1=(n + 1) / (3 * (n - 1))
  b2=2 * (n^2 + n + 3) / (9 * n * (n - 1))
  c1=b1 - 1 / a1
  c2= b2 - (n + 2) / (a1 * n) + a2 / a1^2
  e1=c1 / a1
  e2=c2 / (a1^2 + a2)
  
  d_=(k-w) / sqrt(e1 * S + e2 * S * (S - 1))
  return(d_)
}

################################################################################
###################################2############################################

#observed dataset
observed_datasett=readLines("ms_obs_final.out")
observed_list <- list()
for (i in 1:length(observed_datasett)) {
  observed_dataset = as.numeric(strsplit(observed_datasett[i], "")[[1]])
  observed_list[[i]] = observed_dataset
}
observed_dataset=observed_list
#statistics for observed dataset
oK=calculateK(observed_list)
oK
oW=calculateW(observed_list)
oW
oD=calculateD(observed_list)
oD


sK= vector(length = 10000)
sW= vector(length = 10000)
sD= vector(length = 10000)

simulated_datasett= readLines("ms_sim_final.out")
counter=1
k=1
temp_list = list()
simulated_datasets=vector("list", 10000)

for (line in simulated_datasett) {
  if (line != "") {
    # store the vector data in temp list
    vector_data = as.numeric(strsplit(line, "")[[1]])
    temp_list[[k]] = vector_data
    k=k+1
    
  } else if (!is.null(temp_list)) {
    sK[counter]=calculateK(temp_list)
    sW[counter]=calculateW(temp_list)
    sD[counter]=calculateD(temp_list)
    simulated_datasets[[counter]]=temp_list
    counter=counter+1
    
    temp_list=list()
    k=1
  }
}

#statistics for simulated dataset
sK
sW
sD

################################################################################
####################################3###########################################

#normalising the simulated vectors
normK=(sK - mean(sK)) / sd(sK)
normK
normW=(sW - mean(sW)) / sd(sW)
normW
normD=(sD - mean(sD)) / sd(sD)
normD

#normalising the observed vectors
o_normK=(oK - mean(sK)) / sd(sK)
o_normK
o_normW=(oW - mean(sW)) / sd(sW)
o_normW
o_normD=(oD - mean(sD)) / sd(sD)
o_normD

################################################################################
###################################4############################################

#function to calculate eucledean distance
calculateEucledean<- function(oK,oW,oD,sK,sW,sD) {
  dist=sqrt((oK - sK)^2+(oW - sW)^2+(oD - sD)^2)
  return(dist)
}
#get parameter values
parameter_values <- readLines("pars_final.txt")

distances=vector(length = 10000)
for(i in 1:length(simulated_datasets)){
  distances[i]=calculateEucledean(o_normK,o_normW,o_normD,normK[i],normW[i],normD[i])
}
distances
################################################################################
##################################5#############################################


#find 500 smallest distances and their indexes
indexes= order(distances)[1:500]
distances=distances[indexes]

################################################################################
##################################6#############################################

#get values from pars_final.txt (as strings)
parameter_values=parameter_values[indexes]
parameter_values
################################################################################
##################################7#############################################

#mean and median of the parameter values (as.numeric to convert from string)
mean_=mean(as.numeric(parameter_values))
mean_
median_=median(as.numeric(parameter_values))
median_
################################################################################
##################################8#############################################

#histogram
hist(as.numeric(parameter_values), main = "Histogram of Parameter Values", xlab = "Parameter Values", col= "chocolate4")
#density plot
density_plot=density(as.numeric(parameter_values))
plot(density_plot, main = "Density Plot of Parameter Values", xlab = "Parameter Values", col = "brown1")

################################################################################

