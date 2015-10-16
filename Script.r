#Permutation Script on paired data
#Given two matrix, one control and the other treatment that do not follow normal T distribution
#Return a T observed and a T empirical

#Import the data
psi <- read.table("psi_values_cassette.filtered.50%_population.tsv",header=TRUE)

#Eliminate NA
psi <- na.omit(psi)

#Split the matrixes in 2 matrix each one for condition
psi_infected <- psi[,seq(from=1,to=79, by=2)]
psi_naive <- psi[,seq(from=2, to=78, by=2)]
#Eliminate the column of the exons
psi_infected <- psi_infected[,-1]

#Look the distribution of the mean
psi_infected_mean <- colMeans(psi_infected)
psi_naive_mean <- colMeans(psi_naive)
diff_mean <- psi_infected_mean-psi_naive_mean
hist(diff_mean, main="Mean Differences Distribution",xlab="Mean Difference Values",ylab="Samples")

#Compute the T observed
#Calculate the difference for each value of the two matrix
Difference_Matrix <- matrix(, nrow = 126, ncol = 39)
for (i in 1:126){
  for(j in 1:39)
    Difference_Matrix[i,j] = psi_infected[i,j]-psi_naive[i,j]
}

#Calculate the mean value for each exon and Standard Deviation
Mean <- apply(Difference_Matrix, 1, mean)
SD <- apply(Difference_Matrix,1, sd)

#Calculate T observed 
T_obs <- Mean/(SD/sqrt(x = 38))
 

#Compute Differences of conditions in order to perform permutation (same as before) 
psi_diff <- matrix(, nrow = 126, ncol = 39)
for (i in 1:126){
  for(j in 1:39)
    psi_diff[i,j] = psi_infected[i,j]-psi_naive[i,j]
}

###################################################
##Compute permuted T value in one random matrix####
###################################################

#First create a matrix in which the pairs are swapped
mask <- t(replicate(nrow(psi_diff),sample(c(TRUE,FALSE),ncol(psi_diff),replace=T)))
psi_diff[mask] <- -psi_diff[mask]

#Then Calculate the mean value for each exon and Standard Deviation
Mean <- apply(psi_diff, 1, mean)
SD <- apply(psi_diff,1, sd)

#Finally Calculate T empirical 
T_emp <- Mean/(SD/sqrt(x = 39))
 
#####################################
#Compute 15000 permuted T value######
#####################################
n_sampling <- 15000
n_rows <- 126
t_value_emp_matrix = matrix(, nrow = n_rows, ncol = n_sampling)
for (i in 1:n_sampling){
  print(i)
  mask <- t(replicate(nrow(psi_diff),sample(c(TRUE,FALSE),ncol(psi_diff),replace=T)))
  psi_diff[mask] <- -psi_diff[mask]  
  Mean_t_emp_all <- apply(psi_diff,1, mean)
  SD_t_emp_all <- apply(psi_diff,1, sd)
  t_value_emp_matrix[,i] <- Mean_t_emp_all/(SD_t_emp_all/sqrt(x = 39))
}

#######################################################
##Now Compare the permuted with the observed t-value##
#######################################################

t_emp_cols <- length(t_value_emp_matrix[1,])
T_emp <- vector(mode="numeric", length=n_rows)
for (i in 1:n_rows) {
  T_emp[i] <- 0
}

for (i in 1:length(T_obs)){
  for (j in 1:t_emp_cols) {
    if ( !is.na(t_value_emp_matrix[i,j])) {
      if (abs(T_obs[i]) <= abs(t_value_emp_matrix[i,j])){
        T_emp[i] <- T_emp[i] + 1
      }
    }
  }
  T_emp[i] <- T_emp[i] / 15000
} 

#Check differential expressed values according to 10^-1 and different than  remind that the number of exons and 
#of Tvalues have to be the same, if not this means that something is wrong
exons <- psi$exons
tvalues <- as.numeric(as.vector(T_emp))
df <- cbind(exons,tvalues)
df_subset <- subset(df,df$tvalues<=0.1 & df$tvalues!= 0) #this is your table with name of exons and tvalue significant

