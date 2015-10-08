Permutation test with the t-statistic:

This script is usefull if the user has to compare quantitative measurements from matrices that do not follow a normal distribution. It is not a simple Non-parametric test because do not consider outliers. The procedure is similar of a paired-test but with 15000 permutations. 

Assumption: The dataset do not follow a "T-Normal" distribution

Input Two Matrices: Control and Treatment

Return: a T observed and a T permuted

It is possible to select the treshold by which the T-values observed are significant or not to the respect of the T-values permuted. 

In Script.r it is possible to find an example on how to use the script.  
