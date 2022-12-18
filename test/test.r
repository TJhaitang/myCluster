library(myCluster)
library(microbenchmark)

# Generate a random distance matrix
mat=generateTestMatrix(5000)

# microbenchmark
microbenchmark(myCluster::hclust(mat,"single"),stats::hclust(mat,"single"),times=5)

output=myCluster::hclust(mat,"single")
output[5000,]