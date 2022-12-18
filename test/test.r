library(myCluster)
library(microbenchmark)

# Generate a random distance matrix
mat=generateTestMatrix(5000)

# microbenchmark
microbenchmark(myCluster::hclust(mat,"single"),stats::hclust(mat,"single"),times=5)

# Unit: milliseconds
#                              expr      min       lq     mean   median       uq      max neval
#  myCluster::hclust(mat, "single") 306.5232 306.8167 308.4409 306.8522 309.3882 312.6243     5
#      stats::hclust(mat, "single") 555.2470 555.4061 568.9930 561.4374 569.1808 603.6938     5




output=myCluster::hclust(mat,"single")
output[5000,]


library(myCluster)
mat=generateTestMatrix()
out=myCluster::hclust(mat,"single")