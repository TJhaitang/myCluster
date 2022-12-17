hclust <- function(d,method="single",members=NULL){
    METHODS=c("single", "complete", "average", "mcquitty", "ward.D", "centroid", "median", "ward.D2")
    method=pmatch(method,METHODS)
}