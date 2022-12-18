hclust <- function(d,method="single",members=NULL){
    METHODS=c("single", "complete", "average", "mcquitty", "ward.D", "centroid", "median", "ward.D2")
    method=pmatch(method,METHODS)
    if(method==-1){
        stop("method must be one of ",paste(METHODS,collapse=", "))
    }
    #输入一个对称矩阵
    tryCatch({
        d=as.matrix(d)
    },error=function(e){
        stop("d must be a distance matrix")
    })
    if(!is.matrix(d)){
        stop("d must be a distance matrix")
    }

    output=myCluster(d,method)
    #output:
    #output$merge1
    #output$merge2
    #output$height
    #output$order--?

    #这里需要再处理一下输出
    return(output)
}