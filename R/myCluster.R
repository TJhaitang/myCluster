hclust <- function(d,method="single",members=NULL){
    METHODS=c("single", "complete", "average", "mcquitty", "ward.D", "centroid", "median", "ward.D2")
    method=pmatch(method,METHODS)
    if(method==-1){
        stop("method must be one of ",paste(METHODS,collapse=", "))
    }
    #输入一个对称矩阵
    # s=proc.time()
    # tryCatch({
    #     dd=as.matrix(d)#若不是回头狠狠测试，谁会想到这一句占据了75%的运行时间呢
    # },error=function(e){
    #     stop("d must be a distance matrix")
    # })
    # print(proc.time()-s)
    # if(!is.matrix(dd)){
    #     stop("d must be a distance matrix")
    # }
    # s=proc.time()
    n=attr(d,"Size")
    output=myCluster(n,d,method)
    #output:
    #前两列为merge
    #第三列为height
    #第四列为order

    #将输出转换成hclust的输出
    result=list()
    #取output的前两列n-1行作为result$merge
    result$merge=output[1:(n-1),1:2]#可以这么写吧
    #取output的第三列n-1行作为result$height
    result$height=output[1:(n-1),3]
    #取output的第四列n行作为result$order
    result$order=output[,4]

    result=c(result,
        list(
            labels=attr(d,"labels"),
            method=METHODS[method],
            call=match.call(),
            dist.method=attr(d,"method")
        )
    )
    class(result)="hclust"

    # print(proc.time()-s)

    return(result)
}