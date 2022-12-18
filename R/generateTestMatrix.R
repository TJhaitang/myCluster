generateTestMatrix = function(n=10)
{
    x=runif(n)
    y=runif(n)
    s=cbind(x,y)
    out=dist(s)
    return(out)
}
