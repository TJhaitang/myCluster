generateTestMatrix = function()
{
    x=runif(10)
    y=runif(10)
    s=cbind(x,y)
    out=dist(s)
    return(out)
}
