# 基于prim算法的最短距离法层次聚类算法

## 函数说明

### functions: 
* myCluster::hclust(d,method="single")

### input:
* d: 一个由dist()生成的距离矩阵
* method: 聚类方法，可选值为"single",其他方法敬请期待

### output:
* result: 一个list，包含stats::hclust()函数的返回结果
* result$merge: 一个n-1行2列的矩阵，每一行表示合并的两个簇的编号，若编号k小于0，则-k为样本点的编号；若编号k大于0，则其为第k行合并的簇的编号
* result$height: 一个n-1行1列的矩阵，每一行表示合并的两个簇的距离
* result$order: 一个n行1列的矩阵，每一行表示样本点的编号，按此顺序排列可使得聚类树不交叉
* result$labels: 一个n行1列的矩阵，每一行表示样本点的标签
* result$call: 一个list，包含调用hclust()函数时的参数
* result$method: 一个字符向量，表示聚类方法
* result$dist.method: 一个字符向量，表示距离计算方法

## 算法描述

### 1.最短距离法层次聚类算法
对于给定的n个样本点，首先将每个样本点看成一个簇，然后计算任意两个簇之间的距离，将距离最近的两个簇合并，形成一个新的簇，然后再计算新簇与其他簇之间的距离，将距离最近的两个簇合并，形成一个新的簇，如此重复，直到所有的样本点都被合并为一个簇为止。   
在实际应用时，该算法需要维护一个n×n的距离矩阵，每次聚类都会遍历并更新此矩阵。因此，该算法的时间复杂度较高，为O(n^3)。

### 2.kruskal算法
Kruskal算法是典型的最小生成树算法，它的基本思想是：每次从所有的边中选择一条最小的边，如果这条边连接的两个顶点不在同一个连通分量中，则将这条边加入到最小生成树中，否则舍弃这条边，直到所有的顶点都在同一个连通分量中为止。   
可以发现krustal算法的生成最小生成树的过程即为最短距离法层次聚类算法的过程，然而krustal算法的时间复杂度为O(eloge)，在聚类问题中距离矩阵代表一个完全图，因此krustal算法的时间复杂度为O(n^2logn)。

### 3.prim算法
Prim算法是典型的最小生成树算法，它的基本思想是：每次从所有的顶点中选择一个最小的顶点，如果这个顶点与已经在最小生成树中的顶点相连，则将这条边加入到最小生成树中，否则舍弃这条边，直到所有的顶点都在同一个连通分量中为止。   
由于最短距离法层次聚类算法的过程即为krustal算法的过程，因此同样作为最小生成树算法的prim算法也可以进行最短距离法层次聚类算法的实现，而prim算法的时间复杂度为O(n^2)，远远优于上述两个方法，所以可以使用此方法进行最短距离法层次聚类算法的实现。

### 4.深度优先搜索算法
深度优先搜索算法是一种用于遍历或搜索树或图的算法。沿着树的深度遍历树的节点，尽可能深的搜索树的分支。当节点v的所在边都己被探寻过，搜索将回溯到发现节点v的那条边的起始节点。这一过程一直进行到已发现从源节点可达的所有节点为止。如果还存在未被发现的节点，则选择其中一个作为源节点并重复以上过程，整个进程反复进行直到所有节点都被访问为止。   
由于stats::hcluster()函数中的返回值包含order信息，此信息使得将元素按order排列画出的聚类图像将不会出现交叉，因此可以使用深度优先搜索算法对生成树进行遍历，从而得到order信息。

## 代码实现

### 1. prime算法
```{cpp}
clock_t clusterMethod1(int n,arma::mat * d,clusterChain *chain,int method=0){
    clock_t test = clock();
    clock_t sum1 = 0;
    dubLinkList *list = new dubLinkList(n);
    //建立一个数组存储每个点的最短距离
    double *tempDist = new double[n];
    //第一次迭代
    int i1 = 0;
    int i2 = 1;
    double minDist = std::numeric_limits<double>::infinity();
    for(int i=1;i<n;i++) {
        tempDist[i] = (*d)(i1,i);
        if(tempDist[i] < minDist) {
            minDist = tempDist[i];
            i2 = i;
        }
    }
    chain->addNode(i1,i2,minDist);
    list->remove(i2);
    // clock_t median = clock();
    // //开始迭代
    for(int i=2;i<n;i++) {//Prim算法
        i1 = i2;
        i2 = list->succ(0);
        minDist = tempDist[i2];
        for(int j=i2;j!=-1;j=list->succ(j)) {
            if((*d)(i1,j) < tempDist[j]) {
                tempDist[j] = (*d)(i1,j);
            }
            if(tempDist[j] < minDist) {
                minDist = tempDist[j];
                i2 = j;
            }
        }
        chain->addNode(i1,i2,minDist);
        list->remove(i2);
    }
    // delete[] tempDist;
    // delete list;
    chain->isOrdered = false;
    sum1+=clock()-test;
    return test;
}
```

### 2. 由边集合生成聚类结果
由于上述算法仅保存了最小生成树的边集合，因此需要将边集合转化为聚类结果，即将边集合转化为聚类树   
由于将边按照权重从小到大排序后即聚类顺序，此时若一条边连接了一个已经出现过的点和一个未出现过的点，这意味着这一步是将这未出现过的点与已出现过的点所在的类进行合并。从而将最小生成树转化为聚类二叉树。
```{cpp}
class parentFinder{//一颗省内存和时间的树
private:
    int * parent;
    int n;
public:
    parentFinder(int n):n(n) {
        parent = new int[2*n-1];
        for(int i=0;i<2*n-1;i++) {
            parent[i] = 0;
        }
    }
    ~parentFinder() {
        delete[] parent;
    }

    int find(int i) {
        if(parent[i] == 0) {
            return i;
        } else {
            int idx= find(parent[i]);
            parent[i] = idx;
            return idx;
        }
    }

    void merge(int i1,int i2) {
        parent[i1]=parent[i2]=n++;
    }
};
parentFinder *pf = new parentFinder(n);
for(int i=0;i<n-1;i++) {
    clusterNode *node = (*chain)[i];
    int i1 = pf->find(node->idx1);//将元素间的合并转化为类之间的合并
    int i2 = pf->find(node->idx2);
    if(i1 > i2) {//画图好看一点
        int temp = i1;
        i1 = i2;
        i2 = temp;
    }
    result(i,0) = i1<n?(-i1-1):(i1-n+1);
    result(i,1) = i2<n?(-i2-1):(i2-n+1);
    result(i,2) = node->dist;
    pf->merge(i1,i2);
    // nodeSize[i]=(i1<n?1:nodeSize[i-n])+(i2<n?1:nodeSize[i-n]);
}
```

### 3. 由聚类二叉树生成order信息
```{cpp}
//DFS
int pos=0;
std::stack<int> s;
s.push(n-1);//此节点的编号为n-1，位置为n-2
while(!s.empty()) {
    int node = s.top();
    s.pop();
    if(node < 0) {
        result(pos,3) = -node;
        pos++;
    } else {
        s.push(result(node-1,1));
        s.push(result(node-1,0));
    }
}
```

### 4. R结果处理
```{r}
output=myCluster(n,d,method)
#output:
#前两列为merge
#第三列为height
#第四列为order

#将输出转换成hclust的输出
result=list()
#取output的前两列n-1行作为result$merge
result$merge=output[1:(n-1),1:2]
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

return(result)
```

## 速度测试
```{r}
library(myCluster)
library(microbenchmark)

# Generate a random distance matrix
mat=generateTestMatrix(5000)

# microbenchmark
microbenchmark(myCluster::hclust(mat,"single"),stats::hclust(mat,"single"),times=5)
```
```{r}
Unit: milliseconds
                             expr      min       lq     mean   median       uq      max neval
 myCluster::hclust(mat, "single") 306.5232 306.8167 308.4409 306.8522 309.3882 312.6243     5
     stats::hclust(mat, "single") 555.2470 555.4061 568.9930 561.4374 569.1808 603.6938     5
```

## 两次速度提升
### 1.不在R中进行数据格式转换(-75%)
```{R}
dd=as.matrix(d)
```

### 2.将传递实例改为传递指针(-40%)
```{cpp}
clock_t clusterMethod1(int n,arma::mat * d,clusterChain *chain,int method=0)
```

## 示例
```{R}
library(myCluster)
library(graphics)

x=runif(10)
y=runif(10)
d=dist(rbind(x,y))
result=hclust(d,"single")

plot(result)
```

## 优化空间
### 1.使用更快的深度优先搜索算法
由于生成树边保证了左结点的值小于右结点的值，因此若右节点为叶子节点，则左节点也为叶子节点；若左结点不为叶子节点，则右节点也不为叶子节点。当前搜索并未考虑这一信息

### 2.不再进行数据格式转换
在当前实现中我将dist数据转换为了arma::mat，而dist数据其实是以arma::vec存储的下三角矩阵，因此可以直接使用dist数据，不再进行转换

## 项目github地址
https://github.com/TJhaitang/myCluster