// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <algorithm> // for std::sort
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

class clusterNode {
public:
    int idx1;
    int idx2;
    double dist;
    clusterNode(int i1,int i2,double d):idx1(i1),idx2(i2),dist(d) {}
    bool operator<(const clusterNode &node) const {//用来排序
        return this.dist < node.dist;
    }
};

class clusterChain {
public:
    clusterNode *nodes;
    int n;
    int nNodes;
    bool isOrdered;
    clusterChain(int n):n(0),nNodes(n),isOrdered(false) {
        nodes = new clusterNode[nNodes];
    }
    ~clusterChain() {
        delete[] nodes;
    }
    void addNode(int i1,int i2,double d) {
        nodes[n].idx1 = i1;
        nodes[n].idx2 = i2;
        nodes[n].dist = d;
        n++;
    }
    clusterNode * operator[](int i) {//这样写可能有问题，不太记得了
        return &nodes[i];
    }

    void sort(){
        std::sort(nodes,nodes+nNodes-1);
        // isOrdered = true;//这里还得写到其他方法再说
    }
};

class dubLinkList {//维护活跃结点序号
private:
    int * succ;
    int * pred;
public:
    dubLinkList(int n) {
        succ = new int[n];
        pred = new int[n];
        for(int i=0;i<n;i++) {
            succ[i] = i+1;
            pred[i] = i-1;
        }
        succ[n-1] = -1;
        pred[0] = -1;
    }
    ~dubLinkList() {
        delete[] idx;
    }
    int succ(int i) {
        return succ[i];
    }
    int pred(int i) {
        return pred[i];
    }
    void remove(int i) {
        if(pred[i] != -1) {
            succ[pred[i]] = succ[i];
        }
        if(succ[i] != -1) {
            pred[succ[i]] = pred[i];
        }
    }
};

class parentFinder{
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
}

// 计算两个向量的距离

// 使用最短距离法进行聚类
//     The basis of this algorithm is an algorithm by Rohlf:
//     F. James Rohlf, Hierarchical clustering using the minimum spanning tree,
//     The Computer Journal, vol. 16, 1973, p. 93–95.
void clusterMethod1(int n,arma::mat d,clusterChain *chain,int method=0,arma::vec members=NULL){
    //建立双向链表存储未聚类的点
    dubLinkList *list = new dubLinkList(n);
    //建立一个数组存储每个点的最短距离
    double *tempDist = new double[n];
    //第一次迭代
    int i1 = 0;
    int i2 = 1;
    double minDist = std::numeric_limits<double>::infinity();
    for(int i=1;i<n;i++) {
        tempDist[i] = d(i1,i);
        if(tempDist[i] < minDist) {
            minDist = tempDist[i];
            i2 = i;
        }
    }
    chain->addNode(i1,i2,minDist);
    list->remove(i2);
    //开始迭代
    for(int i=2;i<n;i++) {
        i1 = i2;
        i2 = list->succ(0);
        minDist = tempDist[i2];
        for(int j=i2;j!=-1;j=list->succ(j)) {
            if(d(i1,j) < tempDist[j]) {
                tempDist[j] = d(i1,j);
            }
            if(tempDist[j] < minDist) {
                minDist = tempDist[j];
                i2 = j;
            }
        }
        chain->addNode(i1,i2,minDist);
        list->remove(i2);
    }
    delete[] tempDist;
    delete list;
    chain->isOrdered = false;
}

//根据合并顺序生成聚类结果
void generateResult(int n,clusterChain *chain,arma::mat result) {
    if(!chain->isOrdered) {
        chain->sort();
    }
    parentFinder *pf = new parentFinder(n);
    for(int i=0;i<n-1;i++) {
        clusterNode *node = (*chain)[i];
        int i1 = pf->find(node->idx1);
        int i2 = pf->find(node->idx2);
        if(i1 > i2) {
            int temp = i1;
            i1 = i2;
            i2 = temp;
        }
        result(i,0) = i1<n?(-i1-1):(i1-n+1);
        result(i,1) = i2<n?(-i2-1):(i2-n+1);
        result(i,2) = node->dist;
        pf->merge(i1,i2);
    }
}


// 返回一个包含merge,height的mat
// [[Rcpp::export]]
arma::mat myCluster(arma::mat d,int method=0,arma:vec members=NULL) {
    // 初始化
    int n = d.n_rows;
    clusterChain *chain = new clusterChain(n-1);
    switch(method) {
        case 0:
            clusterMethod1(n,d,chain,method,members);
            break;
        default:
            break;
    }
    // 生成结果
    arma::mat result(n-1,3);
    generateResult(n,chain,result);
    delete chain;
    return result;
}