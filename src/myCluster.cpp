// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <algorithm> // for std::sort
#include <cstddef>
#include <limits> // for std::numeric_limits
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

class clusterNode {
public:
    int idx1;
    int idx2;
    double dist;
    clusterNode(int i1,int i2,double d):idx1(i1),idx2(i2),dist(d) {}
    clusterNode(){idx1=0;idx2=0;dist=0;}
    bool operator<(const clusterNode &node) const {//用来排序
        return this->dist < node.dist;
    }
};

class clusterChain {
public:
    clusterNode *nodes;
    int n;
    int nNodes=0;
    bool isOrdered;
    // int test;
    clusterChain(int n):n(0),nNodes(n),isOrdered(false) {
        nodes = new clusterNode[nNodes];
        // test=0;
    }
    ~clusterChain() {
        delete[] nodes;
    }
    void addNode(int i1,int i2,double d) {
        nodes[n].idx1 = i1;
        nodes[n].idx2 = i2;
        nodes[n].dist = d;
        // test=12138;
        n++;
    }
    clusterNode * operator[](int i) {//这样写可能有问题，不太记得了
        return &nodes[i];
    }

    void sort(){
        std::sort(nodes,nodes+nNodes);//last并不是指最后一个元素而是其再之后的地址？
        // isOrdered = true;//这里还得写到其他方法再说
    }

    void order(int * nodeSize,int num){//这个方法是为了将聚类结果按照层次结构排列，具体方法可以再理解一下
        struct pos_node{
            int pos;
            int node;
        };
        pos_node *pos_nodes = new pos_node[num];
        pos_nodes[0].pos = 0;
        pos_nodes[0].node = num-2;
        int idx=1,parent=0,child=0,pos=0;
        int orderList[num]={0};
        do{
            idx-=1;
            parent = pos_nodes[idx].node;
            pos = pos_nodes[idx].pos;
            child = nodes[parent].idx1;
            // if(nodeSize[child] > 1) {
            //     pos_nodes[idx].pos = pos+1;
            //     pos_nodes[idx].node = child;
            //     idx+=1;
            // } else {
            //     nodes[child].idx1 = pos;
            // }
            if(child<0){
                orderList[pos] = -child;
                pos+=1;
            }
            else{
                pos_nodes[idx].pos = pos;
                pos_nodes[idx].node = child-1;
                idx+=1;
                pos+=nodeSize[child-1];
            }
            child = nodes[parent].idx2;
            if(child<0){
                orderList[pos] = -child;
            }
            else{
                pos_nodes[idx].pos = pos;
                pos_nodes[idx].node = child-1;
                idx+=1;
            }
        }while(idx>0);
        pos=0;
        int mask[num]={0};
        for(int i=0;i<num;i++){
            if(mask[orderList[i]]==0)
            for(int j=0;j<num;j++){
                if(orderList[i]==nodes[j].idx1){//可能需要复制构造函数
                    mask[orderList[i]]=1;
                    if(nodes[j].idx2<0){
                        mask[nodes[j].idx2]=1;
                    }
                    clusterNode temp = nodes[j];
                    nodes[j] = nodes[pos];
                    nodes[pos] = temp;
                    pos+=1;
                }
            }
        }
    }
};

class dubLinkList {//维护活跃结点序号
private:
    int * succList;
    int * predList;
public:
    dubLinkList(int n) {
        succList = new int[n];
        predList = new int[n];
        for(int i=0;i<n;i++) {
            succList[i] = i+1;
            predList[i] = i-1;
        }
        succList[n-1] = -1;
        predList[0] = -1;
    }
    ~dubLinkList() {
        delete[] succList;
        delete[] predList;
    }
    int succ(int i) {
        return succList[i];
    }
    int pred(int i) {
        return predList[i];
    }
    void remove(int i) {
        if(predList[i] != -1) {
            succList[predList[i]] = succList[i];
        }
        if(succList[i] != -1) {
            predList[succList[i]] = predList[i];
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
};

// 计算两个向量的距离

// 使用最短距离法进行聚类
//     算法来源:
//     F. James Rohlf, Hierarchical clustering using the minimum spanning tree,
//     The Computer Journal, vol. 16, 1973, p. 93–95.
void clusterMethod1(int n,arma::mat d,clusterChain *chain,int method=0){
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
    // //开始迭代
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
    // delete[] tempDist;
    // delete list;
    chain->isOrdered = false;
}

//根据合并顺序生成聚类结果
arma::mat generateResult(int n,clusterChain *chain) {
    arma::mat result(n,4);
    if(!chain->isOrdered) {//其他方法好像不太一样，暂时先这么写
        chain->sort();
    }
    parentFinder *pf = new parentFinder(n);
    int * nodeSize = new int[n-1];
    for(int i=0;i<n-1;i++) {
        clusterNode *node = (*chain)[i];
        int i1 = pf->find(node->idx1);//将元素间的合并转化为类之间的合并
        int i2 = pf->find(node->idx2);
        if(i1 > i2) {//画图好看一点
            int temp = i1;
            i1 = i2;
            i2 = temp;
        }
        node->idx1 = i1<n?(-i1-1):(i1-n+1);
        node->idx2 = i2<n?(-i2-1):(i2-n+1);
        pf->merge(i1,i2);
        nodeSize[i]=(i1<n?1:nodeSize[i-n])+(i2<n?1:nodeSize[i-n]);
    }
    // chain->order(nodeSize,n);
    for(int i=0;i<n-1;i++) {
        clusterNode *node = (*chain)[i];
        result(i,0) = node->idx1;
        result(i,1) = node->idx2;
        result(i,2) = node->dist;
    }
    struct pos_node{
        int pos;
        int node;
    };
    pos_node *pos_nodes = new pos_node[n];
    pos_nodes[0].pos = 0;
    pos_nodes[0].node = n-2;
    int idx=1,parent=0,child=0,pos=0;
    do{
        idx-=1;
        parent = pos_nodes[idx].node;
        pos = pos_nodes[idx].pos;
        child = chain->nodes[parent].idx1;
        // if(nodeSize[child] > 1) {
        //     pos_nodes[idx].pos = pos+1;
        //     pos_nodes[idx].node = child;
        //     idx+=1;
        // } else {
        //     nodes[child].idx1 = pos;
        // }
        if(child<0){
            result(pos,3) = -child;
            pos+=1;
        }
        else{
            pos_nodes[idx].pos = pos;
            pos_nodes[idx].node = child-1;
            idx+=1;
            pos+=nodeSize[child-1];
        }
        child = chain->nodes[parent].idx2;
        if(child<0){
            result(pos,3) = -child;
        }
        else{
            pos_nodes[idx].pos = pos;
            pos_nodes[idx]. node = child-1;
            idx+=1;
        }
    }while(idx>0);
    delete chain;
    return result;
}


// 返回一个包含merge,height的mat
// [[Rcpp::export]]
arma::mat myCluster(arma::mat d,int method=0) {//members应该设置默认为空
    // 初始化
    int n = d.n_rows;
    clusterChain *chain = new clusterChain(n-1);
    switch(method) {
        case 1://这里de了2个小时bug,R下标从1开始
            clusterMethod1(n,d,chain,method);
            break;
        default:
            break;
    }
    // 生成结果
    return generateResult(n,chain);

} 