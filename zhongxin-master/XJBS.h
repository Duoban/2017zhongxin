#ifndef __ROUTE_H__
#define __ROUTE_H__

#include <deque>
#include <queue>
#include <stack>
#include <vector>
#include <algorithm>
#include <cstring>
#include <utility>
#include <cstdio>
#include <climits>
#include <iostream>
#include <cstdlib>
#include <random>
#include <cmath>
#include <ctime>
#include <cmath>
#include "XJBS.h"
#include "lib_io.h"
#include "lib_time.h"

using namespace std;

const int INF = 0x33333333;
const int vbound = 2050;

struct n_nEdge {
    int v;
    int c;
    double pro; //每条边被访问的概率
    n_nEdge(int _v, int _c): v(_v), c(_c) {}
};

struct limit_Edge {
    int t, u, c;           
    limit_Edge(){}
    limit_Edge(int _t,int _u, int _c): t(_t), u(_u), c(_c) {}
 };

struct reward_Edge {
    int t, u, c; 
    reward_Edge(){}
    reward_Edge(int _t,int _u, int _c): t(_t), u(_u), c(_c) {}
 };

class Graph {
public:
    void init_graph(char * topo[MAX_EDGE_NUM], int line_num);   //读取信息
    
    void SPFA_weight();                         //SPFA算法计算带权最短路径，强制约束必过奖励边，存放route
    void SPFA_For_SG();                         //SPFA算法计算简化图的带权最短路径，约束奖励边，和最大跳数，不保存route

    int get_cost(int v, int id);                //获取v到id节点的开销
    void set_cost(int v, int id, int cost);     //设置v到id节点的开销 
    void Show_best_route(int s,int e);          //显示s到e在SPFA_weight之后产生在route_weight中的最短路径

    vector<vector<n_nEdge> > graph;             //图信息
    vector<limit_Edge> limits;                  //害虫的位置
    vector<reward_Edge> rewards;                //奖励边的位置

    vector<vector<int> > dist_weight;           //带权最短路径权值
    vector<vector<int> > dist_w_j;              //带权最短路径对应跳数 
    vector<vector<int> > route_weight;          //带权路径

    vector<vector<int> > SG_dist_w;             //简化图的最短路径权值
    vector<vector<int> > SG_dist_j;             //对应的跳数

    vector<int > must_node;                     //放必须经过的节点

    int vis[vbound];                            //一个bool数组用来一般做节点是不是经过的处理
   
    int limit_depth;                            //最大跳数
    int node_num, edge_num;                     //图节点数 ， 边数 
    int must_node_num;                          //必须经过节点的个数
    int reward_edge_num , forbidden_edge_num;   //奖励边数，禁止边数 

    int start,end;                              //起点终点
};

class DFS_Search {
    public:
        DFS_Search(Graph & graph);

        void DFS_withcutoff(int start);     //DFS+剪枝查找

        int sum_stack(stack<int> in);       //求路径栈权值总和
        bool isok(stack<int> in);           //判断路径是否满足条件

        Graph *g;                           //导入的原图
        stack<int> route;                   //存放当前路径
        stack<int> best_route;              //存放最优路径
        int best_cost;                      //最优代价
        int best_jump;                      //最优跳数
};

class Simplified_Graph_Search{
    public:
        Simplified_Graph_Search(Graph & graph);
        void Simplified_Graph();

        Graph *g;                       //导入原图
        Graph simplified_g;             //构造出的简化图（完全图）
        vector<int > g_node;            //所有的约束点集合

        /**********DFS查找************/
        void DFS_search(int start);     //DFS遍历+剪枝方法遍历查找缩略图

        int sum_stack(stack<int> in);   //求路径栈的权值和
        int sum_jump(stack<int> in);    //求路径栈在原图上的跳数总和
        bool isok(stack<int> in);       //判断当前路径栈路径是否满足约束条件
        void Show_Route();              //显示缩略图在原图上的的路径
        /*--------------------------*/
        stack<int> route;               //存放当前路径
        stack<int> best_route;          //存放最优路径
        int best_cost;                  //最优路径开销
        int best_jump;                  //最优路径跳数
        /*****************************/

        /***************TSP***********/
        void TSP_taboo_search();        //禁忌算法
        
        int MAX_GEN;              //迭代次数
        int N;                          //每次搜索邻居个数
        int ll;                         //禁忌长度
        int cityNum;                    //城市数量，编码长度
        
        vector<int > Ghh;               //初始路径编码
        
        int bestT;                      //最佳出现代数
        vector<int > bestGh;            //最好的路径编码
        int bestEvaluation;             //最佳权值
        int bestJump;                   //最佳跳数

        vector<int > LocalGhh;          //当前最好编码
        int localEvaluation;            //当前权值
        int localJump;                  //当前跳数
        vector<int > tempGhh;           //存放临时编码
        int tempEvaluation;             //临时权值
        int tempJump;                   //临时跳 

        vector<vector<int >> jinji;     //禁忌表
        int t;                          //当前代数

        void init_taboo(int g,int c,int m);              //初始化taboo算法所需变量
        void initGroup();               //初始化编码Ghh
        void copyGh(vector<int> &Gha,vector<int> &Ghb); //复制编码体，Gha->Ghb
        int evaluate(vector<int> &chr);         //计算权值
        int evaluate_J(vector<int> &chr);       //计算跳数
        void Linju(vector<int> &Gh,vector<int> &tempGh);//领域交换，得到邻居
        int panduan(vector<int> &tempGh);       //判断编码是否在禁忌表中
        void jiejinji(vector<int> &tempGh);     //解禁忌或者加入禁忌
        int isokGhh(vector<int> &tempGh);       //判断此时编码是否满足经过奖励边
        void print();                           //打印最优解
        /*****************************/
};


void XJBS(char * topo[MAX_EDGE_NUM], int line_num,char * filename);
#endif
