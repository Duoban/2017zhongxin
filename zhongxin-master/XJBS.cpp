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
#include "lib_io.h"
#include "XJBS.h"
#include "antsimulate.h"
#include <sys/time.h>

using namespace std;

void XJBS(char * topo[MAX_EDGE_NUM], int line_num,char * filename)
{
    //char * topo_file;
    
    Graph G;                                  
    G.init_graph(topo, line_num);              //读数据
    double t1 = (double)clock() / CLOCKS_PER_SEC;
    printf("\n初始化图用时：%f \n",t1);

    /***********************************************/
    G.SPFA_weight();                                ////利用SPFA将所有点的最短路径找出来/这个里面的队列有优化的方式
    double t2 = (double)clock() / CLOCKS_PER_SEC; 
    printf("任意节点之间带权最短路径（费用）SPFA算法用时：%f \n",t2-t1);
    /***********************************************/

    /***********************************************/
    G.SPFA_For_SG();
    double t4 = (double)clock() / CLOCKS_PER_SEC;
    printf("任意节点带权最短路径并且满足最大跳数约束的SPFA算法（改）用时：%f \n",t4-t2);
    /***********************************************/

    /***********************************************/
    DFS_Search d_search = DFS_Search(G);
   // d_search.DFS_withcutoff(G.start);
    printf("\n利用DFS+裁剪算法,求解多约束条件下的带权最短路径路径：\n");
    int f = 0;
    if(d_search.best_route.empty()){
        f = 1;
        goto no_answer;
    }    
    printf("DFS算法找的最短开销：%d\n", d_search.best_cost);
    printf("DFS算法找的最短跳数；%d\n", d_search.best_jump);
    printf("DFS算法找到的路径为：\n");
    while(!d_search.best_route.empty()){
        printf("%d ",d_search.best_route.top());
        d_search.best_route.pop();
    }
    printf("\n");
no_answer:
    if (f)
        printf("当前图无解\n");
    double t5 = (double)clock() / CLOCKS_PER_SEC;
    printf("DFS+剪枝算法全局搜索最优解用时：%f \n",t5-t4);
    /***********************************************/

    /***********************************************/
    printf("\n利用化简图算法思路简化问题复杂度："); 
    Simplified_Graph_Search sg_search = Simplified_Graph_Search(G);
    sg_search.Simplified_Graph();
    double t6 = (double)clock() / CLOCKS_PER_SEC;
    printf("\n化简图形成完全图用时：%f \n",t6-t5);

    printf("\n利用DFS+裁剪算法,在简化的完全图上，求解多约束条件下的带权最短路径路径：");
    sg_search.DFS_search(G.start);
    sg_search.Show_Route();
    double t7 = (double)clock() / CLOCKS_PER_SEC;
    printf("简化图DFS全局搜索最优解用时：%f \n",t7-t6);
    
    printf("\n利用禁忌算法,在简化的完全图上，将约束条件下的带权最短路径路径转化为TSP问题：");
    printf("\n");
    sg_search.init_taboo(1000,200,20);
    sg_search.TSP_taboo_search();
    sg_search.print();
    printf("\n");
    double t8 = (double)clock() / CLOCKS_PER_SEC;
    printf("简化图TSP搜索最优解用时：%f \n",t8-t7);
    /***********************************************/

    printf("\n利用启发式蚁群算法，对原图直接进行最优解的搜索：\n");
    ant_simulate(G); //蚁群算法
    double t9 = (double)clock() / CLOCKS_PER_SEC;
    printf("蚁群算法用时: %f \n",t9-t8);
    
    //write_result(topo_file, filename);
    //delete []topo_file;
}

void Graph::init_graph(char * topo[MAX_EDGE_NUM], int line_num) { 
    int line = 0;
    int u, v, c;
    if (line < line_num)
        sscanf(topo[line], "%d %d %d %d %d %d", &node_num, &edge_num, &reward_edge_num, &forbidden_edge_num, &limit_depth, &must_node_num);

    graph.resize(node_num, vector<n_nEdge>());
    limits.resize(forbidden_edge_num, limit_Edge());
    rewards.resize(reward_edge_num, reward_Edge());
    
    dist_weight.resize(node_num, vector<int>(node_num, INF));
    dist_w_j.resize(node_num, vector<int>(node_num, INF));

    route_weight.resize(node_num, vector<int>(node_num,-1));

    SG_dist_w.resize(node_num, vector<int>(node_num, INF));
    SG_dist_j.resize(node_num, vector<int>(node_num, INF));

    must_node.resize(must_node_num);

    //printf("%d %d %d %d\n",node_num,edge_num,reward_edge_num,forbidden_edge_num);

    line += 2;
    sscanf(topo[line], "%d %d", &start , &end);

    line += 2;
    for (int i = 0; i<must_node_num ; ++i,++line)
    sscanf(topo[line], "%d", &must_node[i]);
    //printf("%d %d\n",must_node0,must_node1);
    
    ++line;
    for (int i = 0; i < edge_num; ++i, ++line) {
        sscanf(topo[line], "%d %d %d", &u, &v, &c);//链路起始节点、终止节点、总带宽、单位租用费
       // printf("%d %d %d\n",u,v,c);
        graph[u].emplace_back(v, c);
        graph[v].emplace_back(u, c);
    }
    
    //printf("finish edge input\n");

    ++line;
    for (int i = 0; i < forbidden_edge_num; ++i, ++line) {
        sscanf(topo[line], "%d %d %d", &u, &v, &c);
        //printf("%d %d %d \n",u,v,c);
        limits[i].t = v;
        limits[i].u = u;
        limits[i].c = c;
        set_cost(u,v,INF);
        set_cost(v,u,INF);
    }

    //printf("finish forbidden edge input\n");

    ++line;
    for (int i = 0; i < reward_edge_num; ++i, ++line) {
        sscanf(topo[line], "%d %d %d", &u, &v, &c);
        //printf("%d %d %d\n",u,v,c);
        rewards[i].t = v;
        rewards[i].u = u;
        rewards[i].c = c;
    }
    //printf("finish reward_edge_num\n");
}

int Graph::get_cost(int v, int id){
    for (vector<n_nEdge>::iterator it = graph[v].begin(); it !=graph[v].end(); ++it)
        if ( it->v == id )
            return it->c;
    return -1;
}

void Graph::set_cost(int v,int id,int cost){
    for (vector<n_nEdge>::iterator it = graph[v].begin(); it !=graph[v].end(); ++it)
        if( it->v == id )
            it->c = cost;
}

void Graph::Show_best_route(int s,int e){
    for(int i = e; route_weight[s][i] >= 0; i = route_weight[s][i]){
        printf("%d ", i);
    }
}

#define for_each_node for(int s = 0; s < node_num; ++s)
void Graph::SPFA_weight() {   //要把奖励边重新赋值一下，利用SPFA将所有点的最短路径找出来  带权
    for_each_node {             
        dist_weight[s][s] = 0;
        dist_w_j[s][s] = 1;
        deque<int> q;                
        memset(vis,0,sizeof(vis));
        vis[s] = 1;
        q.push_back(s);
        while (!q.empty()) {
            int u = q.front();
            q.pop_front();
            vis[u] = 0;
            for (unsigned int i = 0; i < graph[u].size(); ++i) {
                int v = graph[u][i].v;
                int dis =  dist_weight[s][u] + graph[u][i].c;
                int jump = dist_w_j[s][u] + 1;
                if (dis < dist_weight[s][v] || (dis == dist_weight[s][v] && jump < dist_w_j[s][v]))
                { 
                    dist_weight[s][v] = dis;
                    dist_w_j[s][v] = jump;
                    route_weight[s][v] = u;
                    if (!vis[v]) {
                        vis[v] = 1;
                        if (q.size () && dist_weight[s][v] < dist_weight[s][q[0]])
                            q.push_front(v);
                        else
                            q.push_back(v);
                    }
                }
            }
        }
    }
    for(vector<reward_Edge>::iterator it = rewards.begin(); it !=rewards.end();++it){
        int a = it->t;
        int b = it->u;
        dist_weight[a][b] = it->c;
        dist_weight[b][a] = it->c;
        dist_w_j[a][b] = 2;
        dist_w_j[b][a] = 2;
        route_weight[a][b] = a;
        route_weight[b][a] = b;
    }
}

void Graph::SPFA_For_SG() {
    for_each_node {             
        SG_dist_w[s][s] = 0;
        SG_dist_j[s][s] = 1;
        deque<int> q;                
        memset(vis,0,sizeof(vis));
        vis[s] = 1;
        q.push_back(s);
        while (!q.empty()) {
            int u = q.front();
            q.pop_front();
            vis[u] = 0;
            for (unsigned int i = 0; i < graph[u].size(); ++i) {
                int v = graph[u][i].v;
                int dis =  SG_dist_w[s][u] + graph[u][i].c;
                int jump = SG_dist_j[s][u] + 1;
                if ((dis < SG_dist_w[s][v] || (dis == SG_dist_w[s][v] && jump < SG_dist_j[s][v])) && (jump) <= limit_depth)
                { 
                    SG_dist_w[s][v] = dis;
                    SG_dist_j[s][v] = jump;
                    if (!vis[v]) {
                        vis[v] = 1;
                        if (q.size () && SG_dist_w[s][v] < SG_dist_w[s][q[0]])
                            q.push_front(v);
                        else
                            q.push_back(v);
                    }
                }
            }
        }
    }
    for(vector<reward_Edge>::iterator it = rewards.begin(); it !=rewards.end();++it){
        int a = it->t;
        int b = it->u;
        SG_dist_w[a][b] = it->c;
        SG_dist_w[b][a] = it->c;
        SG_dist_j[a][b] = 2;
        SG_dist_j[b][a] = 2;
    }
}
#undef for_each_node

DFS_Search::DFS_Search(Graph & graph){
    g = &graph;
    best_cost = INF;
    best_jump = INF;
    route.push(g->start);
    memset(g->vis,0,sizeof(g->vis));
}

void DFS_Search::DFS_withcutoff(int start){
    for (unsigned int i = 0; i < g->graph[start].size(); ++i) {
        if(route.size() > g->limit_depth-1) {continue;}
        if(sum_stack(route) > best_cost) {continue;}
        //if(1) {}
        int v = g->graph[start][i].v;
        route.push(v);
        if(v == g->end){
            if(isok(route)){
                if(sum_stack(route) < best_cost){
                    best_cost = sum_stack(route);
                    best_jump = route.size();
                    best_route = route;
                }
            }
        }
        else{
            DFS_withcutoff(v);
        }
        route.pop();
    }
}

int DFS_Search::sum_stack(stack<int> in){
    int sum = 0;
    while(1){
        int front = in.top();
        in.pop();
        if(in.empty()){break;}
        int behind = in.top();
        sum  = sum + g->get_cost(behind,front);
    }
    return sum;
}

bool DFS_Search::isok(stack<int> in){
    bool E_flag = 0;
    bool V_flag = 0;
    
    stack<int> E_temp = in;
    stack<int> V_temp = in;

    int V_vis[g->must_node_num];
    memset(V_vis,0,sizeof(V_vis));
    while(!V_temp.empty()){
        int t;
        t = V_temp.top();
#define for_eachmustnode for(int i = 0; i < g->must_node_num; ++i)
        for_eachmustnode
            if(t == g->must_node[i])
               V_vis[i]=1;
        V_temp.pop();
    }
    V_flag = 1;
    for_eachmustnode
        V_flag = V_flag & V_vis[i];
#undef for_eachmustnode
    
    int E_vis[g->reward_edge_num];
    memset(E_vis,0,sizeof(E_vis));
    while(!E_temp.empty()){
        int t;
        t = E_temp.top();
#define for_eachmustedge for(int i = 0;i < g->reward_edge_num; ++i)
        for_eachmustedge{
            if(t == g->rewards[i].t){
                E_temp.pop();
                int next = E_temp.top();
                if(next == g->rewards[i].u)
                    E_vis[i] = 1;
                E_temp.push(t);
            }
            if(t == g->rewards[i].u){
                E_temp.pop();
                int next = E_temp.top();
                if(next == g->rewards[i].t)
                    E_vis[i] = 1;
                E_temp.push(t);
            }
        }
        E_temp.pop();
    }
    E_flag = 1;
    for_eachmustedge
        E_flag = E_flag & E_vis[i];
#undef for_eachmustedge
    return V_flag & E_flag;
    //    return 1;
} 

Simplified_Graph_Search::Simplified_Graph_Search(Graph & graph){
    g = &graph;    

    best_cost = INF;         
    best_jump = INF;
    
    route.push(g->start);
    memset(g->vis,0,sizeof(g->vis));
}

void Simplified_Graph_Search::Simplified_Graph(){
    vector<int > new_node;
    new_node = g->must_node;
    for(vector<reward_Edge>::iterator it = g->rewards.begin(); it != g->rewards.end(); it++){
        new_node.emplace_back(it->t);
        new_node.emplace_back(it->u);
    }
    new_node.emplace_back(g->start);
    new_node.emplace_back(g->end);

    //去重
    sort(new_node.begin(),new_node.end());
    new_node.erase(unique(new_node.begin(), new_node.end()), new_node.end());
    g_node = new_node;

    simplified_g.graph.resize(g->node_num, vector<n_nEdge>());

    for(int i = 0;i < new_node.size(); ++i)
        for(int j = 0;j < new_node.size(); ++j)
        {   
            if( i != j )
                simplified_g.graph[new_node[i]].emplace_back(new_node[j],g->SG_dist_w[new_node[i]][new_node[j]]);
        }

    for(vector<reward_Edge>::iterator it = g->rewards.begin(); it != g->rewards.end(); ++it){
        int a = it->t;
        int b = it->u;
        int c = it->c;
        simplified_g.set_cost(a,b,c);
        simplified_g.set_cost(b,a,c);
    }
/*
    for(int i=0; i < simplified_g.graph.size();i++)
        for(int j=0; j<simplified_g.graph[i].size();j++)
        {
            printf("%d %d %d\n",i,simplified_g.graph[i][j].v,simplified_g.graph[i][j].c);
        }
*/
}

void Simplified_Graph_Search::TSP_taboo_search(){
    int nn;
    initGroup();                    //初始化编码Ghh
    copyGh(Ghh, bestGh);            //复制当前编码Ghh到最好编码bestGh

    /*
    printf("\n");
    for(int i=0;i<cityNum;i++){
        printf("%d ",Ghh[i]);
    }
    printf("\n");
    */

    bestEvaluation = evaluate(Ghh);
    bestJump = evaluate_J(Ghh);

    while(t < MAX_GEN){
        nn = 0;
        localEvaluation = INF;
        while( nn < N ){
            Linju(Ghh,tempGhh);             //得到当前编码Ghh的领域编码tempGhh
            if(panduan(tempGhh) == 0){      //判断编码是否在禁忌表中
                //不在
                tempEvaluation = evaluate(tempGhh);
                tempJump = evaluate_J(tempGhh);
                if( tempEvaluation < localEvaluation){
                    copyGh(tempGhh,LocalGhh);
                    localEvaluation = tempEvaluation;
                    localJump = tempJump;
                }
                nn++;
                //printf("%d\n",nn);
            }
            //printf("TEST\n");
        }
        //printf("***************\n");
        if(localEvaluation < bestEvaluation && isokGhh(LocalGhh) && localJump <= g->limit_depth){
            bestT = t;
            copyGh(LocalGhh,bestGh);
            bestEvaluation = localEvaluation;
            bestJump = localJump;
        }
        copyGh(LocalGhh,Ghh);

        jiejinji(LocalGhh);     //解禁忌表，LocalGhh加入禁忌表
        t++;
//        printf("TEST\n");
    }
}

void Simplified_Graph_Search::DFS_search(int start){
    for (unsigned int i = 0; i < simplified_g.graph[start].size(); ++i) {
        if(sum_jump(route) > g->limit_depth-1) {continue;}
        if(sum_stack(route) > best_cost) {continue;}
        //if(1) {}
        int v = simplified_g.graph[start][i].v;
        route.push(v);
        if(v == g->end){
            if(isok(route) && sum_jump(route)<=g->limit_depth){
                if(sum_stack(route) < best_cost){
                    best_cost = sum_stack(route);
                    best_jump = route.size();
                    best_route = route;
                }
            }
        }
        else{
            DFS_search(v);
        }
        route.pop();
    }
}

int Simplified_Graph_Search::sum_jump(stack<int> in){
    int sum = 1;
    while(1){
        int front = in.top();
        in.pop();
        if(in.empty()){break;}
        int behind = in.top();
        sum = sum + g->SG_dist_j[front][behind]-1;
    }
    return sum;
}

int Simplified_Graph_Search::sum_stack(stack<int> in){
    int sum = 0;
    while(1){
        int front = in.top();
        in.pop();
        if(in.empty()){break;}
        int behind = in.top();
        sum  = sum + simplified_g.get_cost(behind,front);
    }
    return sum;
}

bool Simplified_Graph_Search::isok(stack<int> in){
    bool E_flag = 0;
    bool V_flag = 0;
    
    stack<int> E_temp = in;
    stack<int> V_temp = in;

    int V_vis[g->must_node_num];
    memset(V_vis,0,sizeof(V_vis));
    while(!V_temp.empty()){
        int t;
        t = V_temp.top();
#define for_eachmustnode for(int i = 0; i < g->must_node_num; ++i)
        for_eachmustnode
            if(t == g->must_node[i])
               V_vis[i]=1;
        V_temp.pop();
    }
    V_flag = 1;
    for_eachmustnode
        V_flag = V_flag & V_vis[i];
#undef for_eachmustnode
    
    int E_vis[g->reward_edge_num];
    memset(E_vis,0,sizeof(E_vis));
    while(!E_temp.empty()){
        int t;
        t = E_temp.top();
#define for_eachmustedge for(int i = 0;i < g->reward_edge_num; ++i)
        for_eachmustedge{
            if(t == g->rewards[i].t){
                E_temp.pop();
                int next = E_temp.top();
                if(next == g->rewards[i].u)
                    E_vis[i] = 1;
                E_temp.push(t);
            }
            if(t == g->rewards[i].u){
                E_temp.pop();
                int next = E_temp.top();
                if(next == g->rewards[i].t)
                    E_vis[i] = 1;
                E_temp.push(t);
            }
        }
        E_temp.pop();
    }
    E_flag = 1;
    for_eachmustedge
        E_flag = E_flag & E_vis[i];
#undef for_eachmustedge
    return V_flag & E_flag;
    //    return 1;
}

void Simplified_Graph_Search::Show_Route(){
    printf("\n利用DFS寻找简化图满足条件路径：\n");
    int f = 0;
    stack<int> simple_route;
    simple_route = best_route;

    if(best_route.empty()){
        f = 1;
        goto no_answer1;
    }    
    printf("简化图dfs找的最短开销：%d\n", best_cost);     //best-_cost少了1
//    printf("简化图dfs找的最短跳数；%d\n", best_jump);
//    printf("简化图dfs找到的路径为：\n");
    
/*
    while(!simple_route.empty()){
        int a = simple_route.top();
        printf("%d ",simple_route.top());
        simple_route.pop();
        if(!simple_route.empty()){
            int b = simple_route.top();
            printf("(%d) ",g->SG_dist_j[a][b]);
        }
    }
*/    
    printf("原图的最短跳数为：%d\n",sum_jump(best_route));
    printf("转化为原图的路径为：\n");
    while(1){
        int behind = best_route.top();
        best_route.pop();
        if(best_route.empty()){ 
            printf("%d",behind);
            break;
        }
        int front = best_route.top();  //总的来说就是总费用+路径两个问题，总费用差了1 是reward边没
        g->Show_best_route(front,behind);  //因为这里放的是SPFA找的路径 可能导致会跳，但是权值为什么出现1的差别，暂时没有理解
    }
    printf("\n");

no_answer1:
    if (f)
        printf("当前图无解\n");
}

/*g:运行代数 c:每次搜索邻居个数 m:禁忌长度*/
void Simplified_Graph_Search::init_taboo(int g,int c,int m){
    cityNum = g_node.size();
    MAX_GEN = g;
    N = c;
    ll = m;
    
    Ghh.resize(cityNum);
    bestGh.resize(cityNum);
    bestEvaluation = INF;
    bestJump = INF;
    LocalGhh.resize(cityNum);
    localEvaluation = INF;
    localJump = INF;
    tempGhh.resize(cityNum);
    tempEvaluation = INF;
    tempEvaluation = INF;

    jinji.resize(ll,vector<int>(cityNum));
    bestT = 0;
    t = 0;
}

void Simplified_Graph_Search::initGroup(){      //初始化的时候其实可以考虑一些骚操作
    int i,j;
    //起点和重点固定
    Ghh[0] = g->start;
    Ghh[cityNum-1] = g->end;
    for(i = 1;i<cityNum-1;){       //编码长度，首位不管，这里要好好琢磨一下cityNum初值的问题
        struct timeval tv;
        gettimeofday(&tv, NULL);
        srand(tv.tv_usec);
        int random = rand() % cityNum; //这里要思考一下一个问题，mustnode会不会被装满的问题
        //printf("%d\n",g_node[random]);
        Ghh[i] = g_node[random];
        if(Ghh[i] == Ghh[0] || Ghh[i] == Ghh[cityNum-1])    //出去开头和结尾
            continue;
        for(j = 0;j < i; j++)       //去重
            if(Ghh[i] == Ghh[j])
                break;
        if(j == i)
            i++;
    }
}

void Simplified_Graph_Search::copyGh(vector<int> &Gha,vector<int> &Ghb){
    for(int i = 0;i < cityNum;i++)
        Ghb[i] = Gha[i];
}

void Simplified_Graph_Search::print(){
    printf("最佳长度出现迭代数：%d\n",bestT);
    printf("最佳长度：%d\n",bestEvaluation);
    printf("当前跳数：%d\n",bestJump);
    printf("原图最短路径：\n");
    for(int i=0;i< cityNum-1; i++){
        int front = bestGh[i];
        int behind = bestGh[i+1];
        g->Show_best_route(behind,front); 
    }
    printf("%d",bestGh[cityNum-1]);
}

int Simplified_Graph_Search::isokGhh(vector<int> &tempGh){
    int flag = 1;
    int flags[g->reward_edge_num];
    memset(flags,0,sizeof(flags));
    int n=0;
    for(vector<reward_Edge>::iterator it = g->rewards.begin();it!=g->rewards.end();it++,++n){
        int a = it->t;
        int b = it->u;
        for(int i = 0;i < cityNum;i++){
            if(tempGh[i] == a){
                if(tempGh[i-1] == b || tempGh[i+1] == b){
                    flags[n]=1;
                    break;
                }
            }
            if(tempGh[i] == b){
                if(tempGh[i-1] == a || tempGh[i+1] == a){
                    flags[n]=1;
                    break;
                }
            }
        }
    }
    for(int i=0;i < g->reward_edge_num;i++)
    {
        flag = flag && flags[i];
    }
    //printf("out:%d\n",flag);
    return flag;
}

int Simplified_Graph_Search::evaluate(vector<int> &chr){
    int len = 0;
    //0,1 + 1,2 ……+cityNum-2,cityNum-1+cityNum-1,cityNum
    //在初始复制的时候要考虑一下最初节点和最末节点的需求
    for(int i = 1; i < cityNum; i++)
        len += simplified_g.get_cost(chr[i-1],chr[i]);
    return len;
}

int Simplified_Graph_Search::evaluate_J(vector<int> &chr){
    int jump = 1;
    for(int i = 0; i < cityNum-1; i++){
        int front = chr[i];
        int behind = chr[i+1];
        jump = jump + g->SG_dist_j[front][behind]-1;
    }
    return jump;
}

void Simplified_Graph_Search::Linju(vector<int> &Gh,vector<int> &tempGh){
    int i, temp;
    int ran1 = 0, ran2 = 0;

    for(i = 0;i < cityNum; i++)
        tempGh[i] = Gh[i];

    struct timeval tv;
    gettimeofday(&tv, NULL);
    srand(tv.tv_usec);
    while(ran1 == 0 || ran1 == cityNum-1){
        gettimeofday(&tv, NULL);
        srand(tv.tv_usec);
        ran1 = rand() % cityNum;
    }
    while(ran2 == 0 || ran2 == cityNum-1){
        gettimeofday(&tv, NULL);
        srand(tv.tv_usec);
        ran2 = rand() % cityNum;
    }
    while(ran1 == ran2 || ran2 ==0 || ran2 == cityNum-1){ 
        gettimeofday(&tv, NULL);
        srand(tv.tv_usec);
        ran2 = rand() % cityNum;
    }
    //printf("out:%d %d\n",ran1,ran2);
    temp = tempGh[ran1];
    tempGh[ran1] = tempGh[ran2];
    tempGh[ran2] = temp;
}

int Simplified_Graph_Search::panduan(vector<int> &tempGh){
    int i, j;
    int flag = 0;
    for (i = 0; i < ll; i++) {
        flag = 0;
        for (j = 0; j < cityNum; j++) {
             if (tempGh[j] != jinji[i][j]) {
                flag = 1;// 不相同
                break;
             }
        }
        if (flag == 0)//相同，返回存在相同
            break;
    }
    if (i == ll)// 不等
        return 0;// 不存在
    else
        return 1;// 存在
}

void Simplified_Graph_Search::jiejinji(vector<int> &tempGh){
    int i,j,k;
    //删除禁忌表第一个编码，后面编码往前移动
    for(i = 0;i < ll-1;i++){
        for(j = 0;j < cityNum;j++)
            jinji[i][j] = jinji[i+1][j];
    }

    //新的编码加入禁忌表
    for(k = 0;k < cityNum;k++){
        jinji[ll - 1][k] = tempGh[k];
    }
}
