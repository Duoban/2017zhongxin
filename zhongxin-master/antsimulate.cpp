#include "antsimulate.h"
#include <set>
#include <map>
#include <utility>
#include <assert.h>

#define ProNode (1.0)
#define alpha   (0.005)
#define minpha  (0.2)
#define maxpha  (0.8)
#define minix   (0.000001)
#define beta    80
#define FULL_GRADE 100

using namespace std;

Graph *a_Graph;
int limit_path_depth;
set<int> must_node_set;
set<pair<int,int> > must_edge_set;

void Init_a_Graph(Graph &G) //初始化本地数据
{
    a_Graph = &G;
    for(int u = 0;u < G.node_num; u++){
        double cost_sum = 0.0;
        for(vector<n_nEdge>::iterator it = a_Graph->graph[u].begin();it != a_Graph->graph[u].end(); it++){
            //it->pro = ProNode/a_Graph->graph[u].size();
            cost_sum += 1.0/(it->c+1);
        }
        for(vector<n_nEdge>::iterator it = a_Graph->graph[u].begin();it != a_Graph->graph[u].end(); it++){
            it->pro = (1.0/(it->c+1))/cost_sum;
        }
    }
    must_node_set.insert(a_Graph->start);
    must_node_set.insert(a_Graph->end);
    for(vector<int>::iterator it = a_Graph->must_node.begin();it != a_Graph->must_node.end(); it++){
        must_node_set.insert(*it);
    }
    for(vector<reward_Edge>::iterator it = a_Graph->rewards.begin(); it != a_Graph->rewards.end(); it++){
        must_node_set.insert(it->t);
        must_node_set.insert(it->u);
        must_edge_set.insert(make_pair(it->t,it->u));
        must_edge_set.insert(make_pair(it->u,it->t));
    }
    limit_path_depth = a_Graph->limit_depth;
}

int GetNextNode(int node,int &cost) //根据概率值选择下一跳
{
    double rand_pro = (double(rand()%1000))/1000.0;
    for(vector<n_nEdge>::iterator it = a_Graph->graph[node].begin();it != a_Graph->graph[node].end();it++){
        rand_pro -= it->pro;
        if(rand_pro<0.000001){
            cost += it->c;
            return it->v;
        }
    }
    assert(rand_pro<0);
    return 0;
}

void UpdateEdge(int u,int v)
{
    double sum_pro = 0.0;
    int outdegree = a_Graph->graph[u].size();
    for(vector<n_nEdge>::iterator it = a_Graph->graph[u].begin(); it != a_Graph->graph[u].end(); it++){
        if(it->v != v){  //其他边减少概率值
            if(it->pro  - alpha /(outdegree-1) >  minpha / outdegree){
               it->pro -= alpha /(outdegree-1);
            }
            sum_pro += it->pro;
        }
        else {
            if(it->pro + alpha < maxpha ){
               it->pro += alpha; //对应边增加概率值
            }
            sum_pro += it->pro;
        }
    }
    for(vector<n_nEdge>::iterator it = a_Graph->graph[u].begin(); it != a_Graph->graph[u].end(); it++){
        it->pro = it->pro / sum_pro;
    }
    
}

void UpdatePath(vector<int> &path) //更新权值
{
    int u,v;
    u = *(path.begin());
    for(vector<int>::iterator it = path.begin() + 1;it != path.end(); it++){
        v = *it;
        UpdateEdge(u,v);
        u = v;
    }
}
double GetPro(int u,int v) //获得路径权值
{
    for(vector<n_nEdge>::iterator it = a_Graph->graph[u].begin(); it != a_Graph->graph[u].end();it++){
        if(it->v == v){
            return it->pro;
        }
    }
    assert(0);
    return 0;
}

void PrintPath(vector<int> &path)//打印路径
{
    for(vector<int>::iterator it = path.begin();it != path.end()-1; it++){
        printf("%d < %.3lf> ",*it,GetPro(*it,*(it+1)));
    }
    printf("%d\n",*(path.end()-1));
}

int GetPathCost(vector<int> &path) //计算路径代价
{
    int sum_cost = 0;
    int u,v;
    u = *path.begin();
    for(vector<int>::iterator it = path.begin()+1;it != path.end(); it++){
        v = *it;
        sum_cost += a_Graph->get_cost(u,v);
        u = v;
    }
    return sum_cost;
}
int  mincost = 2000000;
int  min_path_len = 100;
vector<int> SavedPath;
int EvaluatePath(vector<int> &path,int cur_set,int cost)
{
    set<int> nodeset;
    set<pair<int,int> > edgeset;
    int u,v;
    u = *path.begin();
    nodeset.insert(u);
    for(vector<int>::iterator it = path.begin()+1; it != path.end(); it++){
       v = *it;
       edgeset.insert(make_pair(u,v));
       edgeset.erase(make_pair(v,u)); 
       u = v;
    }
    int nodes_count = cur_set;
    int path_cost = cost;
    int edges_count = 0;
    for(set<pair<int,int> >::iterator it = edgeset.begin(); it != edgeset.end(); it++){ //路径边和必须边去重取交集
        if(must_edge_set.find(*it) != must_edge_set.end()){
            edges_count ++;
        }
    }
    if(nodes_count == must_node_set.size() && edges_count * 2 == must_edge_set.size()&&path_cost < 2 * mincost && path.size() <= limit_path_depth ){   //满足经过必须点和必须边且长度满足
        if(path_cost < mincost){
           mincost = path_cost;
           SavedPath.clear();
            SavedPath.resize(path.size());
           std::copy(path.begin(),path.end(),SavedPath.begin());
        }
        if(path.size() < min_path_len){
           min_path_len = path.size();
        }
        return FULL_GRADE;
    }
    return 0;
}
void FindPath(int src,int dst)
{
    int cur,next,cost;
    int cur_set,rest_node;
    int min_len = must_node_set.size();
    int min_set = (must_node_set.size() % 5==0) ? (must_node_set.size() * 0.8) : ( must_node_set.size() * 0.8 + 1);
    bool prune_flag;
    vector<int> path;
    set<int> path_node;
    while(true){
        cur = src;
        path.push_back(cur);
        path_node.insert(cur);
        cost = 0;
        prune_flag = false;
        cur_set = 1;
        rest_node = limit_path_depth + 1;
        while(rest_node-- > 0){  //寻找不超过限制跳数的路径（剪枝），松弛约束为限制跳数+1
              next = GetNextNode(cur,cost);
              cur = next;
              if(path_node.find(next) == path_node.end() && must_node_set.find(next) != must_node_set.end()){
                cur_set++; //增加新的必过点
              }
              path_node.insert(next);
              path.push_back(next);
              if(cur_set + rest_node < min_set ){ //剪枝条件：无法满足必过点要求
                   prune_flag = true;
                   break;
              }
              if(cost > mincost * 1.5){ //剪枝条件: 无法满足过小开销
                   prune_flag = true;
                   break;
              }
              if(cur == dst && path.size() >= min_len){//还需满足长度不小于必须经过点的个数
                 break;
              }
        }
        if(prune_flag == true){ //因为剪枝退出
            path.clear();
            path_node.clear();
            continue;
        }
        if(cur == dst && cur_set >= min_set){ //找到合法路径
               break;
        }
        else {  //路径过长没找到合法路径
            path.clear();
            path_node.clear();
            continue;
        }
    }
    int grade = EvaluatePath(path,cur_set,cost); //对路径进行评估
    //printf("cost: %d cur_set: %d min_set: %d\n",cost,cur_set,min_set);
    //PrintPath(path);
    UpdatePath(path);
    if(grade == FULL_GRADE){ //若满足所有要求，再次增加路径权值
        //printf("cost: %d  ",cost);
        //PrintPath(path);
        UpdatePath(path);
    }
}
void ant_simulate(Graph &G){
    Init_a_Graph(G); //初始化参数
    const int ant_num = 1000; //100万蚂蚁
    int id = 0;
    while(id++ < ant_num){  //寻找100万路径
        FindPath(G.start,G.end);
    }
    printf("over!\n");
    if(SavedPath.size() != 0){
    printf("最少开销: %d\n",GetPathCost(SavedPath));
    printf("最短路径：%d\n",min_path_len);
    printf("原图最少开销路径: ");
    for(vector<int>::iterator it = SavedPath.begin(); it != SavedPath.end(); it++){
        std::cout<<*it<<" ";
    }
    std::cout<<endl;
    }
    else {
        std::cout<<"no solution!\n"<<std::endl;
    }
}
