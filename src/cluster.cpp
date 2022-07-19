/**
 * Author      : Masrul Huda (mail2masrul@gmail.com)
 * Host        : iMacPro@Swalm
 * Created     : Tuesday Jul 05, 2022 10:01:02 CDT
 */

#include "gmx_traj.hpp"
#include <vector>
#include <map>
#include <cmath>
#include <iostream>
#include <limits>
#include <iomanip>

class Graph {
public:
    std::map<int,std::vector<int>>adjList;
    std::map<int,bool> visited; 
    std::vector<int>nodes;
    std::vector<std::vector<int>>connected_components; 
    std::vector<int>subgraph; 

    Graph(std::vector<int>nodes):nodes(nodes){
        for (int a: nodes){
            std::vector<int> _temp{};
            adjList[a] = _temp;
            visited[a] = false; 
        }
    }

    void add_edges(int a, int b){
        adjList[a].push_back(b);
        adjList[b].push_back(a);
    }

    void find_connected_component(){
        // set visited false 
        for (auto node : nodes)
            visited[node] = false;

        for (auto node : nodes){
            if (!visited[node]){
                subgraph.clear();
                _DFS(node);
                connected_components.push_back(subgraph);
            }
        }
    }

private: 
    void _DFS(int node){
        visited[node]=true; 
        subgraph.push_back(node);

        for (auto neigh : adjList[node])
            if (!visited[neigh]) _DFS(neigh);
    }
};


class Cluster {
public:
    std::string molecule_name;
    float rcut=0.35; // nm 
    int nClusters; 
    int max_cluster_size;
    std::vector<int> polymer_counts; // no of mono, di, trimer etc
    std::vector<std::vector<int>> clusters; 

    Cluster(std::string _molecule_name):molecule_name(_molecule_name){};
    ~Cluster();
    void set_rcut(float);
    void run(GMXTraj&);

private:
    std::vector<int> nodes;
    std::vector<std::vector<int>> edges; 
    void _create_edges(GMXTraj&);
    bool _is_neighbour(GMXTraj&,int,int); 
};

void Cluster::run(GMXTraj& traj){

    // Creates nodes for Graph 
    nodes.clear();
    auto& mol_trackers = traj.molecule_trackers;
    for (int i =0; i< traj.nmolecules; ++i)
        if (mol_trackers[i].name == molecule_name) 
            nodes.push_back(i);

    // Create edges for Graph 
    edges.clear();
    _create_edges(traj);

    // find connected components/clusters using Graph/DFS   
    Graph graph(nodes); 
    for (int i=0; i< edges.size();++i)
        graph.add_edges(edges[i][0],edges[i][1]);
    graph.find_connected_component();
    clusters = std::move(graph.connected_components);
    nClusters=clusters.size();

    // find max cluster size 
    max_cluster_size = 1; 
    for (auto clust : clusters)
        if (clust.size() > max_cluster_size)
            max_cluster_size  = clust.size();

    // Upto 5 mer counts 
    polymer_counts.clear();
    for (int i=0; i<=5; ++i)
        polymer_counts.push_back(0);

    for (auto clust : clusters)
        if (clust.size() <= 5)
            ++polymer_counts[clust.size()]; 
   
}

Cluster::~Cluster(){
    clusters.clear();
    nodes.clear();
    edges.clear();
}

void Cluster::set_rcut(float _rcut){
    rcut = _rcut;
}

void Cluster::_create_edges(GMXTraj& traj){
    auto& trackers = traj.molecule_trackers;

    #pragma omp  parallel 
    {
        // private variables to threads 
        std::vector<std::vector<int>> edges_private; 

        #pragma omp for schedule(dynamic) 
        for (int i =0; i< traj.nmolecules; ++i){
            if (trackers[i].name != molecule_name) continue;
            for (int j=i+1; j<traj.nmolecules; ++j){
                if (trackers[j].name != molecule_name) continue;
                if (_is_neighbour(traj,i,j))
                    edges_private.push_back(std::vector<int>{i,j});
            }
        }
       
        // Merging private copy of edges_private 
        #pragma omp critical  
        edges.insert(
            edges.end(), 
            std::make_move_iterator(edges_private.begin()), 
            std::make_move_iterator(edges_private.end())
        );
    } // end of omp parallel 

}

bool Cluster::_is_neighbour(GMXTraj& traj, int i, int j){
    auto& trackers = traj.molecule_trackers;
    auto& pos = traj.pos; 
    auto& box = traj.box; 

    int i_sIDx = trackers[i].sIDx;
    int i_eIDx = i_sIDx + trackers[i].natoms; 

    int j_sIDx = trackers[j].sIDx;
    int j_eIDx = j_sIDx + trackers[j].natoms; 
    float dr[3],dist;

    for (int ii=i_sIDx; ii<i_eIDx;++ii){
        for (int jj=j_sIDx; jj<j_eIDx; ++jj){
            dist = 0.0;
            for (int k=0;k<3;++k){
                dr[k] = pos[ii][k] - pos[jj][k];
                dr[k] = dr[k] - box[k][k]*round(dr[k]/box[k][k]);
                dist +=dr[k]*dr[k];
            }
            if (dist < rcut*rcut)return true; 
        }
    }
    return false; 
}

void print_help(){
    std::cout << "Cluster analysis\n";
    std::cout << "Author: Masrul Huda\n";
    std::cout << "Version: 0.1\n\n";
    
    std::cout << "    -f            *trajectory file [GRO,XTC,TRR]\n";
    std::cout << "    -s            topology file [gro] (required, if trajectory is XTC/TRR)\n";
    std::cout << "    --mol-name    *name of molecule\n";
    std::cout << "    --sys-info    information about molecule counts [ascii file]\n";
    std::cout << "    -rcut         cutoff distance between atoms\n";
    std::cout << "    -o            output file\n";
    std::cout << "    -dt           every dt ps  do analysis\n";
    std::cout << "    -b            begin analysis [ps]\n";
    std::cout << "    -e            end analysis [ps]\n";
    std::cout << "    --max-frames  maxium number of frames\n";
    std::cout << "    -h/--help     print this message\n";
    std::cout << "\n* required\n";

    exit(0);

}

int main(int argc, char* argv[]){  

    std::string top_file{""};
    std::string traj_file{""};
    std::string mol_name{""};
    std::string sys_info{""};
    std::string log_file{"result_clust.log"};
    float rcut=0.35;
    float dt =1.0;
    float begin = 0.0;
    float end = std::numeric_limits<float>::max();   
    float max_frames = std::numeric_limits<int>::max();   

    // Parse Command Line Argument 
    int i=1;
    while(i<argc && argv[i][0] == '-') {
        std::string opt{argv[i]} ;
        if(opt == "-f") traj_file=argv[++i] ;
        if(opt == "-s") top_file=argv[++i] ;
        if(opt == "--mol-name") mol_name=argv[++i] ;
        if(opt == "--sys-info") sys_info=argv[++i] ;
        if(opt == "-rcut") rcut = atof(argv[++i]);
        if(opt == "-o") log_file=argv[++i] ;
        if(opt == "-dt") dt = atof(argv[++i]);
        if(opt == "-b") begin = atof(argv[++i]);
        if(opt == "-e") end = atof(argv[++i]);
        if(opt == "--max-frames") max_frames=atoi(argv[++i]) ;
        if(opt == "-h" || opt == "--help") print_help();

        ++i ;
    }
    
    if (traj_file == "" || mol_name == "")
        print_help();


    GMXTraj traj{traj_file,top_file};
    traj.create_molecule_tracker(sys_info);

    Cluster cluster{mol_name};
    cluster.set_rcut(rcut);
    
    std::ofstream log; 
    log.open(log_file);
    log<<std::setw(12)<<"#   Time[ps]"
       <<std::setw(10)<<"nClusts"
       <<std::setw(10)<<"max_size"
       <<std::setw(10)<<"1mers"
       <<std::setw(10)<<"2mers"
       <<std::setw(10)<<"3mers"
       <<std::setw(10)<<"4mers"
       <<std::setw(10)<<"5mers"<<"\n";
    
    int nFrames=0; 
    float tol =0.00001; // tolarance for fmod
    while(traj.next()){
        if (nFrames >= max_frames || traj.time > end) break;
        if (traj.time < begin || std::fmod(traj.time,dt)>tol) continue; 

        std::cout << "Time: "<< traj.time << " nFrames: " << nFrames<<"\n"; 

        cluster.run(traj);
        log <<std::fixed<<std::setw(12)<<std::setprecision(3)<<std::right<<traj.time
            <<std::fixed<<std::setw(10)<<std::right<<cluster.nClusters
            <<std::fixed<<std::setw(10)<<std::right<<cluster.max_cluster_size
            <<std::fixed<<std::setw(10)<<std::right<<cluster.polymer_counts[1]
            <<std::fixed<<std::setw(10)<<std::right<<cluster.polymer_counts[2]
            <<std::fixed<<std::setw(10)<<std::right<<cluster.polymer_counts[3]
            <<std::fixed<<std::setw(10)<<std::right<<cluster.polymer_counts[4]
            <<std::fixed<<std::setw(10)<<std::right<<cluster.polymer_counts[5]
            <<"\n";
        ++nFrames;
    }
    
    return 0;
}

