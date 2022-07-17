/**
 * Author      : Masrul Huda (mail2masrul@gmail.com)
 * Host        : iMacPro@Swalm
 * Created     : Tuesday Jul 05, 2022 10:01:02 CDT
 */

#include "gmx_traj.hpp"
#include <vector>
#include<map>
#include <cmath>
#include <iostream>



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
    /*
     Algo for connected component from following source
     https://cp-algorithms.com/graph/search-for-connected-components.html
     */
        // set visited false 
        for (auto node : nodes)
            visited[node] = false;

        for (auto node : nodes){
            if (!visited[node]){
                subgraph.clear();
                DFS(node);
                connected_components.push_back(subgraph);
            }
        }
    }

    void DFS(int node){
        visited[node]=true; 
        subgraph.push_back(node);

        for (auto neigh : adjList[node])
            if (!visited[neigh]) DFS(neigh);
    }
};


class Cluster {
public:
    std::string molecule_name;
    float rcut=0.35; // nm 
    int nclusters; 
    std::vector<std::vector<int>> connected_components; 

    Cluster(std::string _molecule_name):molecule_name(_molecule_name){};
    ~Cluster();
    void set_rcut(float);
    void run(GMXTraj&);


private:
    std::vector<int> nodes;
    std::vector<std::vector<int>> edges; 
    void _create_edges(GMXTraj&);
};

void Cluster::run(GMXTraj& traj){
    _create_edges(traj);

    // find connected components 
    Graph graph(nodes); 
    for (int i=0; i< edges.size();++i)
        graph.add_edges(edges[i][0],edges[i][1]);
    graph.find_connected_component();
    connected_components = std::move(graph.connected_components);
    
    nclusters=connected_components.size();
    std::cout << "Cluster found: "<<nclusters<<"\n";

}

Cluster::~Cluster(){
    connected_components.clear();
    nodes.clear();
    edges.clear();
}

void Cluster::set_rcut(float _rcut){
    rcut = _rcut;
}


void Cluster::_create_edges(GMXTraj& traj){

    int i_sIDx, i_natoms, i_eIDx;
    int j_sIDx, j_natoms, j_eIDx;
    
    #pragma omp parallel for 
    for (int i =0; i< traj.nmolecules; ++i){
        if (traj.molecule_trackers[i].name != molecule_name){
            continue;
        }
        else {
            i_sIDx = traj.molecule_trackers[i].sIDx;
            i_natoms = traj.molecule_trackers[i].natoms; 
            i_eIDx = i_sIDx + i_natoms; 
            nodes.push_back(i);
        }

        for (int j=i+1; j<traj.nmolecules; ++j){
            if (traj.molecule_trackers[j].name != molecule_name){
                continue;
            }
            else {
                j_sIDx = traj.molecule_trackers[j].sIDx;
                j_natoms = traj.molecule_trackers[j].natoms; 
                j_eIDx = j_sIDx + j_natoms; 
            }
            
            bool is_neigh=false; 
            for (int ii=i_sIDx; ii<i_eIDx;++ii){
                for (int jj=j_sIDx; jj<j_eIDx; ++jj){
                    float dx = traj.pos[ii][0] - traj.pos[jj][0];
                    float dy = traj.pos[ii][1] - traj.pos[jj][1];
                    float dz = traj.pos[ii][2] - traj.pos[jj][2];
                    
                    // apply PBC 
                    dx = dx - traj.lx*round(dx/traj.lx);
                    dy = dy - traj.ly*round(dy/traj.ly);
                    dz = dz - traj.lz*round(dz/traj.lz);

                    float dist = dx*dx + dy*dy + dz*dz; 

                    if (dist < rcut*rcut){
                        is_neigh = true; 
                        break; 
                    }
                }
                if (is_neigh) break;
            }
            
            #pragma critical 
            if (is_neigh){
                std::vector<int> neigh{i,j};
                edges.emplace_back(std::move(neigh));
            }
        }
    }

}


void print_help(){
    std::cout << "Cluster analysis\n\n";
    
    std::cout << "    -f  *trajectory file [gro,xtc,trr]\n";
    std::cout << "    -s  topology file [gro]\n";
    std::cout << "    --mol-name  *name of molecule\n";
    std::cout << "    --sys-info  information about molecule counts [ascii file]\n";
    std::cout << "    -h/--help print this message\n";

    exit(0);

}


int main(int argc, char* argv[])
{  
    
    std::string top_file{""};
    std::string traj_file{""};
    std::string mol_name{""};
    std::string sys_info{""};
    float rcut=0.35;
    
    // Parse Command Line Argument 
    int i=1;
    while(i<argc && argv[i][0] == '-') {
        std::string opt{argv[i]} ;
        if(opt == "-f") traj_file=argv[++i] ;
        if(opt == "-s") top_file=argv[++i] ;
        if(opt == "--mol-name") mol_name=argv[++i] ;
        if(opt == "--sys-info") sys_info=argv[++i] ;
        if(opt == "-rcut") rcut = atof(argv[++i]);
        if(opt == "-h") print_help();
        if(opt == "--help") print_help();

        ++i ;
    }
    
    if (traj_file == "" || mol_name == "")
        print_help();


    GMXTraj traj{traj_file,top_file};
    traj.create_molecule_tracker(sys_info);

    Cluster cluster{mol_name};
    cluster.set_rcut(rcut);

    while(traj.next()){;
        cluster.run(traj);
    }
    
    return 0;
}

