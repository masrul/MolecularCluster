#include <vector>
#include <map>
#include <cmath>
#include <iostream>
#include <limits>
#include <iomanip>
#include <queue>
#include "ctime"

#include "gmx_traj.hpp"
#include "distance.hpp"
#include "timer.hpp"

#if defined(_OPENMP)
#include <omp.h>
#endif

class Graph {
public:
    std::map<int,std::vector<int>>adjList;
    std::map<int,bool> visited; 
    std::vector<int>nodes;
    std::vector<std::vector<int>>connected_components; 
    std::vector<std::vector<int>>sequence_ref; 
    std::vector<int>sequence; 

    Graph(){}; 
    void set_nodes(std::vector<int>_nodes){
        nodes = _nodes;
        for (int node : nodes)
            adjList[node] = std::vector<int>{};
    }

    void add_edges(const std::vector<std::vector<int>>& edges){
        adjList.clear();
        for (const auto& edge : edges){
            adjList[edge[0]].push_back(edge[1]);
            adjList[edge[1]].push_back(edge[0]);
        }
    }

    void find_connected_components(){
        // set visited false 
        for (auto node : nodes)
            visited[node] = false;
        sequence.clear();

        for (auto node : nodes){
            if (!visited[node]){
                sequence.clear();
                _DFS(node);
                connected_components.push_back(sequence);
            }
        }
    }
    
    void traverse_bfs(int node){
        // set visited false 
        for (auto node : nodes)
            visited[node] = false;
        sequence_ref.clear();
        
        sequence_ref.clear();
        _BFS(node);
    }

private: 
    void _DFS(int node){
        visited[node]=true; 
        sequence.push_back(node);

        for (auto neigh : adjList[node])
            if (!visited[neigh]) _DFS(neigh);
    }

    void _BFS(int node){
        // Slightly modified BFS with ref node for each travel
        visited[node]=true; 
        std::queue<int> q;
        q.push(node);

        while (!q.empty()) {
            int current = q.front();
            q.pop();

            for (auto neigh : adjList[current]){
                if (!visited[neigh]){
                    sequence_ref.push_back(std::vector{current,neigh});
                    q.push(neigh);
                    visited[neigh] = true; 
                }
            }
        }
    }

};



class Cluster {
public:
    GMXTraj *traj;
    std::string molecule_name;
    float rcut=0.35; // nm 
    int nClusts; 
    int max_clust_size;
    std::vector<int> polymer_counts; // no of mono, di, trimer etc
    std::vector<std::vector<int>> clusters; 
    Graph graph{}; 


    Cluster(GMXTraj* _traj,std::string _molecule_name){
        molecule_name=_molecule_name;
        traj=_traj;
    }

    ~Cluster();
    void set_rcut(float);
    void run();
    void dump_max_clust(std::string); 

private:
    std::vector<int> nodes;
    std::vector<std::vector<int>> edges; 
    int _natoms_per_mol; 
    void _create_edges();
    bool _is_neighbour(int,int); 
    void _dump_clust(int clust_id,std::string dump_file);
    void _make_whole_molecule(int,float*);
};


void Cluster::run(){

    // Creates nodes for Graph first invokation 
    if ((nodes.size()==0)){
        auto& trackers = traj->molecule_trackers;
        for (int i =0; i< traj->nmolecules; ++i)
            if (trackers[i].name == molecule_name) 
                nodes.push_back(i);
        
        if (nodes.size()==0){
            std::cout<<"Molecule with name("<<molecule_name<<") is not found!\n";
            exit(0);
        }
        else {
            graph.set_nodes(nodes);
            _natoms_per_mol = trackers[nodes[0]].natoms; 
            std::cout<<"Molecule name: "<<molecule_name<<"\n";
            std::cout<<"Number of molecules found: "<<nodes.size()<<"\n";
            std::cout<<"Number of atoms per molecule: "<<_natoms_per_mol<<"\n";
#if defined(_OPENMP)
            size_t num_threads = omp_get_max_threads( );
            std::cout<<"Number of threads used: "<<num_threads<<"\n";
#endif 
            std::cout<<"\nStarting analysis ....\n\n";

        }
    }

    // Create edges for Graph 
    _create_edges();
    graph.add_edges(edges);

    // find connected components/clusters using Graph/DFS   
    graph.find_connected_components();
    clusters = std::move(graph.connected_components);
    nClusts=clusters.size();

    // find max cluster size 
    max_clust_size = 1; 
    for (auto clust : clusters)
        if (clust.size() > max_clust_size)
            max_clust_size  = clust.size();

    // Upto 5 mer counts 
    polymer_counts.clear();
    for (int i=0; i<=5; ++i)
        polymer_counts.push_back(0);

    for (auto clust : clusters)
        if (clust.size() <= 5)
            ++polymer_counts[clust.size()]; 
}

void Cluster::dump_max_clust(std::string prefix){
    
    // Find max cluster information 
    int max_clust_id; 
    int max_clust_size = 1; 
    for (int i=0; i<clusters.size();++i){
        if (clusters[i].size() > max_clust_size){
            max_clust_id = i;
            max_clust_size = clusters[i].size();
        }
    }
    
    std::string dump_file_name = prefix + "_"+std::to_string(traj->time) + ".gro";
    _dump_clust(max_clust_id,dump_file_name);
}

void Cluster::_dump_clust(int clust_id, std::string dump_file_name){

    auto& trackers = traj->molecule_trackers;
    auto& pos = traj->pos; 
    auto& box = traj->box;

    // fix first broken molecule in max cluster
    int ref_molid = clusters[clust_id][0];
    graph.traverse_bfs(ref_molid);
    int sIDx = trackers[ref_molid].sIDx;
    _make_whole_molecule(ref_molid,pos[sIDx]);

    // fixing rest of the molecules 
    for (auto& mol_pair : graph.sequence_ref ){
        sIDx = trackers[mol_pair[0]].sIDx; 
        _make_whole_molecule(mol_pair[1],pos[sIDx]);
    }

    std::ofstream dump_file;
    dump_file.open(dump_file_name); 

    // Write to file 
    int natoms_max_clust = _natoms_per_mol*clusters[clust_id].size();
    dump_file << "Created by Masrul Huda"<<std::endl;
    dump_file << std::setw(10) <<std::left<<natoms_max_clust<<std::endl; 

    int atomID = 0;
    for (auto mol_id : clusters[clust_id]){ 
        int sIDx = trackers[mol_id].sIDx;
        int eIDx = sIDx + _natoms_per_mol; 
        for (int i = sIDx; i<eIDx; ++i){
            ++atomID; 
            dump_file << std::setw(5) << std::right << traj->resids[i]; 
            dump_file << std::setw(5) << std::right << traj->resnames[i]; 
            dump_file << std::setw(5) << std::right << traj->symbols[i];
            dump_file << std::setw(5) << std::right << atomID;
            dump_file << std::fixed<< std::setw(8) << std::setprecision(3)<<std::right << pos[i][0];
            dump_file << std::fixed<< std::setw(8) << std::setprecision(3)<<std::right << pos[i][1];
            dump_file << std::fixed<< std::setw(8) << std::setprecision(3)<<std::right << pos[i][2];
            dump_file << std::endl; 
                    
        }
    }

    dump_file << std::fixed<< std::setw(10) << std::setprecision(5)<<std::right << box[0][0];
    dump_file << std::fixed<< std::setw(10) << std::setprecision(5)<<std::right << box[1][1];
    dump_file << std::fixed<< std::setw(10) << std::setprecision(5)<<std::right << box[2][2];
    dump_file << std::endl; 
}

void Cluster::_make_whole_molecule(int mol_id, float* pos_ref){
    auto& trackers = traj->molecule_trackers;
    auto& pos = traj->pos; 
    auto& box = traj->box; 
    float dr[3],dist;

    int sIDx = trackers[mol_id].sIDx;
    int eIDx = sIDx + _natoms_per_mol; 
    
    for (int i=sIDx; i<eIDx; ++i){
        for (int k=0;k<3;++k){
            dr[k] = pos[i][k] - pos_ref[k];
            if (dr[k] > 0.5*box[k][k])
                pos[i][k] -= box[k][k];
            else if (dr[k] < -0.5*box[k][k])
                pos[i][k] += box[k][k]; 
        }
    }
}

Cluster::~Cluster(){
    clusters.clear();
    nodes.clear();
    edges.clear();
}

void Cluster::set_rcut(float _rcut){
    rcut = _rcut;
}

void Cluster::_create_edges(){
    edges.clear();
    auto& trackers = traj->molecule_trackers;

    #pragma omp  parallel 
    {
        // private variables to threads 
        std::vector<std::vector<int>> edges_private; 

        #pragma omp for schedule(dynamic) 
        for (int i =0; i< traj->nmolecules; ++i){
            if (trackers[i].name != molecule_name) continue;
            for (int j=i+1; j<traj->nmolecules; ++j){
                if (trackers[j].name != molecule_name) continue;
                if (_is_neighbour(i,j))
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

bool Cluster::_is_neighbour(int iMol, int jMol){
    auto& trackers = traj->molecule_trackers;
    auto& pos = traj->pos; 
    auto& box = traj->box; 

    int i_sIDx = trackers[iMol].sIDx;
    int i_eIDx = i_sIDx + trackers[iMol].natoms; 

    int j_sIDx = trackers[jMol].sIDx;
    int j_eIDx = j_sIDx + trackers[jMol].natoms; 
    float dr[3],dist;

    // Check based on molecule  
    if (trackers[iMol].natoms > 5){

        // find i_max_dist 
        float i_max_dist = 0.0; 
        int jj = i_sIDx; 
        for (int ii=i_sIDx+1; ii<i_eIDx;++ii){
            for (int k=0;k<3;++k)
                dr[k] = pos[ii][k] - pos[jj][k];
            min_img_dr(dr,box);
            dist =dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
            if (dist > i_max_dist)
                i_max_dist = dist; 
        }
        i_max_dist = std::sqrt(dist);

        // find j_max_dist 
        float j_max_dist = 0.0; 
        jj = j_sIDx; 
        for (int ii=j_sIDx+1; ii<j_eIDx;++ii){
            for (int k=0;k<3;++k)
                dr[k] = pos[ii][k] - pos[jj][k];
            min_img_dr(dr,box);
            dist =dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
            if (dist > j_max_dist)
                j_max_dist = dist; 
        }
        j_max_dist = std::sqrt(dist);

        // find dist between mols 
        for (int k=0;k<3;++k)
            dr[k] = pos[i_sIDx][k] - pos[j_sIDx][k];
        min_img_dr(dr,box);
        float com_dist =dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        com_dist = std::sqrt(com_dist);
        
        if (com_dist-i_max_dist -j_max_dist > rcut)
            return false; 
    }

    // checked based on atoms 
    for (int ii=i_sIDx; ii<i_eIDx;++ii){
        for (int jj=j_sIDx; jj<j_eIDx; ++jj){
            dist = 0.0;
            for (int k=0;k<3;++k)
                dr[k] = pos[ii][k] - pos[jj][k];
            min_img_dr(dr,box);
            dist =dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
            if (dist < rcut*rcut)return true; 
        }
    }
    return false; 
}


void print_help(){
    std::cout << "Cluster analysis\n";
    std::cout << "Author: Masrul Huda\n";
    std::cout << "Version: 0.1\n\n";
    
    std::cout << "    -f                 *trajectory file [GRO,XTC,TRR]\n";
    std::cout << "    -s                  topology file [gro] (required, if trajectory is XTC/TRR)\n";
    std::cout << "    --mol-name          *name of molecule\n";
    std::cout << "    --sys-info          information about molecule counts [ascii file]\n";
    std::cout << "    -rcut               cutoff distance between atoms\n";
    std::cout << "    -o                  output file\n";
    std::cout << "    -dt                 every dt ps  do analysis\n";
    std::cout << "    -b                  begin analysis [ps]\n";
    std::cout << "    -e                  end analysis [ps]\n";
    std::cout << "    --max-frames        maxium number of frames\n";
    std::cout << "    --dump-max-clust    set flag to dump max cluster\n";
    std::cout << "    --max-clust-prefix  Prefix of max clust files\n";
    std::cout << "    -h/--help           print this message\n";
    std::cout << "\n* required\n";

    exit(0);
}


int main(int argc, char* argv[]){  
    //Start a timer 
    Timer timer{};

    std::string top_file{""};
    std::string traj_file{""};
    std::string mol_name{""};
    std::string sys_info{""};
    std::string log_file{"result_clust.log"};
    std::string max_clust_prefix{"max_clust"};
    float rcut=0.35;
    float dt =1.0;
    float begin = 0.0;
    float end = std::numeric_limits<float>::max();   
    float max_frames = std::numeric_limits<int>::max();   
    bool dump_max_clust = false; 

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
        if(opt == "--dump-max-clust") dump_max_clust=true;
        if(opt == "--max-clust-prefix") max_clust_prefix=argv[++i];
        if(opt == "-h" || opt == "--help") print_help();

        ++i ;
    }
    
    if (traj_file == "" || mol_name == ""){
        std::cout<<"Required options missing!\n";
        std::cout<<"Try: cluster -h/--help for details.\n";
        exit(0);
    }

    // Create trajectory 
    GMXTraj traj{traj_file,top_file};
    traj.create_molecule_tracker(sys_info);
    
    // Create cluster analysis instance 
    Cluster cluster{&traj,mol_name};
    cluster.set_rcut(rcut);
    
    // Write to log 
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
    
    // Traverse trajectory and do analysis 
    int nFrames=0; 
    float tol =0.00001; // tolarance for fmod
    while(traj.next()){
        if (nFrames >= max_frames || traj.time > end) break;
        if (traj.time < begin || std::fmod(traj.time,dt)>tol) continue; 

        cluster.run();
        if (dump_max_clust)
            cluster.dump_max_clust(max_clust_prefix);
        std::cout 
            <<"Time: "<< traj.time 
            <<" nFrames: "<< nFrames
            <<" nCluts: "<< cluster.nClusts
            <<"\n"; 

        log <<std::fixed<<std::setw(12)<<std::setprecision(3)<<std::right<<traj.time
            <<std::fixed<<std::setw(10)<<std::right<<cluster.nClusts
            <<std::fixed<<std::setw(10)<<std::right<<cluster.max_clust_size*76
            <<std::fixed<<std::setw(10)<<std::right<<cluster.polymer_counts[1]
            <<std::fixed<<std::setw(10)<<std::right<<cluster.polymer_counts[2]
            <<std::fixed<<std::setw(10)<<std::right<<cluster.polymer_counts[3]
            <<std::fixed<<std::setw(10)<<std::right<<cluster.polymer_counts[4]
            <<std::fixed<<std::setw(10)<<std::right<<cluster.polymer_counts[5]
            <<"\n";
        ++nFrames;
    }

    // Print elapsed time 
    timer.print_elapsed_time();
    
    return 0;
}
 
