# MolecularCluster

This tool kit finds the number of clusters/aggregations form during molecular simulation. 

## Options 

    -f                 *trajectory file [GRO,XTC,TRR];
    -s                  topology file [gro] (required, if trajectory is XTC/TRR)
    --mol-name          *name of molecule
    --sys-info          information about molecule counts [ascii file]
    -rcut               cutoff distance between atoms
    -o                  output file
    -dt                 every dt ps  do analysis
    -b                  begin analysis [ps]
    -e                  end analysis [ps]
    --max-frames        maxium number of frames
    --dump-max-clust    set flag to dump max cluster
    --max-clust-prefix  Prefix of max clust files
    -h/--help           print this message
                        * required
                        
## Installation
```bash 
git clone https://github.com/masrul/MolecularCluster/
cd MolecularCluster
mkdir build && cd MolecularCluster
cmake ..
make && make install
```
