# RecTAD

## 1. Introduction
RecTAD is a single-cell TAD-like domain identification method based on Hi-C contact matrix reconstruction. This method is grounded in graph theory and treats the single-cell Hi-C contact matrix as the adjacency matrix of a graph. First, sparse matrix factorization is employed to obtain preliminary embeddings of nodes in the entire graph. Then, in the embedding space, graph spectral augmentation theory is applied to enhance the representation of community structures in the graph. Based on the enhanced embeddings, the Hi-C contact matrix is reconstructed. Finally, TAD-like domains are identified on the reconstructed Hi-C contact matrix using multi-scale Haar features.

## 2. Example

**2.1 import modules**
```
from RecTAD import edge2adj
from RecTAD import HiCRecon
from RecTAD import TADCaller

import numpy as np
```

**2.2 prepare data**
```
##read original data and convert it to dense format
path_input = "C:/Users/Lenovo/Desktop/testdata.txt"
graph_edge = np.loadtxt(path_input)
chr='chr3'
resolution = 50000
graph_adj = edge2adj(graph_edge, chr = chr, resolution = resolution, reference = "hg19")
```

**2.3 reconstruct Hi-C map**
```
#use HiCRecon to reconstruct Hi-C map
reconstructor = HiCRecon(graph_adj, dimension=128)
matrix_recon = reconstructor.getRecon()
```

**2.4 call TAD-like domains**
```
#use TADCaller to call TAD-like domains

caller = TADCaller(min_window=3, max_window=50)
boundaries, TADs = caller.fit(matrix_recon)

print("TAD boundaries:", boundaries)
```
