from RecTAD import edge2adj
from RecTAD import HiCRecon
from RecTAD import TADCaller

import numpy as np

##读取原始数据并转换为密集矩阵形式
path_input = "Path to testdata/testdata.txt"
graph_edge = np.loadtxt(path_input)
chr='chr3'
resolution = 50000
graph_adj = edge2adj(graph_edge, chr = chr, resolution = resolution, reference = "hg19")

#调用HiCRecon进行重建
reconstructor = HiCRecon(graph_adj, dimension=128)
matrix_recon = reconstructor.getRecon()

#调用TADCaller进行TAD边界检测
caller = TADCaller(min_window=3, max_window=30)
boundaries, TADs = caller.fit(matrix_recon)


print("TAD boundaries:", boundaries)
