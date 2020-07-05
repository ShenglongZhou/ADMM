This matlab code was updated at 05/07/2020 by Shenglong Zhou and aims at solving two types of sparse and low-rank matrices recovery. 

The solver is programmed based on the algorithm proposed in the paper:
S. Zhou, N. Xiu, Z. Luo and L. Kong, (2015),  Sparse and Low-Rank Covariance Matrix Estimation, Journal of the Operations Research Society of China, 3(2): 231-250.

[1] Folders
   (a) 'solver' contains the main solver 'ADMM.m'.
   (b) 'funcs' contains several useful functions including: 
           (b.1) 'Approx_rank' calculates the approximate rank of a matrix defined as in above paper.
           (b.2) 'Examples.m' generates two types examples related to sparse and low-rank matrices.
           (b.3) 'FTRate.m' computes  the false positive rate and the true positive rate.
           (b.4) 'SubGraph.m' helpes 'Graph.m' to generate graphs.
[2] 'Graph.m'  generates the graphs shown in above paper.
[3] 'Demon.m'  tests two examples and generating related data shown in tables in above paper.


Run 'Graph.m' or 'Demon.m' to see the performance of the solver. 
Please give credits to this paper if you use the codes for your research.




