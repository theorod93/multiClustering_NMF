# Multi-clustering via Non-negative Matrix Tri-factorization

Multi-view (or multi-modal) data integration via Non-negative Matrix Tri-factorization based on data restrictions with the purpose of simultaneous clustering of the samples and features of the data.

Suppose that multiple datasets are available, some of which share the same samples (e.g. genomics and transcriptomics on the same samples), and other share the same features (e.g. genomics on two separate studies, or sample sets). The goal of multi-clustering is to simultaneously cluster the samples and features of each dataset based on restrictions and descriptives of the data. Categories of restrictions on the data-views (or modalities) are the following:

1. Same samples, i.e. same row-clusters

2. Same relationships between row-clusters and column-clusters

3. Same features, i.e. same column-clusters

## Restrictive Multi-NMTF

The proposed solution to this problem is the minimisation of the following equation

% h<sub>&theta;</sub>(x) = &theta;<sub>o</sub> x + &theta;<sub>1</sub>x
