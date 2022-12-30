# Multi-clustering via Non-negative Matrix Tri-factorization

Multi-view (or multi-modal) data integration via Non-negative Matrix Tri-factorization based on data restrictions with the purpose of simultaneous clustering of the samples and features of the data.

Suppose that multiple datasets are available, some of which share the same samples (e.g. genomics and transcriptomics on the same samples), and other share the same features (e.g. genomics on two separate studies, or sample sets). The goal of multi-clustering is to simultaneously cluster the samples and features of each dataset based on restrictions and descriptives of the data. Categories of restrictions on the data-views (or modalities) are the following:

1. Same samples, i.e. same row-clusters

2. Same relationships between row-clusters and column-clusters

3. Same features, i.e. same column-clusters

## Non-negative Matrix Factorisation

In order to identify clusters simultaneously for all rows (samples) and columns (features) of the given datasets, the proposed solution is based on Non-negative Matrix Factorisation (NMF), which essentially splits a data matrix into the multiplication of two non-negative matrices that share one dimension. An extension of NMF, called Non-negative Matrix Tri-Factorisation (NMTF) splits a data matrix into the multiplication of three non-negative matrices which allows the second restriction from above, to be considered in the equation. 

Multiple extensions and variations have been proposed for both NMF and NMTF, including multi-view, i.e. splitting multiple data-views (or modalities) simulataneously according to underlying assumptions (e.g. common samples, or features). Some examples are the following:

- **Single latent representation on all views (Collective NMF - CNMF)**:  Akata Z.,  Thurau, C. and Bauckhage C.. (2011) Non-negative Matrix Factorization in Multimodality Data for Segmentation and Label Prediction. 16th Computer Vision Winter Workshop
  > Condition is usually too strict for applications on real data

- **Common consensus view (Joint-NMF)**: Jialu Liu, Chi Wang, Jing Gao, and Jiawei Han, Multi-View Clustering via Joint Nonnegative Matrix Factorization, Proceedings of the 2013 SIAM International Conference on Data Mining. 2013, 252-260
  > Clustering performance is affected by the quality of the views (sparsity, missing data etc.)i.e. consensus clustering is likely to degrade when incorporating low quality views

- **Specified consensus view (Co-regularised NMF-CoNMF)**: He X, Kan M-Y, Xie P, Chen X. Comment-based multi-view clustering of web 2.0 items. , ACM Press; 2014. pp. 771-782
  > Less strict and less sensitive to the quality of the views

## Restrictive Multi-NMTF

The proposed solution to this problem is the minimisation of the following equation:

<img width="891" alt="image" src="https://user-images.githubusercontent.com/45739717/210072155-69e2fde4-d1c7-44bb-a0ab-f57139ecf11d.png">

The optimised updates of the above equation are given below:

<img width="709" alt="image" src="https://user-images.githubusercontent.com/45739717/210073750-60bdfee1-f5de-4853-a48c-6e08ec0febb8.png">

