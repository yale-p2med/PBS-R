# PBS
Defining pathway based distance score to identify patient heterogeneity using gene expression data.

## Overview
Distance based unsupervised clustering of the gene expression data is commonly used to identify heterogeneity in samples. Due to the high noise level in the gene expression data and correlation between genes, traditional distances such as Euclidean distance may not be appropriate at representing the true biological differences between samples. Pre-defined biological pathways have been associated with disease phenotypes in a way that different genes in the same pathway can be perturbed in different subjects that have similar clinical outcomes. Therefore assessing the biological differences between samples based on pathways instead of single genes may enhance the robustness and accuracy of the clustering results. 

PBS is a software package that assesses the biological differences between samples using gene expression data by assuming that ontologically defined biological pathways in biologically similar samples have similar behavior. This package includes three parts: simulation, pathways and pbs_distance.

##References
Yan X, Chu JH, Gomez J, Koenigs M, Holm C, He X, Perez MF, Zhao H, Mane S, Martinez FD, Ober C, Nicolae DL, Barnes KC, London SJ, Gilliland F, Weiss ST, Raby BA, Cohn L, Chupp GL. 2015. Noninvasive analysis of the sputum transcriptome discriminates clinical phenotypes of asthma. Am J Respir Crit Care Med. 191(10):1116-25. 
