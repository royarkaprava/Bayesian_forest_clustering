# Spectral Clustering, Bayesian Spanning Forest, and Forest Process

# Author Contributions Checklist Form

## Data

### Simulations:
Each simulation R script contains the steps to generate the data.

### Data applications:
The "penguins" dataset can be found via the R package "palmerpenguins" from CRAN.

The high-dimensional image clustering application data is stored in the "yaleB10subjects.Rda" file.

The fMRI neuroscience application data is stored in the "fmri_labels.Rda" file. The spatial coordinates for the regions of interest (ROIs) are in the "centroid.txt", and the ROI names are in "regionname.txt".

## Code

### Source code:
All of the source code for the Bayesian spanning forest model is stored in the "forest_class.R" file.

"Forest" is an R reference class defined to store a Forest object. Each forest object encapsulates a copy of data, parameters, one Markov chain, as well as the functions to perform Markov chain updates.

"extractC" is a helper function to extract the clustering point estimate from a spanning forest adjacency matrix.

"sample_eta_star", "sampleZ", "sampleZ_k", and "sampleV" are functions used for multi-view clustering that involves combining and updating information based on a list of Forest objects. The "multiview_fmri_application/multiview_fmri_application.Rmd" file contains a use case.

"matchAtoB" is a helper function that runs the Hungarian algorithm that permutes the cluster index in "clusteringA", so that its match rate to "clusteringB" is maximized.

"clusteringAccu" is a function that calculates the maximized match rate (clustering accuracy).

"getCoAssignmentMat" is function that calculates the coassignment probability matrix.

"getPointEstC" is a function that produces the point estimate for clustering from a pairwise matrix.

"forestClust" is a simplified wrapper function for fitting the Bayesian forest model for one dataset.

## Instructions for use

The code is written completely in R and does not require any additional compilation.

### Computing Time:
Each fitting of a Bayesian spanning forest model (1000 iterations of MCMC) takes about 5 minutes. The multi-view clustering and multiple-chain versions will take longer time, as they involve updating multiple forests.

### Script or Markdown that were used to generate figures and table in the paper.

Figure 1 is an illustrative plot that can be generated using "clusterSim" R package, and does not require fitting of the Bayesian spanning forest model.

Figure 2 can be generated using "compare_eigenvectors.r".

Figure 3 can be generated using "uq_near_manifold.r".

Figure 4 and 5 can be generated using "multiview_fmri_application/multiview_fmri_application.Rmd".

Figure S.1, S.2 and Table S.1 can be generated using "high_dim_application/uq_high_dimesion_clustering.r".

Figure S.3 can be generated using "covariate_dependent_clustering/penguins.Rmd".

Figure S.4 and S.5 can be generated using "uq_mixture_model_data_t_dist.r".

Figure S.6 can be generated using "uq_near_manifold.r".

Figure S.7, S.8, S.9 can be generated using "simulations/imbalanced_vs_clusteringErr.Rmd", "simulations/num_clusters_vs_clusteringErr.Rmd", and "simulations/loosely_connected_graphs_vs_clusteringErr.Rmd". For the convenience of review, the repeated experiment results are stored in the "RDa" binary files. 

Figure S.10, S.11, S.12 can be generated using "simulations/mixing_diagnostics_multiple_initialization.Rmd". 

Figure S.13 can be generated using "multiview_fmri_application/multiview_fmri_application.Rmd".

Figure S.14 can be generated using "uq_mixture_model_data_t_dist.r".

Figure S.15 can be generated using "compare_eigenvectors.r".


