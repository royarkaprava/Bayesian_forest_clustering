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

"getPointEstC" is a function that produces the point estimate for clustering from a pairwise matrix.

"forestClust" is a simplified wrapper function for fitting the Bayesian forest model for one dataset.

## Instructions for use

The code is written in pure R and does not require any additional compilation.

Figure 1 is an illustrative plot that can be generated using "clusterSim" R package, and does not require fitting of the Bayesian spanning forest model.


