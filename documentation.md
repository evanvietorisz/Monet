## Table of contents

- The basic pipeline
- How to change the pipeline to load mutant C. elegans embryos
- How to change the pipeline to load fish embryos
- How to do clustering
- Function documentation
- Key functions in the pipeline
- Other helpful functions
- Aux and helper
- Outdated or obsolete

---

## Basic pipeline overview

The most basic functionality of the pipeline is converting a .zip file containing annotated embryo imaging data into a color-coded visualization of the time-cells in the dataset. This takes place in `celeganstest.m`, the parent script for the project. This section outlines the flow of information through the script.

### Loading data

The user starts by specifying the directory of a .zip file (`emb`) that contains the imaging data. This data is loaded by the script `loadcells_unnamed.m`, starting at the first frame in the raw data and extending through the frame `tend`, which is specified as an input parameter. `tend` should generally refer to the last time where cell positions. `loadcells_unnamed.m` outputs two representations of the data as Matlab objects: `cells`, which is a cell array containing information about each cell in the dataset and is not used in this project, and `embdat`, which is a struct array in which every frame in the raw data is represented as a struct with fields for the positions of cells, names of cells, the succession of those cells, and other information. We call data of this form an "embryo struct". This is the closest thing to the raw data that the user can interact with in the pipeline.

### Alignment

The next step is to align input data to a reference dataset so that different embryos share a common coordinate system. This is done by `coalignNamedEmbryosPerTime2.m`, which computes the best fit rotation between a set of landmark cells in an input embryo (`embdat`) and a reference embryo, which is loaded from the local Matlab path under the name `ref_emb.mat`. This is done at every frame `embdat`. The function outputs `embdat_stabilized`, which is the same as `embdat` except that the positions of the cells are aligned with the reference.

`embdat` and reference embryo structs need not be the same length, but should cover roughly the same portion through development. The alignment matches the different timepoints in the datasets against each other using the time at which both the input and reference are four cells, so if the input dataset starts drastically later than that, some adjustment needs to be made to the alignment script.

If you are interested in aligning an input embryo through a certain point through its development, but not after, you could either make manual alterations to `coalignNamedEmbryosPerTime2.m` or split the original `embdat` into two variables, align one, not align the other, and then concatenate them afterwards. Doing so would still likely require edits to the alignment script/work to make sure the reference embryo was a good match for the truncated part of development.

Note that the current reference embryo is a C. elegans embryo, so do not do this for a Zebrafish embryo, and manually inspect `embdat_stabilized` if applying the process to a mutant embryo whose phenotype differs drastically from the reference embryo.

Note also that if you are interested in simply removing movement within a single embryo, the function `internallyAlignNamedEmbryo.m` is closer to what you need and can probably be adapted to fit your purpose. 

### Feature computation

The next step is to compute features, which is a process of converting an `embdat_stabilized` into a feature space in which rows correspond to cells from many frames in `embdat_stabilized` and columns contain spatiotemporal features expressing those cells’ states.

Not all frames from `embdat_stabilized` make it into the feature space. The ones that do are `[starttime:sampling:endtime]`, where `starttime` is the frame at which the embryo is greater than or equal to `starting_num_cells` and `endtime` is `tend-1` (because `loadcells_unnamed.m` leaves an empty entry at the end of `embdat_stabilized` which cannot have features computed on it). `sampling` is also specified by the user. If you want to adjust which frames from `embdat_stabilized` have feature representations computed for them, adjust these values.

Next comes the feature computation. These are both done by `computeFeatures4.m`.

`computeFeatures4.m` is able to compute the following features:

| Description of feature | Formal name in pipeline | Normalization |
| ---------------------------- | ----------------------------- | ----------------- |
| First derivative of cell density in AP axis† | `first_deriv_x` | Z-score |
| First derivative of cell density in DV axis† | `first_deriv_y` | Z-score |
| First derivative of cell density in LR axis† | `first_deriv_z` | Z-score |
| Laplacian of cell density (2nd derivative of cell density) | `laplacian` | Z-score |
| Local normal vector (unit vector) | `local_normal_vector` | Divide by average magnitude |
| Residual | residual | Nothing |
| Distances to nearest neighbors | `distances_to_k_nns` | Z-score\* |
| Real space vectors to nearest neighbors | `signed_distances` | Divide by average magnitude |
| Time of beginning of current cell cycle | `beginning_of_current_cell_cycle | Scale to [0,1]\*\* |
| Time of end of current cell cycle | `end_of_current_cell_cycle` | Scale to [0,1]\*\* |
| Duration of current cell cycle | `duration_of_current_cell_cycle` | Scale to [0,1]\*\* |
| Time until next division | `time_until_next_division` | Scale to [0,1]\*\* |
| Time since last division | `time_since_last_division` | Scale to [0,1]\*\* |
| Fraction of current cell cycle elapsed | `fraction_of_current_cycle_elapsed` | Nothing |
| Real space vector to position x minutes ago | `vec_to_position_x_mins_ago` ||
| Real space vector to position x minutes in future | `vec_to_position_x_mins_in_future` ||
| Real space vector to site of next division | `vec_to_site_of_last_division` | Divide by average magnitude |
| Real space vector to site of last division | `vec_to_site_of_next_division` | Divide by average magnitude |
| Displacement over cell lifetime | `displacement_vec_over_cell_cycle` | Divide by average magnitude |
| Ratio of own cell cycle to mean of ancestors’ cell cycles | `similarity_to_ancestors_cell_cycle` | Nothing |
| Ratio of own cell cycle to mean of cousins’ cell cycles | `similarity_to_cousins_cell_cycle` | Nothing |
| Proportion of current nns that are descendants of past nns | `proportion_of_present_nns_that_are_descendants_of_past_nns` | Nothing |
| Nuclear size | `filtered_nuc_size` | Z-score |
| Number of cells in embryo | `number_of_cell_in_embryo` | Scale to [0,1]\*\* |
| Time through development (by frame)| `developmental_time` | Scale to [0,1]\*\*\* |

\* Z-score applied to all entries at the same time, not column by column  
\*\* Scaled such that 0 is the time at which there are 4 cells in the embryo (a common reference point between embryos) and 1 is the last time in the loaded dataset  
\*\*\* Scaled such that 0 is 0 and 1 is the number of cells in the embryo at the last time through which it is loaded into the pipeline.  
† The features directly compute derivatives with respect to the x,y,z coordinate axes of the data. These correspond respectively to the AP, DV, and LR axes only because of the alignment procedure 

The features very in terms of what kind of information they require from the embryo struct from which they are computed:

**Local distribution of cells (cell positions only)**

1. First derivative of cell density in AP axis
    - Calculated using a [Sobel kernel](https://en.wikipedia.org/wiki/Sobel_operator)
    - Corresponds to the AP axis because of alignment preprocessing step
2. First derivative of cell density in DV axis
    - Same as above
3. First derivative of cell density in LR axis
    - Same as above
4. Laplacian of cell density (2nd derivative of cell density)
    - Calculated using a difference of [Gaussian Kernel](https://en.wikipedia.org/wiki/Difference_of_Gaussians)
5. Local normal vector
    - Calculated using SVD
6. Residual
    - sum of perpendicular distances from neighboring cells to local tangent plane, as defined by normal
7. Distances to nearest neighbors
8. Real space vectors to nearest neighbors

**Cell migration (require cell tracking)**

9. Time of beginning of current cell cycle
10. Time of end of current cell cycle
11. Duration of current cell cycle
12. Time until next division
13. Time since last division
14. Fraction of current cell cycle elapsed
15. Real space vector to position x minutes ago
16. Real space vector to position x minutes in future
17. Real space vector to site of next division
18. Real space vector to site of last division
19. Real space vector to site of last division

**Migration within cell’s lineage (require tracking and mother/daughter relationships)**

20. Ratio of own cell cycle to mean of ancestors’ cell cycles
21. Ratio of own cell cycle to mean of cousins’ cell cycles
22. Proportion of current nns that are descendants of past nns

**Other (inherent to dataset)**
23. Nuclear Size
24. Number of cells in embryo
25. Time through development (by frame)

The way to specify what features are computed is in the `feats` argument of `computeFeatures4.m`, which is a string array containing the names of the features (middle column of table above). The list need not be in any particular order. 

### Hyperparameters

The way the features are computed is affected by a set of eight hyperparameters. Those parameters are passed in to `computeFeatures4.m`’s `params_to_use` argument, which is formatted as a cell array whose columns contains the names of the hyperparameters and the values assigned to them. They need not be in any particular order. The hyperparameters are:

| Parameter | Formal name in pipeline | Features it affects |
|-------------------|-------------------------|------------------------|
| *radius* of bounding sphere used to calculate cell density features (microns)| radius_boundingsphere_celldensity | 1, 2, 3, 4 |
| *radius* of bounding sphere used to calculate residual (microns) | radius_boundingsphere_residual | 5, 6 |
| *radius* of bounding sphere used to filter nuclear size (microns) | radius_boundingsphere_nucsize | 23 |
| *number* of neighbors for inter-cell distances (`int`) | k_intercelldistance | 7, 8 |
| *number* of generations (degree of cousins) back to look for cell cycle comparison features (int) | degree_cousins | 20, 21 |
| *time window* for migration paths (`float`)\*| timewindow_migrationpaths |14, 15 |
| *time window* for proportion of current nns that are descendants of past nns (`float`)\* | timewindow_migrationcongruency | 22 |
| *number of neighbors* for proportion of current nns that are descendants of past nns (`int`) | k_migrationcongruency | 22 |

\*The exact amount an embryo takes to undergo development varies between embryos, so to account for this, amounts of time are input as a fraction of the amount of time it takes the embryo to go from 4 to twitch (~400 mins). Therefore, a value of `timewindow_migrationpaths` = .0625 represents, in general, ~25 minutes of real time. 

### Normalization

After the specified features are computed, `computeFeatures4.m` also normalizes them according to the scheme above. For more detail about the logic behind the normalization, see the project writeup. To control or change what normalizations are performed for which features, adjust the `ordered_feats` and `ordered_norms` variables defined within `computeFeatures4.m`. These are ordered lists whose indices relate the different possible features with normalization schemes.

### Dimensionality reduction

The normalized feature space is then run through one of two dimensionality reduction algorithms: t-SNE or UMAP. The algorithm it uses is specified by the `mode` variable defined in `celeganstest.m`. The variable `neighborparam` is fed into the dimensionality reduction, specifying either the perplexity in t-SNE or `n_neighbors` in UMAP.


### What is the embeddingStruct?

Throughout the process, the pipeline creates another in-house data type--the `embeddingStruct`--primarily through the function `embeddingStructMaker.m`. The `embeddingStruct` is a struct array whose entries correspond to the rows in the feature space and which contain fields for the name of the cell that was imaged, its terminal tissue fate (if any), its top-level lineage, its real space coordinate from `embdat_stabilized`, its coordinate in the reduced space, the frame in `embdat_stabilized` from which it was obtained, and the number of cells that existed in the embryo at that time. The `embeddingStruct` makes it easy to search for time-cells in the feature space with specific properties and to color-code visualizations 

The `embeddingStruct` identifies the terminal fates specified by the variable `celltypes` and the lineages specified by the variable `lineages`, which are both string arrays. If you want to change what fates and lineages are identified and cataloged, change these variables and some code in `embeddingStructMaker.m`. 

### Visualization

With a reduced dimensionality representation of the feature space computed, the last task in `celeganstest.m` is to create colored scatter plots. The script uses the `embeddingStruct` to create arrays of RGB triplets to color the points in the low-d representation according to terminal fate, lineage, and time through development. The colors of the fates and lineages are specified by the variables `celltypecategories` and `lineagecategories`, which are arrays of RGB triplets whose entries correspond to those in `celltypes` and `lineages`. If you want to change the colors associated with these designations, change these values. Time is expressed as a linear scaling of the progression through development in hsv color space, going from green to red.

### Summary

`celeganstest.m` is the main script of this project, and it converts an embryo .zip file to a set of color-coded scatter plots in low dimensions. It uses `loadcells_unnamed.m` to load the .zip file, `coalignNamedEmbryosPerTime2.m` to align an input C. elegans embryo to a reference, `computeFeatures4.m` to compute and normalize a feature representation of the embryo’s development. In the process, it creates three types of information: an "embryo struct", a feature space representing the development, and an "embeddingStruct". 

---

## How to change the pipeline to load mutant C. elegans embryos

As stated in the Basic Pipeline Overview, loading embryos is handled by `loadcells_unnamed.m`.  The important arguments for the function are `endtime`, which is often passed in through the `celeganstest.m` variable `tend`, and the anisotropy of the embryo to load. `loadcells_unnamed.m` loads the embryo from the first frame in the raw data through the `endtime`. It outputs an embryo struct. See `celeganstest.m` and `loadWTEmbs.m` for examples of loading wild type C. elegans embryos.

In order to load a mutant C. elegans embryo, one has to simply change the `endtime` and `anisotropy` to match the new dataset. Picking the `endtime` usually has to do with identifying when the embryo begins moving, at which point tracking is unreliable. 

One caveat with mutants has to do with alignment. `coalignNamedEmbryosPerTime2.m` aligns an input embryo to a reference embryo, which is a wild type. Examine the result of aligning a mutant embryo to a wild type to using `showEmbMovie.m` to make sure the cell movements look smooth. A problem could arise during alignment if the mutant embryo differs drastically from the WT. 

See `loadMutants.m` for examples of loading mutant embryos. In that function, some of the `tend` values provided were taken from excel spreadsheets showing throughout which point an embryo was edited, which are in the project directory.

---

## How to change the pipeline to load fish embryos

The process is largely the same as for mutant C. elegans embryos. However, `coalignNamedEmbryosPerTime2.m` should be bypassed--because an input embryo to a reference C. elegans embryo using handpicked landmark cells, there is no meaningful way to apply it to a fish embryo in its current state.

For the fish, issues with performance during feature computation can arise because there are many more cells in the fish embryo. For this reason, it is recommended to increase the sampling from 2 (standard for C. elegans embryos) to, say, 10 or 20. 

Since fish are unlikely to have tracking data, one should compute their feature representations using features that require position information only. For a list of these and directions about how to do this, see the feature computation section of the Basic Pipeline Overview.

---

## How to do clustering

This section describes how to perform a 1- and 2- stage clustering on a feature space using the Louvain clustering algorithm from the Seurat package, how to manipulate them, and how to analyze their purity. 

This process uses scripts in both Matlab and R. The basic flow of information is that the feature space is computed in Matlab (`computeFeatures4.m`), then clustered using the Seurat package in `louvain_clustering.R`, then clustered again in `iterative_clustering.R` (if desired), then analyzed in `clusterEval.m`.

### How to perform one round of clustering on a feature space

The clustering is performed in `louvain_clustering.R`, which loads the feature space (a matrix saved as a .txt file) according to the specified `data_name`. The script then computes the clustering for a set of values of the Louvain ‘resolution’ parameter, specified by the variable resolutions. Small values of the resolution parameter (~1) lead fewer clusters being found, while larger values lead to more. 

The script then saves the result as a .txt file. The table contains the cluster assignments of all the time-cells in the feature space. Note: the clusters are zero-indexed. The table’s rows correspond to the time-cells and its columns correspond to distinct clustering trials specified by the different values of `resolutions`. Its column headers are named in the form: “resolution : # of clusters produced”.

### How to perform a second round of clustering

From there, one can choose to do a second round of clustering, which is done in `iterative_clustering.R`. The script loads in the first round of cluster assignments from the file saved by `louvain_clustering.R` (specified by `initial_clustering_name`). The script re-clusters a specific set of clusters (`clusters_to_recluster`) from a specific trial (`trial_of_interest`) in the initial clustering. For simplicity, it does so for one value of the resolution parameter (specified by `res`). 

The script then saves as .txt files the cluster assignments for each reclustered cluster as separate text files, whose names are of the form: “result_family_name” +  “_trialXclusterY”, where result_family_name is an input parameter, X is trial_of_interest, and Y is the number of the cluster from that trial that was reclustered. In these results, the rows correspond to the time-cells from the cluster that was reclustered, and they are in the same order as they are in the entire dataset (feature vector or initial clustering result). Note: this second set of saved files does not by itself encode what time-cells in the original dataset the clusters correspond to; X, Y in their names point the user back to the relevant entries in the first round of clustering.

### Analyzing and manipulating clustering in Matlab

Analysis of the clusters is performed in `clusterEval.m`. It loads in a first pass of clustering as a struct called `result_struct`, whose `data` field contains the matrix of cluster assignments for various values of the resolution. The rows of this matrix correspond to those of the feature space that gave rise to it, and therefore also to those in the corresponding `embeddingStruct`. This makes it easy to assess the identities of a given cluster with a few lines of code. 

The script is able to compute the purity of a cluster from the first pass directly from this matrix. By default, it counts how many clusters are: 1) ‘early’, meaning that their median timepoint is before the time at which the embryo from which the feature vector was computed was 200 cells; 2) ‘pure in sublineage’, meaning that it is more than a certain threshold (`threshold`) composed of cells of a single sublineage at the level ABxxxx, MSxx, Cxx, Dx, E, Z; ‘pure in fate’, meaning that it is more than a certain threshold (`threshold`) composed of cells whose sublineages (above) have the same mode fate among their terminal children; or impure, meaning it does not fit into any other category. Note: the ‘pure in fate’ metric was based off of the sublineage-terminal fate correspondences found in [this resource]( https://www.researchgate.net/figure/Overview-of-Caenorhabditis-elegans-cell-fate-map-and-division-asynchrony-A-A-Nomarski_fig2_278042073) and was used because of time-consuming technical issues arising from using the proper embryo parts list to compute exact terminal fates. To change how purity metrics are calculated, adjust the `clusterPurity2.m` function at the bottom of the script. 

The purity metrics can also be calculated on clusters produced by the second round of clustering. This is done by way of a `reclustering_struct`. The `reclustering_struct` is a struct array that connects the information in the cluster assignments from the first round and second round. It loads each file in a family of re-clustering results produced by `iterative_clustering.R`. Its entries correspond to the re-clustered clusters. For each one the field `new_clusters` gives the cluster assignments from the second round of clustering and the field `old_clustering_inds` contains the indices that those cells correspond to in the entire feature space. This object makes it more convenient to manipulate a second round of clustering in Matlab. It is used to organize the purity calculation for re-clustered clusters.

For concreteness about how to do this, consider a small clustering example. Say you compute and save a feature space under the name `featurespace.txt`. You then load this in `louvain_clustering.R`, and cluster it for resolution = 1.2 and 1.5. If there are 20 cells in the feature space, then the results will be a 20 x 2 matrix containing the cluster assignments for the two trials. You then save this as `first_clustering.txt`. Then, you open that in `iterative_clustering.R`. Say you decide to re-cluster clusters 6, 7, and 11 from the trial 1 (resolution = 1.2) of the first clustering. This produces three .txt files called, say, `second_clustering_trial1cluster6.txt`, `second_clustering_trial1cluster7.txt`, and `second_clustering_trial1cluster11.txt`.  You can then open these in succession in `clusterEval.m`. These will make a `reclustering_struct` whose elements 1, 2, and 3 correspond to these three re-clusterings. Their `new_clusters` fields give the cluster assignments for the time-cells in those clusters, and their `old_clustering_inds` fields give the indices of those time-cells as they appear in the parent feature space.

---

## Function Documentation

### Key Functions 

`celeganstest.m`  
Parent function of the project. Takes in an embryo .zip file and produces a set of color-coded scatter plots, a feature representation of the development, and an `embeddingStruct` for convenience and organizational purposes.

Arguments:  
- `emb`: Full pathname of acetree .zip file to be opened/analyzed
- `input_xy_res`: x and y resolution at which the embryo was imaged; available in corresponding .xml file
- `anisotropy`: Anisotropy at which embryo was imaged; available in corresponding .xml file
- `tend`: Last time (frame) from which to load cell position data from .zip file
- `centering_time`: Time (frame) used to align centroids of input C. elegnans embryo and reference embryo; by default set to tend
- `starting_num_cells`: Defines time (frame) from which to start *computing features*; distinct from loading
- `sampling`: Defines frequency with which frames from loaded embryo (embdat) have their feature representation computed. Ex: sampling = 2 means every other frame will have feat. rep. computed
- `feats`: List of features you want to compute. Order doesn't matter
- `params_to_use`: Pairs of hyperparameter names and desired values. Order doesn't matter. If computing a subset of all features and don't use a certain hyperparameter, just assign it a dummy value
- `normalize_bool`: Whether computeFeature4 should normalize the feature representations.  Normalizations corresponding to the input features are defined within the function itself
- `mode`: Either "tsne" or "umap"; whether to do dimensionality reduction with t-SNE or UMAP. Note: UMAP requires that umapFileExchange (1.4.2 or later) be in the path
- `neighborparam`: Primary parameter of dimensionality reduction; perplexity in t-SNE and n_neighbors in UMAP
- `only2d`: Whether to do dimensionality and reduction and visualization to 2 dimensions, as opposed to do it to both 2 and 3
- `celltypes`: List of celltypes with which to annotate data in embeddingStruct and color points in visualization
- `celltypecategories`: RGB colors assigned to celltypes. Must be in same order as celltypes
- `lineages`: List of lineages with which to annotate data in embeddingStruct and color points in visualization
- `lineagecategories`: RGB colors assigned to lineages. Must be in same order as lineages

Outputs:  
- `embdat_stabilized`: the embryo struct containing the raw information from the embryo .zip file. Each entry is a frame from the imaging data and is a struct containing the (aligned) positions of the cells therein in real space, their names, their succession, and other information; computed by `loadcells_unnames.m`. 
- `normed_percellfeature`: the normalized feature space; computed by `computeFeatures4.m`. 
- `embeddingStruct`: a struct array whose indices correspond to the entries in the `normed_percellfeature` and which contain each entry’s cellID, terminal fate (if any), top level lineage, real space coordinate, reduced dimension coordinate, frame from `embdat_stabilized` from which it was taken, and the number of cells there were in the embryo at that time; computed by `embeddingStructMaker.m`. 

`loadcells_unnamed.m`  
Creates an embryo struct from an unzipped embryo .zip file. Written by AS prior to June 2021. 

`parseCellsFromEmb.m`  
Creates a cells object containing information about each cell ID in a given embryo struct. Not utilized in this pipeline, but helpful for reference. 


`coalignNamedEmbryosPerTime2.m`  
Aligns a C. elegans embryo to a reference by computing the best fit rotation between a user-defined set of landmark cells defined in both embryos. Used as a key preprocessing step. This function is where one could devise a method to align an input embryo through a particular time through development and leave the remainder of the dataset unaligned.

Arguments:  
- `embref`: `embryoStruct` of reference embryo
- `tstartref`: start time of reference emb, usually just set to 1
- `tendref`: last timepoint of reference emb
- `xyscaleref`: xy resolution of reference emb
- `emb`: embryo struct of embryo being aligned
- `tstart`: start time of emb being aligned
- `tend`: last timepoint of emb being aligned
- `anisotropy`: anisotropy of emb being aligned
- `xyscale`: xy resolution of emb being aligned
- `centering_time`: timepoint in emb being aligned used to compute the translation to center it at the origin
Outputs:  
- `transformedemb`: the aligned embryo

`computeFeatures4.m`  
Compute and normalize spatiotemporal feature vector for each cell at each time in an embyro struct. Use this to control what features are computed and how/whether they are normalized.

Arguments:  
- `embryo`: the embryo struct from which to compute features
- `times`: the times (frames) from embryo struct at which cells will have their feature representations computed
- feats`: contains the names of the features you wish to compute; can be in any order
- `param_values`: contains the names (column 1) and values (column 2) of the hyperparameters required to compute features; can be in any order
- `normalize_bool`: whether or not to normalize
- `input_xy_res`: xy resolution of the embryo
Outputs:  
- `featurevector`: matrix in which rows are feature representation of each individual cell, unnormalized
- `normed_featurevector`: same as `featurevector` but normalized according to the scheme below. Features as listed in `ordered_feats` correspond to the normalizations in `ordered_norms`
- `indices_list`: array whose columns contain the starting and ending indices of each feature, respectively. ex. if feature 3 goes from columns 6 to 8 in the `featurevector`, then `indices_list[3,:] = [6,8]`

`embeddingStructMaker.m`  
Used to create `embeddingStruct`. Note: in order to assign terminal fate tissues, requires that `allpharynx.mat`, `allneuron.mat`, `allhypoderm.mat`, `allbodywallmuscle.mat`, `allgut.mat` to be in the path. 

`louvain_clustering.R`  
Performs Louvain clustering on an embryo feature space and saves the result as a .txt file.

Arguments:  
- `data_directory`: The directory to load from and save to
- `data_name`: the name of the .txt file containing the feature space
- `column_name`: the name of the .txt file containing the ordered names of the columns (features) in the feature space
- row_name`: the name of the .txt file containing the ordered names of the rows (cells) in the feature space
- `result_name`: the name of the .txt file that will contain the clustering result
- `resolutions`: an array containing the values of the resolution (key parameter of Louvain algorithm) you want to use

Outputs:  
- `results`: a matrix containing the cluster assignments of all time-cells in the feature space. Rows correspond to the time-cell entries in the uploaded feature space; columns correspond to each clustering done with a different resolution. Headers are of format “resolution:number of clusters it produced”

`iterative_clustering.R`  
Performs the second step of an iterative clustering process. Loads a first clustering pass
as a text file, performs a second round of clustering on specified clusters from a specified trial within it, then saves the clustering results from those to individual text files.

Arguments  
- `data_directory`: The directory to load from and save to
- `data_name`: the name of the .txt file containng the feature space
- `column_name`: the name of the .txt file containing the ordered names of the columns (features) in the feature space
- `initial_clustering_name`: the name of the .txt file containing the results from the first round of clustering
- `result_family_name`: the name you want all cluster asignment output files to start with. names are. “result_family_name_trialXclusterY”, where X is `trial_of_interest` and Y is the cluster from that that was re-clustered.
- `trial_of_interes`t: the trial from the initial clustering that you want to run the second iterative step on
- `clusters_to_recluster`: the indices of the clusters from the `trial_of_interest` that you want to apply the second clustering step to; zero-indexed  
- `res`: the resolution parameter for Louvain clustering
Outputs  
- `results`: a vector containing the cluster assignments of the time-cells in the cluster from the initial clustering that was re-clustered. time-cells (rows) are in the same order as they are in the feature space/initial clustering. Headers are of format “resolution:number of clusters it produced”

`clusterEval.m`  
Opens and manipulates cluster assignment data from R Louvain clustering scripts and analyzes their purity. Opens .txt files containing clustering assignments from `louvain_clustering.R` and `iterative_clustering.R`. For a detailed discussion about how to use the script, see the ‘How to do clustering’ section. 


### Other helpful functions

`loadMutants.m`  
Loads and saves mutant embryo .zip files. Computes embryo structs, features spaces, and embeddingStructs. Works identically to the beginning of `celeganstest.m`.

`fateAndLineageScores.m`  
Computes ‘terminal fate clustering’ and ‘lineage continuity’ scores for a feature space. ‘Terminal fate clustering’ is percent of k nns that are of same tissue type, averaged over all terminal (annotated to be of a terminal fate) cells in the embryo. ‘Lineage continuity’ is percent of k nns that are self, mother or daughter, averaged over all cells in the embryo. In both cases, information about tissue types and cell ID come from `embeddingStruct`, which is input as an argument.

`hyperparameterTuning4.m`  
Runs hyperparameter optimization. Creates a cell array called `optimization_cases`, which contains information about each trial (parameter combination) used in the optimization. Then computes and scores the feature space associated with each one, and saves the result to a new cell array called `optimization_cases_WITHRESULTS`, which contains the scores and feature spaces associated with each trial. Computes the same metrics as `fateAndLineageScores.m`, but with some modifications. 1) the terminal fate clustering score is computed in a feature space where there is only one entry per cell ID, not multiple. That entry is found by taking the entry for a particular cell whose capture frame (from the embryo struct) is the median of all frames for that cell in the dataset. 2) the lineage continuity score only looks at mother and daughter cells within a certain amount of time forward and back from the cell being queried. In all past optimizations, was 7.5 minutes. 

`knnClassifier.m`  
Runs a k-nn classification using three WT C. elegans embryos. Classifies cells by tissue type and cell ID. Computes over all training and validation pairs of the three embryos. Allows user to reduce space to lower dimensions (d <= dimension of feature space) before scoring. Has variables to define the sweeps over d and k, the number of neighbors used in the classification. 

`loadWTEmbs.m`  
Loads and saves diSPIM WT embs, same to `loadMutants.m`. Works identically to the beginning of `celeganstest.m`.

`overlayEmbryos.m`  
Creates a UMAP plot showing a color-coded UMAP reduction of the feature spaces of two embryos. Feature spaces are concatenated and then reduced simultaneously by UMAP. Used for comparisons of two feature spaces. 

`showEmbMovie.m`  
Shows a movie of an embryo developing over time. Convenient for quick visualization. 


#### Exploratory notebooks and scratchwork

`analyzeClusters.m`  
Used to analyze Louvain clustering early fall 2021. Pie charts of clusters by coarse-level sublineages, movies colored by clusters, investigation of terminal fate composition of clusters using parts list. 

`backpropagateLineages.m`  
Used to create lines connecting mother/daughter cells in UMAP visualizations

`graphExperimentation.m`  
`normalExperimentation.m`  
Used to explore normal orientation propagation while designing cohereNormals.m.

`pseudoManifoldTest.m`  
Used to explore how dimensionality reduction affected performance metrics in hyperparameter optimization. Results became part of implementation of `hyperparameterTuning4.m`.

### Aux and helper

`avgMagnitude.m`  
Normalizes an array of real space vectors by their average magnitude. Used in feature computation.

`backTraverse.m`  
Given the frame and index of a cell in an embryo struct, gives the frame and index of an ancestor cell in its lineage a specified amount of time or generations in the past. Used in feature computation. 

Arguments:  
- `start_time`: the frame of the starting cell
- `start_ind`: the index of the starting cell (from the names field of that entry in the embdat struct)
- `mode`: (‘timesteps or generations’)specifies whether time_parameter refers to a number of frames or a number of generations in the cell lineage
- `time_parameter`: the number of timesteps (frames) or generations you go in the past
- `embryo`: the embryo struct (often called `embdat`) being analyzed 
Outputs:   
- `time_pt`: time frame of the inputs cell’s ancestor
- `ind_pt`: index frame of the inputs cell’s ancestor (from the names field of that entry in the `embdat` struct)
- `ancestor_times`/`ancestor_inds`: ordered lists containing the timepoints/indices of all new cells IDs encountered while tracing backward in the lineage. Their indices correspond to each other, i.e. `ancestor_times[2]` corresponds `ancestor_inds[2]`.

`backTraverse_keepall.m`
Same as `backtraverse.m` except that `ancestor_times`/`ancestor_inds` keep the timepoints/indices of every single timecell encountered while traversing back through the lineage, not the timepoints/inds of each new cell ID encountered. Was created for exploratory purposes in order to avoid having to change all calls to `backtraverse.m` throughout the pipeline

`cellCycleMetrics.m`  
Calculates various metrics associated with a given cell’s cell cycle. Used in feature computation.

Arguments: 
- `t`: timepoint of cell in question in embryo struct
- `i`: index of cell in question in embryo struct
- `embryo`: embryo struct being analyzed
Outputs:  
- `beginning_of_current_cycle`: timepoint of beginning of cell’s cell cycle
- `end_of_current_cycle`: timepoint of end of cell’s cell cycle
- `duration_of_current_cell_cycle`: duration of cell’s cell cycle
- `time_until_next_division`: time remaining in cell’s cell cycle
- `time_since_last_division`: time elapsed since previous division
- `fraction_of_current_cycle_elapsed`: `time_since_last_division`/ `duration_of_current_cell_cycle`

`cohereNormals.m`  
Points normal vectors on a surface whose directions are correct whose signs are not in a coherent (outward) direction using algorithm in [this resource](http://hhoppe.com/recon.pdf) Used in feature computation. 

`compareTimesteps.m`  
Plots the cell positions at two timepoints from the same embryo struct in the same scatter plot. 

`embeddingLookup.m`  
Isolates the entries of an `embeddingStruct` whose points are within some user-specified range in the embedded space. Used for quickly seeing which time-cells make up a particular group in dimensionality reduced space.

`featSummary.m`  
Plots raw values of easily visualizable features from a feature space. The columns of the features in question are specified by the user.

`forwardTraverse.m`  
Given the frame and index of a cell (or cells) in an embryo struct, gives the frame and indices of its descendants a specified amount of time or number of generations in the future. Used in feature computation.

Arguments: 
- `start_time`: timepoint (frame) of starting cell(s)
- `cells_to_look_at`: indices of starting cell(s)
- `mode`: ‘timesteps’ or ‘generations’; specifies whether `time_parameter` refers to an amount of time (number of frames in the embryo) or number of generations
- `time_parameter`: either the number of timesteps or generations to advance by
- embryo: the embryo struct

Outputs:  
- `time_pt`: the timepoint (frame) of descendant cells
- `cells_to_look_at`: the indices of the descendant cells

`forwardTraverse_keepall.m`  
Same as `forwardTraverse.m`, but keeps two ordered lists of the timepoints and indices of all descendant cells encountered while moving forward through the lineage, analogously to `backtraverse_keepall.m`.

Arguments:  
- same as `forwardTraverse.m`  
Outputs:  
- same as `forwardTraverse.m`
- `all_later_cell_times`: the times (frames) of all descendant cells encountered 
- `all_later_cell_inds`: the indices of all descendant cells encountered

`gaussians.m`  
Calculates the values of a origin-centered 3-dimensional normal distribution, its first derivatives with respect to the three coordinate axes, or the Laplacian at some point in space. Used in feature computation.

`internallyAlignNamedEmbryo.m`  
Removes an embryo’s rotations in the eggshell during development using affine transformations.

`loadEmbryo_unzipped.m`  
Helper function called by` loadcells_unnamed.m` involved in creating an embryo struct from an unzipped embryo .zip file.

`normalize3.m`  
Normalizes a `featurevector`. Functionality is same as that in `computeFeatures4`, but this function can be used in other exploratory scripts more conveniently.

Arguments:  
- `percellfeature`: the feature space
- `feature_names_list`: a string array of the different features in the feature space
- `sizes_list`: a # of features x 2 array where the columns give the starting and ending columns of each feature, respectively. Ex: if the second feature occupies columns 5 through 8 in the feature space, then `sizes_list[2,:] = [5,8]`. Rows correspond to rows in `feature_names_list`.
- `norm_types_list`: normalizations associated with each feature, chosen from ‘scale_4cell_to_end’, ‘scale_0_to_max_num_cells’, ‘column_zscore’, ‘all_zscore’, ‘divide_avg_mag’. Rows correspond to rows in `feature_names_list`.
- `embryo`: the embryo struct.
Outputs:  
- `P`: the normalized feature space

`normalsPlot.m`  
Plots normal vectors and points they originate from in one 3d scatter. Used while designing `cohereNormals.m`.

`stepback.m`  
Given a cell’s timepoint and index in an embryo struct, gives the index of the cell in the previous timepoint that gives rise to the cell in question. Used as a helper function in `backTraverse.m`.
