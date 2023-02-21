"

Iterative Louvain Clustering

script to perform the second step of an iterative clustering process. Loads a first clustering pass
as a text file, performs a second round of clustering on specified clusters from a specified trial within it,
then saves the clustering results from those to individual text files

PARAMETERS
_________

data_directory
  The directory to load from and save to
  
data_name
  the name of the .txt file containng the feature space
  
column_name
  the name of the .txt file containing the ordered names of the columns (features) in the feature space
  
row_name
  the name of the .txt file containing the ordered names of the rows (cells) in the feature space
  
initial_clustering_name
  the name of the .txt file containing the results from the first round of clustering
  
result_family_name
  the name you want all cluster asignment output files to start with. names are
  result_family_name + 'trialXclusterY', where X is trial_of_interest and
  Y is the cluster from that that was reclustered. 

trial_of_interest
  the trial from the initial clustering that you want to run the second iterative step on
  
clusters_to_recluster
  the indices of the clusters from the trial_of_interest that you want to apply the second clustering step to;
  zero-indexed
  
res
  the resolution for Louvain clustering

  
OUTPUTS
_______

results (array of shape # of rows in feature space x 1)
  a vector containing the cluster assignments of the time-cells in the cluster from the
  initial clustering that was reclustered. time-cells (rows) are in the same order as
  they are in the feature space/initial clustering. Headers are of format  
  'resolution:number of clusters it produced'
  
"

data_directory <- "/Users/evanvietorisz/Documents/R/Data"
data_name <- "dispim_emb2_full_featurespace_Nov6.txt"
column_name <- "dispim_emb2_columnnames_Nov6.txt"
row_name <- "dispim_emb2_rownames_Nov6.txt"
#note: for plain text, last character of data has to be EOL (\n or return)

#the first round of clustering
initial_clustering_name <- "bulk_clustering_jan8.txt"

#parent name of output files (result_family_name + "_trialXclusterY")
result_family_name <- "second_reclustering_of_jan8experiment"

#the clustering from the first round of clustering that you want to do the second clustering round on
trial_of_interest <- 11

#which clusters from that first round you actually want to recluster
clusters_to_recluster <- 0:48 #note: zero-indexed

#resolution of louvain clustering
res <- 2


#load packages
library("Seurat")
library("ggplot2")

#Backend

#set working directory
working_directory <- "/Users/evanvietorisz/Documents/R"
setwd(working_directory)

#Load data

pth <- file.path(data_directory,data_name)
initial_cell_data <- read.table(file=pth)

pth <- file.path(data_directory,column_name)
columnnames <- read.table(file=pth,sep = "\t") #necessary to specify bc entries have spaces

pth <- file.path(data_directory,row_name)
initial_rownames <- read.table(file=pth)

#initial clustering result
pth = file.path(data_directory,initial_clustering_name)
initial_clustering <- read.table(file=pth,header=TRUE)

#pare down

for (clus in clusters_to_recluster){
  
  cluster_of_interest <- clus
  
  inds <- initial_clustering[,trial_of_interest]==cluster_of_interest
  cell_data <- initial_cell_data[inds,]
  rownames <- initial_rownames[inds,]
  
  #Seurat
  
  #create Seurat object
  S_obj.data <- t(cell_data)
  #note: transpose bc Seurat uses rows as features, cols as cells
  
  colnames(S_obj.data) <- rownames(cell_data) #to make the next command work
  
  S_obj <- CreateSeuratObject(counts = S_obj.data, project = "louvain_jan3")
  
  #specify that all features are 'variable features' so that subsequent steps process all of them
  S_obj <- FindVariableFeatures(S_obj, selection.method = "vst", nfeatures = nrow(columnnames))
  
  #populate scale.data field with cell_data
  S_obj[["RNA"]]@scale.data <- t(as.matrix(cell_data))
  #note: transpose
  
  #create neighbor graph from clustering
  S_obj <- FindNeighbors(S_obj,features = VariableFeatures(object=S_obj),dims=NULL)
  #note: dims = NULL tells the function not to use PCA coords for the graph and instead just use the raw data (variable features)
  
  #Do louvain clustering
  
  #create objects for results
  result <- matrix(nrow=nrow(cell_data),1)
  
  S_obj <- FindClusters(S_obj, resolution = res)
  cluster_ids <- as.vector(S_obj@meta.data[["seurat_clusters"]]) #note: is a factor
  cluster_ids_int <- as.numeric(as.character(cluster_ids))#convert to ints
  #note: "seurat_clusters gets overridden for each trial"
  
  result <- cluster_ids_int
  trial_description <- paste(c(as.character(res),":",as.character(max(cluster_ids_int)+1)),collapse="")
  
  #save clustering output
  result_name <- paste(c(result_family_name,"_trial",as.character(trial_of_interest),"cluster",as.character(cluster_of_interest),".txt"),collapse="")
  
  pth <- paste(c(data_directory,result_name),collapse="")
  write.table(result, file = pth, append = FALSE, sep = "\t",
              row.names = FALSE, col.names = trial_description)
  
  #remove S_obj in case there could arise problems with overwriting the fields
  rm(list="S_obj")
  
}


