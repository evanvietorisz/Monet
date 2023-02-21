"

Louvain Clustering

script to perform Louvain clustering on an embryo feature space and save the result as a .txt file

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

result_name
  the name of the .txt file that will contain the clustering result

resolutions (array)
  an array containing the values of the resolution (key parameter of Louvain algorithm) you want to use
  
  
OUTPUTS
_______

results (array of shape # of rows in feature space x length(resolutions))
  a matrix containing the cluster assignments of all time-cells in the feature space.
  Rows correspond to the time-cell entries in the uploaded feature space;
  columns correspond to each clustering done with a different resolution.
  Headers are of format  'resolution:number of clusters it produced'
  
"

  



#Backend

#set working directory
working_directory <- "/Users/evanvietorisz/Documents/R"
setwd(working_directory)

#load packages
library("Seurat")
library("ggplot2")


#Load data

data_directory <- "/Users/evanvietorisz/Documents/R/Data"
data_name <- "dispim_emb2_full_featurespace_Nov6.txt"
column_name <- "dispim_emb2_columnnames_Nov6.txt"
row_name <- "dispim_emb2_rownames_Nov6.txt"
#note: for plain text, last character of data has to be EOL (\n or return)

#set name of result
result_name <- "clustering_result.txt"

#set resolutions on which to do clustering
resolutions <- seq(1, 3, by = 0.1)



pth <- file.path(data_directory,data_name)
cell_data <- read.table(file=pth)

pth <- file.path(data_directory,column_name)
columnnames <- read.table(file=pth,sep = "\t") #necessary to specify bc entries have spaces

pth <- file.path(data_directory,row_name)
rownames <- read.table(file=pth)


#Seurat 

#create Seurat object
S_obj.data <- t(cell_data)
#note: transpose bc Seurat uses rows as features, cols as cells

colnames(S_obj.data) <- rownames(cell_data) #to make the next command work

S_obj <- CreateSeuratObject(counts = S_obj.data, project = "louvain")

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
results <- matrix(nrow=nrow(cell_data),ncol=length(resolutions))
trial_descriptions <- character(length(resolutions))
for (trial in 1:length(resolutions)) {
  
  res <- resolutions[trial]
  
  S_obj <- FindClusters(S_obj, resolution = res)
  cluster_ids <- as.vector(S_obj@meta.data[["seurat_clusters"]]) #note: is a factor
  cluster_ids_int <- as.numeric(as.character(cluster_ids))#convert to ints
  #note: "seurat_clusters gets overridden for each trial"
  
  results[,trial] <- cluster_ids_int
  trial_descriptions[trial] <- paste(c(as.character(res),":",as.character(max(cluster_ids_int)+1)),collapse="") #note: +1 bc zero-indexed

}

#save clustering output

pth <- paste(c(data_directory,result_name),collapse="")
write.table(results, file = pth, append = FALSE, sep = "\t",
            row.names = FALSE, col.names = trial_descriptions)

"
#note: how to open
pth <- file.path(data_directory,result_name)
open_file <- read.table(file=pth)

opened_data <- open_file[2:nrow(open_file),]
opened_columnnames <- as.character(open_file[1,])
"




