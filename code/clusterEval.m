%script to open/manipulate clustering data from R scripts and analyze
%clusters' purity


%prep

%raw data, parameters

%holds initial clustering result
result_struct = importdata('bulk_clustering_jan8.txt');%raw data from R; columns are different experiments "resolution:number of clusters"

%Yichi's old clustering
%{
%result_struct = importdata('Nov6_louvain_dim55_clus.txt');%old = louvain_dim55_clus.txt; new = Nov6_louvain_dim55_clus.txt;
%}

%%
%embryo struct whose feature space yielded the clustering
emb = load('dispim_emb2_full.mat'); %note: dispim_emb2.mat is through comma; dispim_emb2_full.mat is through 1.5. dispim_emb2 is WG emb
featurespace = emb.normed_percellfeature;
embStruct = emb.embeddingStruct;

%%
%aux

%compile the names
all_names = cell(1);
for i = 1:length(emb.embeddingStruct)
    
    all_names{i} = char(emb.embeddingStruct(i).name);

end
all_names = all_names';

%compile all the celltypes
all_types = cell(1);
for i = 1:length(emb.embeddingStruct)
    
    all_types{i} = char(emb.embeddingStruct(i).celltype);

end
all_types = all_types';

%taken from https://www.researchgate.net/figure/Overview-of-Caenorhabditis-elegans-cell-fate-map-and-division-asynchrony-A-A-Nomarski_fig2_278042073
sublins_and_modefates = ["ABalaa","neuron"; %for coloring the text labels
              "ABalap","neuron";
              "ABalpa","pharynx";
              "ABalpp","neuron";
              "ABaraa","pharynx";
              "ABarap","pharynx";%split w neuron
              "ABarpa","hypoderm";%split w neuron
              "ABarpp","hypoderm";
              "ABplaa","hypoderm";%split w neuron
              "ABplap","hypoderm";%split w neuron
              "ABplpa","neuron";%split w excretory
              "ABplpp","neuron";
              "ABpraa","neuron";%split w hypoderm
              "ABprap","neuron";%split w hypoderm
              "ABprpa","neuron";
              "ABprpp","neuron";
              "MSaa","pharynx";
              "MSap","bodywallmuscle";
              "MSpa","pharynx";
              "MSpp","bodywallmuscle";
              "Caa","hypoderm";
              "Cap","bodywallmuscle";
              "Cpa","hypoderm";
              "Cpp","bodywallmuscle";
              "Da","bodywallmuscle";
              "Dp","bodywallmuscle";
              "E","gut";
              "Z","germline";
              "early/other","none"];
          
sublins = sublins_and_modefates(:,1);
modefates = sublins_and_modefates(:,2);

%construct master list [cell name, sublineage, associated modefate]
names_sublins_modefates = repmat(["seed"],length(all_names),3);
for i = 1:length(all_names)
    
    names_sublins_modefates(i,1) = all_names{i};  
        
    
    for j = 1:length(sublins)-1 %-1 bc will never actually equal the last name "early/other"
        
        if startsWith(all_names(i),sublins(j))
            
            names_sublins_modefates(i,2) = sublins(j);
            names_sublins_modefates(i,3) = modefates(j);

            break
            
        end
        
        %if not filled, mark as "early/other"
        if names_sublins_modefates(i,2) == "seed"
            
            names_sublins_modefates(i,2) = "early/other";
            names_sublins_modefates(i,3) = "none";
            
        end
        
        
    end


end

%%

%code to evaulate a single pass of clustering

%computes purity directly from result_struct

dat = result_struct.data;

%original clustering

trial_num = 11;
threshold = .75;

purity_data = [];
for cluster_num = 0:max(dat(:,trial_num))
    
    inds = find(dat(:,trial_num)==cluster_num);
    
    %with early category
    purity = clusterPurity2(inds,threshold,names_sublins_modefates,embStruct,early_threshold);%names_sublins_modefates defined earlier in this script

    purity_data = [purity_data;purity];
    
end


%%

%code to make reclustering_struct, used for evaluating iterative clustering

%what the reclustering_struct is
%{

reclustering_struct is struct array whose entries are the differnet
clusters that were redone, in order. built from the individual files
containing the second round clustering information

.new_clusters is ordered (relative to all_names) clustering data for the cells in that cluster
.old_clustering_inds are the indices in the original embeddingStruct that
the entries in new_clusters correspond to

%}


%key things to change when importing new data

%clusters_that_were_redone
%   refers to the clusters from the original clustering that were
%   reclustered. starts at 0 because clustering is done in R, which is
%   0-indexed

%filename_root
%   the names of the files that contain the cluster data from the
%   reclustered clusters

%settings on existing data
%{
iterative clustering 

trial 1 (recluster largest 9 clusters):

clusters_that_were_redone = [0:8]
filename_root = reclustering_of_jan8experiment_trial11cluster

trial 2 (recluster largest 49 clusters):

clusters_that_were_redone = [0:48]
filename_root = second_reclustering_of_jan8experiment_trial11cluster

%}

clusters_that_were_redone = [0:48]; %these numbers should be the same as
                                    %every number that will come after
                                    %'cluster' in the file names below
filename_root = 'second_reclustering_of_jan8experiment_trial11cluster';


%create the reclustering struct
reclustering_struct = struct();
for i = 1:length(clusters_that_were_redone)
    
    clus = clusters_that_were_redone(i);
    
    filename = [filename_root,num2str(clus),'.txt'];
    temp = importdata(filename);

    
    reclustering_struct(i).new_clusters = temp.data;
    reclustering_struct(i).old_clustering_inds = find(dat(:,11)==clusters_that_were_redone(i));
    
end

%%

%iterative clustering cont

%compute purity metrics across a set of both first stage and second stage
%clusters
purity_data = [];
for cluster_num = 0:max(dat(:,trial_num))

    %description of what the below does
    %{
    
    if the current cluster is in the list of ones that were redone
    
        figure out which entry of the result struct corresponds to that cluster
    
        for all the clusters in that .data field
            compute the purity
            append it to the master list
        end
    end
    
    %}
    if ~isempty(find(clusters_that_were_redone == cluster_num))
        
        ind_to_use = find(clusters_that_were_redone == cluster_num);
        
        new_clusters = reclustering_struct(ind_to_use).new_clusters;
        
        for j = 0:max(new_clusters)
            
            local_indices_of_cluster_j = find(new_clusters==j);
            inds = reclustering_struct(ind_to_use).old_clustering_inds(local_indices_of_cluster_j);
            
            purity = clusterPurity2(inds,threshold,names_sublins_modefates,embStruct,early_threshold);%names_sublins_modefates defined earlier in this script

            purity_data = [purity_data;purity];
            
        end
        
    else
        
        inds = find(dat(:,trial_num)==cluster_num);
       
        purity = clusterPurity2(inds,threshold,names_sublins_modefates,embStruct,early_threshold);%names_sublins_modefates defined earlier in this script

        purity_data = [purity_data;purity];

    end
    
end


%%
function [purity] = clusterPurity2(inds,threshold,names_sublins_modefates,embStruct,early_threshold)
%{ 
tells you whether a particular cluster meets purity criteria (greater than
threshold percent:

%note: early_threshold is the number of cells before which qualifies as
early

-1) early
all others not early
0) impure 
1) pure in coarse-level sublineage
2) pure in fate corresponding to coarse-level sublineage


inds is the indices of the cells in the cluster, either in 1s and 0s or the
indices themselves

threshold is the minimum percent composition of a single category a cluster must have to be considered pure
 
names_sublins_modefates is a list of all names, sublineages, corresponding
mode fates of all time-cells in the dataset

%}

%default
purity = 0;

%aux
entries = names_sublins_modefates(inds,:);


%check if early

%find last time in embryo at which it's 200 cells big
ct = 1;
while embStruct(ct).numCellsAtCapture < early_threshold + 1 
    ct = ct+1;
end
early_threshold_frame = embStruct(ct).captureFrame;

%if median capture frame from selection of cells is earlier than that, it's
%early, it's early
if median(unique([embStruct(inds).captureFrame])) < early_threshold_frame 
    purity = -1; 
end

%check if pure in coarse-level sublineage (all pure pure in sublineage will be pure in corresponding mode fate) 
mode_count = length(find(entries(:,2)==string(mode(categorical(entries(:,2))))));%number of times the mode SUBLIN entry shows up

if mode_count/size(entries,1) >= threshold
    purity = 1;
    return
end


%if not pure in sublineage, then could still be pure in mode fate
if purity ~= 1
    
    mode_count = length(find(entries(:,3)==string(mode(categorical(entries(:,3))))));%number of times the mode MODE FATE entry shows up
    
    
    if mode_count/size(entries,1) >= threshold
        purity = 2;
    return
    
end
    
end


end









