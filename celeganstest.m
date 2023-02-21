%{
%CELEGANSTEST
    Compute spatiotemporal feature representation of a C. elegans embryo
    and visualize it in low dimensions using UMAP/t-SNE. Parent script of
    project.

DEPENDENCIES
____________
Matlab Packages:
    Statistics and Machine Learning Toolbox
    Curve Fitting Toolbox
    Bioinformatics Toolbox

Third Party Packages:
    StarryNite-master
    umapFileExchange (1.4.2 or later)

Functions:
    Loading raw data
        loadEmbryo_unzipped.m
        loadcells_unnamed.m
        parseCellsFromEmb.m
    Embryo alignment
        coalignNamedEmbryosPerTime2.m
    Feature Computation
        computeFeatures4.m
        avgMagnitude.m; 
        backTraverse.m
        cellCycleMetrics.m
        cohereNormals.m
        forwardTraverse.m
        gaussians.m
        stepback.m
    EmbeddingStruct
        embeddingStructMaker.m

Files:
    ref_emb.mat (in directory); used for C. elegans embryo alignment in coalignNamedEmbryosPerTime2
    allpharynx.mat; used to annotate terminal tissues in embeddingStructMaker
    allneuron.mat; "" "" 
    allhypoderm.mat; "" "" 
    allbodywallmuscle.mat; "" "" 
    allgut.mat; "" "" 


IMPORTANT OUTPUTS
_______
embdat_stabilized (struct array)
    A struct where each entry is a frame loaded from the embryo.zip file;
    each frame contains fields for the (post-alignment) real space
    positions of cells, their names, their daughters, and other info.
    Closest thing to the raw data in the pipleine
    
normed_percellfeature (array, shape depends)
    An array whose rows are feature representations of each cell in each
    frame of embdat_stabilized. Size depends on which features are computed
    and which hyperparameters are used. The main idea of the project
    
embeddingStruct (struct array)
    A struct array where the entries are the same as normed_percellfeature
    and contain each time-cell's name, fate, lineage, real space
    coordinate, frame it came from in embdat_stabilized, and other info.
    Used for easily creating scatter plot color schemes, searching for
    time-cells cells with certain properties
    


PARAMETERS
__________
emb (char vector)
    Full pathname of acetree .zip file to be opened/analyzed

input_xy_res (scalar)
    x and y resolution at which the embryo was imaged; available in
    corresponding .xml file

anisotropy (scalar)
    Anisotropy at which embryo was imaged; available in
    corresponding .xml file
    
tend (scalar)
    Last time (frame) from which to load cell position data from .zip file

centering_time (scalar)
    Time (frame) used to align centroids of input C. elegnans embryo and
    reference embryo; by default set to tend
    
starting_num_cells (scalar)
    Defines time (frame) from which to start *computing features*; distinct from loading
 
sampling (scalar)
    Defines frequency with which frames from loaded embryo (embdat) have
    their feature representation computed. Ex: sampling = 2 means every
    other frame will have feat. rep. computed

feats (string array)
    List of features you want to compute. Order doesn't matter

params_to_use (cell array of shape 8 x 2)
    Pairs of hyperparameter names and desired values. Order doesn't matter.
    If computing a subset of all features and don't use a certain
    hyperparameter, just assign it a dummy value

normalize_bool (bool)
    Whether computeFeature4 should normalize the feature representations.
    Normalizations corresponding to the input features are defined within
    the function itself

mode (string, either "tsne" or "umap")
    Whether to do dimensionality reduction with t-SNE or UMAP. Note: UMAP
    requires that umapFileExchange (1.4.2 or later) be in the path

neighborparam (scalar)
    Primary parameter of dimensionality reduction; perplexity in t-SNE and
    n_neighbors in UMAP

only2d (bool)
    Whether to do dimensionality and reduction and visualization to 2
    dimensions, as opposed to do it to both 2 and 3

celltypes (string array)
    List of celltypes with which to annotate data in embeddingStruct and
    color points in visualization

celltypecategories (array of shape length(celltypes) x 3)
    RGB colors assigned to celltypes. Must be in same order as celltypes

lineages (string array)
    List of lineages with which to annotate data in embeddingStruct and
    color points in visualization

lineagecategories (array of shape length(lineages) x 3)
    RGB colors assigned to lineages. Must be in same order as lineages




%}

%%

%%% Inputs and Parameters %%%

%embryo loading

%WG embryo
emb = '/Users/evanvietorisz/Documents/MATLAB/WT_dispim_aug22/Decon_emb1_edited_v13_11202018namingfix.zip';%embryo directory
input_xy_res = .16; % resolution of xy axes at which emb was imaged (in microns)
anisotropy = 1; %note: Wormguides embryo is anisotropy = 1,inut_xy_res = .16
%note: find both in corresponding .xml file

%code to load other C. elegans embryos
%{
%other WT
%{

anisotropy = 5;
inut_xy_res = .15;
tend = 271;
emb = '/Users/evanvietorisz/Documents/MATLAB/WT_embryos_Aug7/20140407_JIM113_SiO-0.15_1_s1_emb_mirroredcorrect_edited.zip';
%}

%5 homeotic mutants

%{
%APX-1
% anisotropy = 1,inut_xy_res = .254
emb = '/Users/evanvietorisz/Documents/MATLAB/Mutants_Feb_12/ZD_RW10348_APX-1_20111104_2_s4_emb1_edited/ZD_RW10348_APX-1_20111104_2_s4_emb1_edited_xy_qc.zip';
input_xy_res = .254; % resolution of xy axes in microns
anisotropy = 1/input_xy_res;

tstart=1;
tend=180; %from excel spreadsheet saying what the last trustworthy time is (corresponds to 350 cell stage)
%}

%{
%MOM-1
% anisotropy = 1,inut_xy_res = .254
emb = '/Users/evanvietorisz/Documents/MATLAB/Mutants_Feb_12/ZD_RW10348_MOM-2_20110917_3_s3_emb1_edited/ZD_RW10348_MOM-2_20110917_3_s3_emb1_edited_xy_qc.zip';
input_xy_res = .254; % resolution of xy axes in microns
anisotropy = 1/input_xy_res;

tstart=1;
tend=190; %from excel spreadsheet saying what the last trustworthy time is (corresponds to 350 cell stage)
%}

%{
%PIE-1
% anisotropy = 1,inut_xy_res = .254
emb = '/Users/evanvietorisz/Documents/MATLAB/Mutants_Feb_12/ZD_RW10348_PIE-1_20110315_1_s3_emb1_edited/ZD_RW10348_PIE-1_20110315_1_s3_emb1_edited_xy_qc.zip';
input_xy_res = .254; % resolution of xy axes in microns
anisotropy = 1/input_xy_res;

tstart=1;
tend=164; %from excel spreadsheet saying what the last trustworthy time is (corresponds to 350 cell stage)
%}

%{
%SKN-1
% anisotropy = 1,inut_xy_res = .254
emb = '/Users/evanvietorisz/Documents/MATLAB/Mutants_Feb_12/ZD_RW10348_SKN-1_20110422_1_s3_emb1_edited/ZD_RW10348_SKN-1_20110422_1_s3_emb1_edited_xy_qc.zip';
input_xy_res = .254; % resolution of xy axes in microns
anisotropy = 1/input_xy_res;

tstart=1;
tend=170; %from excel spreadsheet saying what the last trustworthy time is (corresponds to 350 cell stage)
%}

%{
%WWP-1
% anisotropy = 1,inut_xy_res = .254
emb = '/Users/evanvietorisz/Documents/MATLAB/Mutants_Feb_12/ZD_RW10348_WWP-1_20111208_2_s3_emb1_edited/ZD_RW10348_WWP-1_20111208_2_s3_emb1_edited_xy_qc.zip';
input_xy_res = .254; % resolution of xy axes in microns
anisotropy = 1/input_xy_res;

tstart=1;
tend=168; %from excel spreadsheet saying what the last trustworthy time is (corresponds to 350 cell stage)
%}
%}

tend=411; %time through which raw data is loaded into embdat; should be when worm starts moving (determined by inspection)

%alignment
centering_time = tend; %timept used to align centroids of C. elegans embs in alignment step; usually same as tend

%computing features
starting_num_cells = 6; %number of cells at which to start *computing features*
sampling=2; % frequency with which to compute features on raw emb dataset; ex. sampling = 2 --> every other frame gets feature rep. computed

feats =      ["beginning_of_current_cell_cycle";
              "end_of_current_cell_cycle";       
              "duration_of_current_cell_cycle";  
              "time_until_next_division";        
              "time_since_last_division";        
              "fraction_of_current_cycle_elapsed";
              "similarity_to_ancestors_cell_cycle";          
              "similarity_to_cousins_cell_cycle";
              "first_deriv_x";                    
              "first_deriv_y";                    
              "first_deriv_z";                    
              "laplacian";                        
              "local_normal_vector";              
              "residual";                         
              "distances_to_k_nns"                  
              "signed_distances"                  
              "vec_to_position_x_mins_ago";                             
              "vec_to_position_x_mins_in_future";                             
              "vec_to_site_of_last_division";                            
              "vec_to_site_of_next_division";                             
              "displacement_vec_over_cell_cycle";                             
              "proportion_of_present_nns_that_are_descendants_of_past_nns"; 
              "filtered_nuc_size";                
              "number_of_cell_in_embryo";         
              "developmental_time"];   

%hyperparameters
jul2_5nns_lineage = {"radius_boundingsphere_celldensity",25;
    "radius_boundingsphere_residual",   10;
    "radius_boundingsphere_nucsize",    20;
    "k_intercelldistance",              5;
    "degree_cousins",                   1;
    "timewindow_migrationpaths",        .0625;
    "timewindow_migrationcongruency",   .0750;
    "k_migrationcongruency",            15};
%from optimization July 2 2021, optimized for lineage continuity, 
%scored using k = 5 nns

%other hyperparameter sets
%{
%other param sets from July 2 2021 optimization
%{
%optimizations 7/2
jul2_20nns_fate = {"radius_boundingsphere_celldensity",15;
    "radius_boundingsphere_residual",   10;
    "radius_boundingsphere_nucsize",    25;
    "k_intercelldistance",              5;
    "degree_cousins",                   1;
    "timewindow_migrationpaths",        .15;
    "timewindow_migrationcongruency",   .0375;
    "k_migrationcongruency",            5};

jul2_20nns_lineage = {"radius_boundingsphere_celldensity",25;
    "radius_boundingsphere_residual",   10;
    "radius_boundingsphere_nucsize",    5;
    "k_intercelldistance",              5;
    "degree_cousins",                   1;
    "timewindow_migrationpaths",        .15;
    "timewindow_migrationcongruency",   .0625;
    "k_migrationcongruency",            5};


jul2_5nns_fate = {"radius_boundingsphere_celldensity",25;
    "radius_boundingsphere_residual",   10;
    "radius_boundingsphere_nucsize",    25;
    "k_intercelldistance",              5;
    "degree_cousins",                   1;
    "timewindow_migrationpaths",        .1250;
    "timewindow_migrationcongruency",   .0625;
    "k_migrationcongruency",            15};
%}

% old sets of hyperparameters
%{
%hyperparameters 3/14
default_param_values = {"radius_boundingsphere_celldensity",20;
                        "radius_boundingsphere_residual",   20;
                        "radius_boundingsphere_nucsize",    20;
                        "k_intercelldistance",              50;
                        "degree_cousins",                   1;
                        "timewindow_migrationpaths",        .125;
                        "timewindow_migrationcongruency",   .125;
                        "k_migrationcongruency",            50};
  

%hyperparamters 3/24
lineage_optimized_param_values = {"radius_boundingsphere_celldensity",20;
    "radius_boundingsphere_residual",   10;
    "radius_boundingsphere_nucsize",    5;
    "k_intercelldistance",              5;
    "degree_cousins",                   1;
    "timewindow_migrationpaths",        .125;
    "timewindow_migrationcongruency",   .125;
    "k_migrationcongruency",            5};


%}
%}

params_to_use = jul2_5nns_fate;%actual assignment

normalize_bool = true;%whether or not to normalize


%dimensionality reduction and plotting

mode = "umap"; % "tsne" or "umap", note: umap requires umapFileExchange to be in path
neighborparam = 30; %perplexity in tsne and n_neighbors in umap

only2d = true;

%celltypes recognized for embeddingStruct and plotting
celltypes = ["pharynx","neuron","bodywallmuscle","hypoderm","gut"];

%RGB colors used to color different cell types in scatter plots
celltypecategories = [[0,1,0];[0,0,1];[1,0,0];[1,0,1];[0,1,1]];

%analogous to above
lineages = ["ABa","ABp","C","D","E","MS","P"];
lineagecategories = [[1,0,0];[0,1,0];[0,0,1];[1,0,1];[1,.5,0];[0,1,1];[.5,0,1]];


%%

%%% Load Raw Data %%%

templocation='temp_unzip/';
%unzip zipfile to temp file
if ~exist(emb,'file')
    errors.zipfilemissing=1;
    return
end
disp('the file exists')
try
    unzip(emb,templocation);
catch exception
    errors.zipfilecorrupted=1;
    return
end
disp('the file was unzipped')

%load embryo struct and cells struct
[cells,embdat] = loadcells_unnamed(templocation,tend,anisotropy,false); %note: cells data structure not used in this pipeline
rmdir(templocation,'s');

%%

%%% Alignment to Reference C. elegans Dataset %%%

%align input to reference
load('ref_emb.mat')
embdat_stabilized = coalignNamedEmbryosPerTime2(ref_emb,1,size(ref_emb,2),.16,embdat,1,tend,anisotropy,input_xy_res,centering_time);
%note: .16 is  the input_xy_res of the reference of the reference

disp('the input was transformed')

%update cells object
cells_stabilized = parseCellsFromEmb(embdat_stabilized,tend);

%backup
cells_backup=cells;
cells=cells_stabilized;
disp('the data was prepared')


%%

%%% Prepare Data and Make embeddingStruct %%%


%times through which to compute feature representations for all cells
starttime = 1;
while size(embdat_stabilized(starttime).finalpoints,1) < starting_num_cells
    starttime = starttime + 1;
end
endtime=tend-1; %should be one less than last timepoint loaded because loading code leaves an empty entry at the end of the Struct array

%make embeddingStruct
embeddingStruct = embeddingStructMaker(embdat_stabilized,starttime,sampling,endtime,lineages);

%remove any false nucs
nuc_inds = [];
for i = 1:length(embeddingStruct)
    if startsWith(embeddingStruct(i).name,'Nuc')
        nuc_inds = [nuc_inds;i];
    end
end
embeddingStruct(nuc_inds) = [];

numPointsInEmbedding = size(embeddingStruct,2);

%%

%%% Compute Features %%% 

times = [starttime:sampling:endtime];

[percellfeature,normed_percellfeature,~] = computeFeatures4(embdat_stabilized,times,feats,params_to_use,normalize_bool,input_xy_res);
disp('feature extraction worked')


%load into embeddingStruct
for i = 1:numPointsInEmbedding
    embeddingStruct(i).featurevector = percellfeature(i,:);
    embeddingStruct(i).normed_featurevector = normed_percellfeature(i,:);
end



%%

%%% Dimension Reduction %%%

if mode == "tsne"
    %2d
    reducedfeature=tsne(normed_percellfeature,'Perplexity',neighborparam);
    disp('tsne2d worked')
    %3d
    if only2d == false
        reducedfeature3=tsne(normed_percellfeature,'NumDimensions',3,'Perplexity',neighborparam);
        disp('tsne3d worked')
    end
elseif mode == "umap"
    %2d
    reducedfeature=run_umap(normed_percellfeature,'n_neighbors',neighborparam);
    disp('umap2d worked')
    %3d
    if only2d == false
        reducedfeature3=run_umap(normed_percellfeature,'n_components',3,'n_neighbors',neighborparam);
        disp('umap3d worked')
    end
else
    disp('mode must be "tsne" or "umap"')
end

%load into embeddingStruct
for i = 1:numPointsInEmbedding
    embeddingStruct(i).twoDpoint = reducedfeature(i,:);
    if only2d == false
        embeddingStruct(i).threeDpoint = reducedfeature3(i,:);
    end
end

%%

%%% Color Scheme %%%

%temporal evolution of embryo
timecolors=zeros(numPointsInEmbedding,3);
for i = 1:numPointsInEmbedding
    hsvcolor = [embeddingStruct(i).captureFrame/endtime,1,1];
    timecolors(i,:)=hsv2rgb(hsvcolor);
end

%cell type 

%default color gray
celltypecolors = repmat([.5,.5,.5],numPointsInEmbedding,1);

for i = 1:numPointsInEmbedding
    for j = 1:length(celltypes)
        if embeddingStruct(i).celltype == celltypes(j)
            celltypecolors(i,:) = celltypecategories(j,:);
        end
    end

end

%lineage 

%default color gray
lineagecolors = repmat([.25,.25,.25],numPointsInEmbedding,1);

for i = 1:numPointsInEmbedding
    for j = 1:length(lineages)
        if embeddingStruct(i).lineage == lineages(j)
            lineagecolors(i,:) = lineagecategories(j,:);
        end
    end
end

%%

%%% Plotting %%%

figure('Name','Embeddings','NumberTitle','off')
figtitle = 'XXXX';
sgtitle(figtitle)

startcells = embeddingStruct(1).numCellsAtCapture;
endcells = embeddingStruct(numPointsInEmbedding).numCellsAtCapture;

startcellsString = num2str(startcells);
endcellsString = num2str(endcells);

timetitle = 'Time Evolution';
timetitle = [timetitle newline 'starting number of cells = ' startcellsString ' | ending number of cells = ' endcellsString];

difftitle = 'Cell Differentiation';
difftitle = [difftitle newline 'Pharynx (g), Neuron (b), Body Wall Muscle (r), Hypoderm (f), Gut (c)'];

lintitle = 'Lineages';
lintitle = [lintitle newline 'ABa (r), ABp (g), C (b), D (f), E (o), M (c), P (p)'];


%time evolution
if only2d == false
    subplot(2,3,1)
else
    subplot(1,3,1)
end
scatter(reducedfeature(:,1),reducedfeature(:,2),10,timecolors(1:length(reducedfeature),:));
pbaspect([1,1,1])
title(timetitle)
if mode == "tsne"
    xlabel('tsne-1')
    ylabel('tsne-2')
elseif mode == "umap"
    xlabel('umap-1')
    ylabel('umap-2')
end

if only2d == false
    subplot(2,3,4)
    scatter3(reducedfeature3(:,1),reducedfeature3(:,2),reducedfeature3(:,3),10,timecolors(1:length(reducedfeature),:));
    pbaspect([1,1,1])
    title(timetitle)
    if mode == "tsne"
        xlabel('tsne-1')
        ylabel('tsne-2')
        zlabel('tsne-3')
    elseif mode == "umap"
        xlabel('umap-1')
        ylabel('umap-2')
        zlabel('umap-3')
    end
end

%cell differentiation

if only2d == false
    subplot(2,3,2)
else
    subplot(1,3,2)
end
hold on
scatter(reducedfeature(:,1),reducedfeature(:,2),10,celltypecolors(1:length(reducedfeature),:));
pbaspect([1,1,1])

title(difftitle)
if mode == "tsne"
    xlabel('tsne-1')
    ylabel('tsne-2')
elseif mode == "umap"
    xlabel('umap-1')
    ylabel('umap-2')
end
hold off

if only2d == false
    subplot(2,3,5)
    scatter3(reducedfeature3(:,1),reducedfeature3(:,2),reducedfeature3(:,3),10,celltypecolors(1:length(reducedfeature),:));
    pbaspect([1,1,1])
    title(difftitle)
    if mode == "tsne"
        xlabel('tsne-1')
        ylabel('tsne-2')
        zlabel('tsne-3')
    elseif mode == "umap"
        xlabel('umap-1')
        ylabel('umap-2')
        zlabel('umap-3')
    end
end

if only2d == false
    subplot(2,3,3)
else
    subplot(1,3,3)
end
scatter(reducedfeature(:,1),reducedfeature(:,2),10,lineagecolors(1:length(reducedfeature),:));
pbaspect([1,1,1])
title(lintitle)
if mode == "tsne"
    xlabel('tsne-1')
    ylabel('tsne-2')
elseif mode == "umap"
    xlabel('umap-1')
    ylabel('umap-2')
end

if only2d == false
    subplot(2,3,6)
    scatter3(reducedfeature3(:,1),reducedfeature3(:,2),reducedfeature3(:,3),10,lineagecolors(1:length(reducedfeature),:));
    pbaspect([1,1,1])
    title(lintitle)
    if mode == "tsne"
        xlabel('tsne-1')
        ylabel('tsne-2')
        zlabel('tsne-3')
    elseif mode == "umap"
        xlabel('umap-1')
        ylabel('umap-2')
        zlabel('umap-3')
    end
end


return
