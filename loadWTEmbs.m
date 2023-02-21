% script to load the three dispim embryos with one command


%Embryo 1
%{

% anisotropy = 1,inut_xy_res = .16
emb = '/Users/evanvietorisz/Documents/MATLAB/WT_dispim_aug22/Decon_emb1_edited_v13_11202018namingfix.zip';
input_xy_res = .16; % resolution of xy axes in microns
anisotropy = 1;

tstart=1;
tend=277; %cutoff should be when worm starts moving around (determined manually)
%}

%Embryo 2
%{

%Embryo 2

% anisotropy = 1,inut_xy_res = .16
emb = '/Users/evanvietorisz/Documents/MATLAB/WT_dispim_aug22/Decon_emb1_MGedits.zip';
input_xy_res = .16; % resolution of xy axes in microns
anisotropy = 1;

tstart=1;
tend=305; %cutoff should be when worm starts moving around (determined manually)

%}

%Embryo 3
%{

% anisotropy = 1,inut_xy_res = .16
emb = '/Users/evanvietorisz/Documents/MATLAB/WT_dispim_aug22/Decon_emb1_updated.zip';
input_xy_res = .16; % resolution of xy axes in microns
anisotropy = 1;

tstart=1;
tend=261; %cutoff should be when worm starts moving around (determined manually)

%}


%%


%1

close all
clear all

% anisotropy = 1,inut_xy_res = .16
emb = '/Users/evanvietorisz/Documents/MATLAB/WT_dispim_aug22/Decon_emb1_edited_v13_11202018namingfix.zip';
input_xy_res = .16; % resolution of xy axes in microns
anisotropy = 1;

tstart=1;
tend=379; %cutoff should be when worm starts moving around (determined manually)


try
    [percellfeature,normed_percellfeature,embeddingStruct,embdat_stabilized,starttime,sampling,endtime] = pipeline(emb,input_xy_res,anisotropy,tstart,tend);

    
    save('dispim_emb1_full.mat','anisotropy','embdat_stabilized','embeddingStruct','input_xy_res','starttime','sampling','endtime','tstart','tend','normed_percellfeature','percellfeature')
    
catch
    disp('error on emb1')
end

%%
%2

close all
clear all

% anisotropy = 1,inut_xy_res = .16
emb = '/Users/evanvietorisz/Documents/MATLAB/WT_dispim_aug22/Decon_emb1_MGedits.zip';
input_xy_res = .16; % resolution of xy axes in microns
anisotropy = 1;

tstart=1;
tend=359; %cutoff should be when worm starts moving around (determined manually)


try
    [percellfeature,normed_percellfeature,embeddingStruct,embdat_stabilized,starttime,sampling,endtime] = pipeline(emb,input_xy_res,anisotropy,tstart,tend);

    
    save('dispim_emb2_full.mat','anisotropy','embdat_stabilized','embeddingStruct','input_xy_res','starttime','sampling','endtime','tstart','tend','normed_percellfeature','percellfeature')
    
catch
    disp('error on emb2')
end
%%

%3

close all
clear all

% anisotropy = 1,inut_xy_res = .16
emb = '/Users/evanvietorisz/Documents/MATLAB/WT_dispim_aug22/Decon_emb1_updated.zip'; 
input_xy_res = .16; % resolution of xy axes in microns
anisotropy = 1;

tstart=1;
tend=330; %cutoff should be when worm starts moving around (determined manually)


try
    [percellfeature,normed_percellfeature,embeddingStruct,embdat_stabilized,starttime,sampling,endtime] = pipeline(emb,input_xy_res,anisotropy,tstart,tend);

    
    save('dispim_emb3_full.mat','anisotropy','embdat_stabilized','embeddingStruct','input_xy_res','starttime','sampling','endtime','tstart','tend','normed_percellfeature','percellfeature')
    
catch
    disp('error on emb3')
end



%%
%{
figure

subplot(3,1,1)

embdat = embdat1;
num_cells = zeros(length(embdat),1);
for i = 1:length(num_cells)  
    num_cells(i) = size(embdat(i).finalpoints,1);
end

plot(num_cells,'*')
xlabel('frame')
ylabel('num cells')
pbaspect([1,1,1])

subplot(3,1,2)

embdat = embdat2;
num_cells = zeros(length(embdat),1);
for i = 1:length(num_cells)  
    num_cells(i) = size(embdat(i).finalpoints,1);
end

plot(num_cells,'*')
xlabel('frame')
ylabel('num cells')
pbaspect([1,1,1])

subplot(3,1,3)

embdat = embdat3;
num_cells = zeros(length(embdat),1);
for i = 1:length(num_cells)  
    num_cells(i) = size(embdat(i).finalpoints,1);
end

plot(num_cells,'*')
xlabel('frame')
ylabel('num cells')
pbaspect([1,1,1])


sgtitle('number of cells with time in 3 dispim embryos')

%}


function [percellfeature,normed_percellfeature,embeddingStruct,embdat_stabilized,starttime,sampling,endtime] = pipeline(emb,input_xy_res,anisotropy,tstart,tend)

%extra

centering_time = tend; %last trustable timestep in dataset--should be same as tend
sampling=2;
endtime=tend-1; %should be one less than last timepoint loaded because of implementation

%features 2/15
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
          
          
%{          
%features that don't require tracking mar 4
feats =      ["first_deriv_x";                    
              "first_deriv_y";                    
              "first_deriv_z";                    
              "laplacian";                        
              "local_normal_vector";              
              "residual";                         
              "distances_to_k_nns";                  
              "signed_distances";                                              
              "filtered_nuc_size";                
              "number_of_cell_in_embryo";         
              "developmental_time"];
%}

              
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


%old fish
%{
%jank fish 10/11/20
emb = '/Users/evanvietorisz/Documents/MATLAB/other/stack_emb_forcedfancythresh2_6001_edited.zip';
tend = 599;
centering_time = 599;
anisotropy = 5.18;
%set normalize to false
%}
%old mutants
%{
%other WT all anisotropy = 5, inut_xy_res = .15
%emb = '/Users/evanvietorisz/Documents/MATLAB/WT_embryos_Aug7/20140407_JIM113_SiO-0.15_1_s1_emb_mirroredcorrect_edited.zip';
%emb = '/Users/evanvietorisz/Documents/MATLAB/WT_embryos_Aug7/20140407_JIM113_SiO-0.15_1_s2_emb1_edited.zip';
%emb = '/Users/evanvietorisz/Documents/MATLAB/WT_embryos_Aug7/20140407_JIM113_SiO-0.15_1_s3_emb1_edited.zip';

%pal-1 mutants anisoropy = 1/.254, xy_res = .254
%emb = '/Users/evanvietorisz/Documents/MATLAB/pal-1_mutants_Aug18/ZD_RW10348_PAL-1_20110318_2_s3_emb1_edited_xy.zip'; %155, 366 cells
%emb ='/Users/evanvietorisz/Documents/MATLAB/pal-1_mutants_Aug18/ZD_RW10348_PAL-1_20110318_2_s3_emb2_edited_xy.zip'; %%167, 359 cells*edit this

%DID YOU CHANGE THE ANISOTROPY,RES,CENTERING TIME?

%mutant 1 is reliable until t = 155, at which point it has 366 cells
%reference has 364 cells at its t = 208

%mutant 2 is reliable until t = 167, at which point it has 372 cells
%reference has 364 cells at its t = 213
%}


%starttime=101; %note: for ref emb t=50 is 24 points, t=100 is 87 points,t=230 is 410 points
%made this just come out of specifying how many cells you want to start
%with

starting_num_cells = 6; %number of cells at which to start computing features--must be one greater than k


normalize_bool = true;

mode = "umap"; % "tsne" or "umap", note: umap requires umapFileExchange to be in path
neighborparam = 30; %perplexity in tsne and n_neighbors in umap

only2d = true;

%first item in title of saved figure
%note: don't use underscores
feature_description = 'fish distfeat last entry test';

%load and prepare data
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
[cells,embdat] = loadcells_unnamed(templocation,tend,anisotropy,false);
rmdir(templocation,'s');

%related to old fish
%{
%jank take out anisotropy for fish embryo 10/11/20
embdat_stabilized = embdat;

%note: last frame in struct generally contains empty point array
for frame  = 1:length(embdat_stabilized)-1
    %remove anisotropy
    embdat_stabilized(frame).finalpoints(:,3) = embdat_stabilized(frame).finalpoints(:,3).*anisotropy;
end

%center input 
cent = [mean(embdat_stabilized(centering_time).finalpoints(:,1)),mean(embdat_stabilized(centering_time).finalpoints(:,2)),mean(embdat_stabilized(centering_time).finalpoints(:,3))];
for i = 1:tend
    embdat_stabilized(i).finalpoints = embdat_stabilized(i).finalpoints - cent;
end
%}

%%

%align input to reference (reference is NOT internally aligned)
%note the .16 is  scale of the reference and is just a property of the data
load('ref_emb.mat')
embdat_stabilized = coalignNamedEmbryosPerTime2(ref_emb,1,size(ref_emb,2),.16,embdat,1,tend,anisotropy,input_xy_res,centering_time);

%note: anisotropy appears not to be necessary for the alignment-- no diff
%difference
disp('the input was transformed')

%update cells object
cells_stabilized = parseCellsFromEmb(embdat_stabilized,tend);

%backup
cells_backup=cells;
cells=cells_stabilized;
disp('the data was prepared')


%%
%start at specified number of cells***abstract this away
starttime = 1;
while size(embdat_stabilized(starttime).finalpoints,1) < starting_num_cells
    starttime = starttime + 1;
end


%embeddingStruct
celltypes = ["pharynx","neuron","bodywallmuscle","hypoderm","gut"];
lineages = ["ABa","ABp","C","D","E","MS","P"];

embeddingStruct = embeddingStructMaker(embdat_stabilized,starttime,sampling,endtime,lineages);


%remove any debris/false nucs
nuc_inds = [];
for i = 1:length(embeddingStruct)
    if startsWith(embeddingStruct(i).name,'Nuc')
        nuc_inds = [nuc_inds;i];
    end
end
embeddingStruct(nuc_inds) = [];


numPointsInEmbedding = size(embeddingStruct,2);


%%

%aug 8 
jul2_20nns_lineage = {"radius_boundingsphere_celldensity",25;
    "radius_boundingsphere_residual",   10;
    "radius_boundingsphere_nucsize",    5;
    "k_intercelldistance",              5;
    "degree_cousins",                   1;
    "timewindow_migrationpaths",        .15;
    "timewindow_migrationcongruency",   .0625;
    "k_migrationcongruency",            5};





%%% compute featured %%% 
times = [starttime:sampling:endtime];
disp(starttime)
disp(sampling)
disp(endtime)
params_to_use = jul2_20nns_lineage;

[percellfeature,normed_percellfeature,~] = computeFeatures4(embdat_stabilized,times,feats,params_to_use,normalize_bool,input_xy_res);
disp('feature extraction worked')



end











