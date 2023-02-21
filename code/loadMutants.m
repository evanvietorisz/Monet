%function for loading and saving mutant embryos onto local disk as .mat
%files

%note: throughout this, mom-2 is incorrectly called mom-1

%%
%APX-1

close all
clear all

% anisotropy = 1,inut_xy_res = .16
emb = '/Users/evanvietorisz/Documents/MATLAB/WT_embryos_Aug7/Decon_emb1_MGedits.zip';
input_xy_res = .16; % resolution of xy axes in microns
anisotropy = 1;

tstart=1;
tend=305; %cutoff should be when worm starts moving around (determined manually)


try
    [percellfeature,normed_percellfeature,embeddingStruct,embdat_stabilized,starttime,sampling,endtime] = pipeline(emb,input_xy_res,anisotropy,tstart,tend);

    
    save('WT_aug8.mat','anisotropy','embdat_stabilized','embeddingStruct','input_xy_res','starttime','sampling','endtime','tstart','tend','normed_percellfeature','percellfeature')
    
catch
    disp('error on WT')
end





%%
%APX-1

close all
clear all

% anisotropy = 1,inut_xy_res = .254
emb = '/Users/evanvietorisz/Documents/MATLAB/Mutants_Feb_12/ZD_RW10348_APX-1_20111104_2_s4_emb1_edited/ZD_RW10348_APX-1_20111104_2_s4_emb1_edited_xy_qc.zip';
input_xy_res = .254; % resolution of xy axes in microns
anisotropy = 1/input_xy_res;

tstart=1;
tend=180; %from excel spreadsheet saying what the last trustworthy time is (corresponds to 350 cell stage)
%}


try
    [percellfeature,normed_percellfeature,embeddingStruct,embdat_stabilized,starttime,sampling,endtime] = pipeline(emb,input_xy_res,anisotropy,tstart,tend);

    
    save('APX1_aug8.mat','anisotropy','embdat_stabilized','embeddingStruct','input_xy_res','starttime','sampling','endtime','tstart','tend','normed_percellfeature','percellfeature')
    
catch
    disp('error on APX1')
end
    
%%


%MOM-2

close all
clear all

% anisotropy = 1,inut_xy_res = .254
emb = '/Users/evanvietorisz/Documents/MATLAB/Mutants_Feb_12/ZD_RW10348_MOM-2_20110917_3_s3_emb1_edited/ZD_RW10348_MOM-2_20110917_3_s3_emb1_edited_xy_qc.zip';
input_xy_res = .254; % resolution of xy axes in microns
anisotropy = 1/input_xy_res;

tstart=1;
tend=190; %from excel spreadsheet saying what the last trustworthy time is (corresponds to 350 cell stage)
%}

try
    [percellfeature,normed_percellfeature,embeddingStruct,embdat_stabilized,starttime,sampling,endtime] = pipeline(emb,input_xy_res,anisotropy,tstart,tend);
    
    save('MOM2_aug8.mat','anisotropy','embdat_stabilized','embeddingStruct','input_xy_res','starttime','sampling','endtime','tstart','tend','normed_percellfeature','percellfeature')
    
catch
    disp('error on MOM2')
end

%%

%PIE-1

close all
clear all


% anisotropy = 1,inut_xy_res = .254
emb = '/Users/evanvietorisz/Documents/MATLAB/Mutants_Feb_12/ZD_RW10348_PIE-1_20110315_1_s3_emb1_edited/ZD_RW10348_PIE-1_20110315_1_s3_emb1_edited_xy_qc.zip';
input_xy_res = .254; % resolution of xy axes in microns
anisotropy = 1/input_xy_res;

tstart=1;
tend=164; %from excel spreadsheet saying what the last trustworthy time is (corresponds to 350 cell stage)
%}


try
    [percellfeature,normed_percellfeature,embeddingStruct,embdat_stabilized,starttime,sampling,endtime] = pipeline(emb,input_xy_res,anisotropy,tstart,tend);

    save('PIE1_aug8.mat','anisotropy','embdat_stabilized','embeddingStruct','input_xy_res','starttime','sampling','endtime','tstart','tend','normed_percellfeature','percellfeature')
    
catch
    disp('error on PIE1')
end

%%

%SKN-1

close all
clear all

% anisotropy = 1,inut_xy_res = .254
emb = '/Users/evanvietorisz/Documents/MATLAB/Mutants_Feb_12/ZD_RW10348_SKN-1_20110422_1_s3_emb1_edited/ZD_RW10348_SKN-1_20110422_1_s3_emb1_edited_xy_qc.zip';
input_xy_res = .254; % resolution of xy axes in microns
anisotropy = 1/input_xy_res;

tstart=1;
tend=170; %from excel spreadsheet saying what the last trustworthy time is (corresponds to 350 cell stage)
%}

try
    [percellfeature,normed_percellfeature,embeddingStruct,embdat_stabilized,starttime,sampling,endtime] = pipeline(emb,input_xy_res,anisotropy,tstart,tend);

    save('SKN1_aug8.mat','anisotropy','embdat_stabilized','embeddingStruct','input_xy_res','starttime','sampling','endtime','tstart','tend','normed_percellfeature','percellfeature')
    
catch
    disp('error on SKN1')
end



%%

%WWP-1

close all
clear all

% anisotropy = 1,inut_xy_res = .254
emb = '/Users/evanvietorisz/Documents/MATLAB/Mutants_Feb_12/ZD_RW10348_WWP-1_20111208_2_s3_emb1_edited/ZD_RW10348_WWP-1_20111208_2_s3_emb1_edited_xy_qc.zip';
input_xy_res = .254; % resolution of xy axes in microns
anisotropy = 1/input_xy_res;

tstart=1;
tend=168; %from excel spreadsheet saying what the last trustworthy time is (corresponds to 350 cell stage)
%}

try
    [percellfeature,normed_percellfeature,embeddingStruct,embdat_stabilized,starttime,sampling,endtime] = pipeline(emb,input_xy_res,anisotropy,tstart,tend);

    save('WWP1_aug8.mat','anisotropy','embdat_stabilized','embeddingStruct','input_xy_res','starttime','sampling','endtime','tstart','tend','normed_percellfeature','percellfeature')
    
catch
    disp('error on WWP1')
end
%%


%PIE-1 %%note: to load this properly, you need to go to computeFeatures4
%and change the threefiftycelltime threshold to 348

close all
clear all


% anisotropy = 1,inut_xy_res = .254
emb = '/Users/evanvietorisz/Documents/MATLAB/Mutants_Feb_12/ZD_RW10348_PIE-1_20110315_1_s3_emb1_edited/ZD_RW10348_PIE-1_20110315_1_s3_emb1_edited_xy_qc.zip';
input_xy_res = .254; % resolution of xy axes in microns
anisotropy = 1/input_xy_res;

tstart=1;
tend=164; %from excel spreadsheet saying what the last trustworthy time is (corresponds to 350 cell stage)
%}

try
    [percellfeature,normed_percellfeature,embeddingStruct,embdat_stabilized,starttime,sampling,endtime] = pipeline(emb,input_xy_res,anisotropy,tstart,tend);

    save('PIE1_aug8.mat','anisotropy','embdat_stabilized','embeddingStruct','input_xy_res','starttime','sampling','endtime','tstart','tend','normed_percellfeature','percellfeature')
    
catch
    disp('error on PIE1')
end
