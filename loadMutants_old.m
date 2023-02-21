
%loading and saving mutant datasets
%in function below, start at 21 cells



%%%note: throughout this, mom-2 is incorrectly called mom-1




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
    [percellfeatureAPX1,normed_percellfeatureAPX1,embeddingStructAPX1] = pipeline(emb,input_xy_res,anisotropy,tstart,tend);

    save('percellfeatureAPX1_mar25','percellfeatureAPX1','normed_percellfeatureAPX1','embeddingStructAPX1')
catch
    disp('error on APX1')
end
    



%MOM-1

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
    [percellfeatureMOM1,normed_percellfeatureMOM1,embeddingStructMOM1] = pipeline(emb,input_xy_res,anisotropy,tstart,tend);

    save('percellfeatureMOM1_mar25','percellfeatureMOM1','normed_percellfeatureMOM1','embeddingStructMOM1')
catch
    disp('error on MOM1')
end



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
    [percellfeaturePIE1,normed_percellfeaturePIE1,embeddingStructPIE1] = pipeline(emb,input_xy_res,anisotropy,tstart,tend);

    save('percellfeaturePIE1_mar25','percellfeaturePIE1','normed_percellfeaturePIE1','embeddingStructPIE1')
catch
    disp('error on PIE1')
end



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
    [percellfeatureSKN1,normed_percellfeatureSKN1,embeddingStructSKN1] = pipeline(emb,input_xy_res,anisotropy,tstart,tend);

    save('percellfeatureSKN1_mar25','percellfeatureSKN1','normed_percellfeatureSKN1','embeddingStructSKN1')
catch
    disp('error on SKN1')
end

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
    [percellfeatureWWP1,normed_percellfeatureWWP1,embeddingStructWWP1] = pipeline(emb,input_xy_res,anisotropy,tstart,tend);

    save('percellfeatureWWP1_mar25','percellfeatureWWP1','normed_percellfeatureWWP1','embeddingStructWWP1')
catch
    disp('error on WWP1')
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



[percellfeaturePIE1,normed_percellfeaturePIE1,embeddingStructPIE1] = pipeline(emb,input_xy_res,anisotropy,tstart,tend);

save('percellfeaturePIE1_mar25','percellfeaturePIE1','normed_percellfeaturePIE1','embeddingStructPIE1')

%%

%plot the mutants and WT



%load the WT, cut down to 350 cell
%load mutants
%create color arrays for each wt

%do the 5 umaps, plot


%load percellfeatures, embeddingStructs
load('percellfeatureWT_mar25.mat')

%pare down to 350 cells
t = 1;
while embeddingStruct(t).numCellsAtCapture <= 350
    t = t + 1;
end
percellfeature = percellfeature(1:t,:);
normed_percellfeature = normed_percellfeature(1:t,:);
embeddingStruct = embeddingStruct(1:t);


load('percellfeatureAPX1_mar25.mat')
load('percellfeatureMOM1_mar25.mat')
load('percellfeaturePIE1_mar25.mat')
load('percellfeatureSKN1_mar25.mat')
load('percellfeatureWWP1_mar25.mat')

celltypes = ["pharynx","neuron","bodywallmuscle","hypoderm","gut"];
lineages = ["ABa","ABp","C","D","E","MS","P"];

%generate colors for UMAPS
[celltypecolorsWT,lineagecolorsWT,timecolorsWT] = colorArrays(embeddingStruct,celltypes,lineages); %embeddingStruct just called embeddingStruct for WT
[celltypecolorsAPX1,lineagecolorsAPX1,timecolorsAPX1] = colorArrays(embeddingStructAPX1,celltypes,lineages);
[celltypecolorsMOM1,lineagecolorsMOM1,timecolorsMOM1] = colorArrays(embeddingStructMOM1,celltypes,lineages);
[celltypecolorsPIE1,lineagecolorsPIE1,timecolorsPIE1] = colorArrays(embeddingStructPIE1,celltypes,lineages);
[celltypecolorsSKN1,lineagecolorsSKN1,timecolorsSKN1] = colorArrays(embeddingStructSKN1,celltypes,lineages);
[celltypecolorsWWP1,lineagecolorsWWP1,timecolorsWWP1] = colorArrays(embeddingStructWWP1,celltypes,lineages);


%create umaps
%define two lineages to put in
%%
feature1 = normed_percellfeature; %note: percellfeature doesn't have WT in its name
name1 = 'WT';
celltypecolors1 = celltypecolorsWT;
lineagecolors1 = lineagecolorsWT;
timecolors1 = timecolorsWT;


feature2 = normed_percellfeatureWWP1;
celltypecolors2 = celltypecolorsWWP1;
lineagecolors2 = lineagecolorsWWP1;
timecolors2 = timecolorsWWP1;
name2 = 'WWP-1';

[reduced1,reduced2] = codecompose(feature1,feature2,neighborparam);


%plot

figure

%fate
subplot(1,3,1)
hold on
scatter(reduced1(:,1),reduced1(:,2),20,celltypecolors1)
scatter(reduced2(:,1),reduced2(:,2),6,celltypecolors2,'filled','MarkerEdgeColor',[0,0,0])
hold off
xlabel('UMAP-1')
ylabel('UMAP-2')
title('Pharynx (g), Neuron (b), Body Wall Muscle (r), Hypoderm (f), Gut (c)')

%time
subplot(1,3,2)
hold on
scatter(reduced1(:,1),reduced1(:,2),20,timecolors1)
scatter(reduced2(:,1),reduced2(:,2),6,timecolors2,'filled','MarkerEdgeColor',[0,0,0])
hold off
xlabel('UMAP-1')
ylabel('UMAP-2')
title('time: green -> red (51 cells to 571 cells)')

%lineage
subplot(1,3,3)
hold on
scatter(reduced1(:,1),reduced1(:,2),20,lineagecolors1)
scatter(reduced2(:,1),reduced2(:,2),6,lineagecolors2,'filled','MarkerEdgeColor',[0,0,0])
hold off
xlabel('UMAP-1')
ylabel('UMAP-2')
title('lineages: ABa (r), ABp (g), C (b), D (f), E (o), M (c), P (p)')

tt = ['co-decomposition of ',name1,' (rings) and ',name2,' (dots)'];
sgtitle(tt)

%%

%would probably make sense to make embeddings of the WT and mutants by
%themseselves

%uses things defined above

reducedWT = run_umap([normed_percellfeature],'n_neighbors',neighborparam);
reducedAPX1 = run_umap([normed_percellfeatureAPX1],'n_neighbors',neighborparam);
reducedMOM1 = run_umap([normed_percellfeatureMOM1],'n_neighbors',neighborparam);
reducedPIE1 = run_umap([normed_percellfeaturePIE1],'n_neighbors',neighborparam);
reducedSKN1 = run_umap([normed_percellfeatureSKN1],'n_neighbors',neighborparam);
reducedWWP1 = run_umap([normed_percellfeatureWWP1],'n_neighbors',neighborparam);


%%
figure

subplot(2,6,1)
scatter(reducedWT(:,1),reducedWT(:,2),20,celltypecolorsWT)
xlabel('UMAP-1')
ylabel('UMAP-2')
title('WT')

subplot(2,6,2)
scatter(reducedAPX1(:,1),reducedAPX1(:,2),20,celltypecolorsAPX1)
xlabel('UMAP-1')
ylabel('UMAP-2')
title('APX-1')

subplot(2,6,3)
scatter(reducedMOM1(:,1),reducedMOM1(:,2),20,celltypecolorsMOM1)
xlabel('UMAP-1')
ylabel('UMAP-2')
title('MOM-1')

subplot(2,6,4)
scatter(reducedPIE1(:,1),reducedPIE1(:,2),20,celltypecolorsPIE1)
xlabel('UMAP-1')
ylabel('UMAP-2')
title('PIE-1')

subplot(2,6,5)
scatter(reducedSKN1(:,1),reducedSKN1(:,2),20,celltypecolorsSKN1)
xlabel('UMAP-1')
ylabel('UMAP-2')
title('SKN-1')

subplot(2,6,6)
scatter(reducedWWP1(:,1),reducedWWP1(:,2),20,celltypecolorsWWP1)
xlabel('UMAP-1')
ylabel('UMAP-2')
title('WWP-1')

%

subplot(2,6,7)
scatter(reducedWT(:,1),reducedWT(:,2),20,lineagecolorsWT)
xlabel('UMAP-1')
ylabel('UMAP-2')


subplot(2,6,8)
scatter(reducedAPX1(:,1),reducedAPX1(:,2),20,lineagecolorsAPX1)
xlabel('UMAP-1')
ylabel('UMAP-2')


subplot(2,6,9)
scatter(reducedMOM1(:,1),reducedMOM1(:,2),20,lineagecolorsMOM1)
xlabel('UMAP-1')
ylabel('UMAP-2')


subplot(2,6,10)
scatter(reducedPIE1(:,1),reducedPIE1(:,2),20,lineagecolorsPIE1)
xlabel('UMAP-1')
ylabel('UMAP-2')


subplot(2,6,11)
scatter(reducedSKN1(:,1),reducedSKN1(:,2),20,lineagecolorsSKN1)
xlabel('UMAP-1')
ylabel('UMAP-2')


subplot(2,6,12)
scatter(reducedWWP1(:,1),reducedWWP1(:,2),20,lineagecolorsWWP1)
xlabel('UMAP-1')
ylabel('UMAP-2')

tt = 'individual decompositions of WT and 5 mutants (time = 21 cells to 350 cells)';
tt = [tt newline 'Pharynx (g), Neuron (b), Body Wall Muscle (r), Hypoderm (f), Gut (c)'];
tt = [tt newline 'lineages: ABa (r), ABp (g), C (b), D (f), E (o), M (c), P (p)'];

sgtitle(tt)


%%
%comparing overlap of 2 WTs


percellfeatureWT2 = normed_percellfeature;
embeddingStrutWT2 = embeddingStruct;

load('percellfeatureWT_mar25.mat')

percellfeatureWT1 = normed_percellfeature;
embeddingStrutWT1 = embeddingStruct;

[celltypecolors2,lineagecolors2,timecolors2] = colorArrays(embeddingStrutWT2,celltypes,lineages);
[celltypecolors1,lineagecolors1,timecolors1] = colorArrays(embeddingStrutWT1,celltypes,lineages);

[reduced1,reduced2] = codecompose(percellfeatureWT1,percellfeatureWT2,neighborparam);


figure


%fate
subplot(1,3,1)
hold on
scatter(reduced1(:,1),reduced1(:,2),20,celltypecolors1)
scatter(reduced2(:,1),reduced2(:,2),6,celltypecolors2,'filled','MarkerEdgeColor',[0,0,0])
hold off
xlabel('UMAP-1')
ylabel('UMAP-2')
title('Pharynx (g), Neuron (b), Body Wall Muscle (r), Hypoderm (f), Gut (c)')

%time
subplot(1,3,2)
hold on
scatter(reduced1(:,1),reduced1(:,2),20,timecolors1)
scatter(reduced2(:,1),reduced2(:,2),6,timecolors2,'filled','MarkerEdgeColor',[0,0,0])
hold off
xlabel('UMAP-1')
ylabel('UMAP-2')
title('time: green -> red (21 cells to 571 cells)')

%lineage
subplot(1,3,3)
hold on
scatter(reduced1(:,1),reduced1(:,2),20,lineagecolors1)
scatter(reduced2(:,1),reduced2(:,2),6,lineagecolors2,'filled','MarkerEdgeColor',[0,0,0])
hold off
xlabel('UMAP-1')
ylabel('UMAP-2')
title('lineages: ABa (r), ABp (g), C (b), D (f), E (o), M (c), P (p)')

sgtitle('codecomposition of 2 WTs to check temporal overlap')


%%
%removing weird outlier cells from PIE1
%{

% removing weird cell from PIE1
weird_inds = [];
for i = 1:size(normed_percellfeaturePIE1,1)
    
    if reducedPIE1(i,1) < -120
        weird_inds = [weird_inds;i];
    end
    
end



ebstr = embeddingStructPIE1(weird_inds);



not_weird_inds = setdiff([1:length(normed_percellfeaturePIE1)]',weird_inds);

not_weird_percellfeaturePIE1 = normed_percellfeaturePIE1(not_weird_inds,:);
not_weird_reducedPIE1 = run_umap([not_weird_percellfeaturePIE1],'n_neighbors',neighborparam);


figure 

subplot(1,2,1)
scatter(not_weird_reducedPIE1(:,1),not_weird_reducedPIE1(:,2),20,celltypecolorsPIE1(not_weird_inds,:))
xlabel('UMAP-1')
ylabel('UMAP-2')


subplot(1,2,2)
scatter(not_weird_reducedPIE1(:,1),not_weird_reducedPIE1(:,2),20,lineagecolorsPIE1(not_weird_inds,:))
xlabel('UMAP-1')
ylabel('UMAP-2')

sgtitle('redoing PIE-1 without outlier cells')






[celltypecolorsWT,lineagecolorsWT,timecolorsWT] = colorArrays(embeddingStruct,celltypes,lineages);


feature1 = normed_percellfeature; %note: percellfeature doesn't have WT in its name
name1 = 'WT';
celltypecolors1 = celltypecolorsWT;
lineagecolors1 = lineagecolorsWT;
timecolors1 = timecolorsWT;


feature2 = normed_percellfeaturePIE1;
celltypecolors2 = celltypecolorsPIE1;
lineagecolors2 = lineagecolorsPIE1;
timecolors2 = timecolorsPIE1;
name2 = 'PIE-1';

[reduced1,reduced2] = codecompose(feature1,feature2,neighborparam);


%plot

figure

%fate
subplot(1,2,1)
hold on
scatter(reduced1(:,1),reduced1(:,2),20,celltypecolors1)
scatter(reduced2(:,1),reduced2(:,2),6,celltypecolors2,'filled','MarkerEdgeColor',[0,0,0])
hold off
xlabel('UMAP-1')
ylabel('UMAP-2')
title('Pharynx (g), Neuron (b), Body Wall Muscle (r), Hypoderm (f), Gut (c)')

%lineage

subplot(1,2,2)
hold on
scatter(reduced1(:,1),reduced1(:,2),20,lineagecolors1)
scatter(reduced2(:,1),reduced2(:,2),6,lineagecolors2,'filled','MarkerEdgeColor',[0,0,0])
hold off
xlabel('UMAP-1')
ylabel('UMAP-2')
title('Pharynx (g), Neuron (b), Body Wall Muscle (r), Hypoderm (f), Gut (c)')


sgtitle('comparison of WT and PIE-1 after removing outlier cells')
%}
%%

%decomposing mutants and WT together

% removing weird cell from PIE1
reducedPIE1 = run_umap([normed_percellfeaturePIE1],'n_neighbors',neighborparam);
weird_inds = [];
for i = 1:size(normed_percellfeaturePIE1,1)
    
    if reducedPIE1(i,1) < -120
        weird_inds = [weird_inds;i];
    end
    
end

normed_percellfeaturePIE1_edited = normed_percellfeaturePIE1(setdiff([1:size(normed_percellfeaturePIE1,1)],weird_inds),:);


normed_percellfeature_all = [normed_percellfeature;
            normed_percellfeatureAPX1;
            normed_percellfeatureMOM1;
            normed_percellfeaturePIE1_edited;
            normed_percellfeatureSKN1;
            normed_percellfeatureWWP1];
        
reduced_all = run_umap([normed_percellfeature_all],'n_neighbors',neighborparam);

%%
%plot
WT_color = [.5,.5,.5];
APX1_color = [1,0,0];
MOM1_color = [0,1,0];
PIE1_color = [.5,0,1];
SKN1_color = [1,0,1];
WWP1_color = [0,0,1];

WT_colors = repmat(WT_color,size(normed_percellfeature,1),1);
APX1_colors = repmat(APX1_color,size(normed_percellfeatureAPX1,1),1);
MOM1_colors = repmat(MOM1_color,size(normed_percellfeatureMOM1,1),1);
PIE1_colors = repmat(PIE1_color,size(normed_percellfeaturePIE1_edited,1),1);
SKN1_colors = repmat(SKN1_color,size(normed_percellfeatureSKN1,1),1);
WWP1_colors = repmat(WWP1_color,size(normed_percellfeatureWWP1,1),1);


%%
%getting lineage colors for mutants


%remove weird cells from Pie1 
embeddingStructPIE1_edited = embeddingStructPIE1(setdiff([1:size(normed_percellfeaturePIE1,1)],weird_inds));

%WT colors already exist
[~,APX1_colors,~] = colorArrays(embeddingStructAPX1,celltypes,lineages);
[~,MOM1_colors,~] = colorArrays(embeddingStructMOM1,celltypes,lineages);
[~,PIE1_colors,~] = colorArrays(embeddingStructPIE1_edited,celltypes,lineages);
[~,SKN1_colors,~] = colorArrays(embeddingStructSKN1,celltypes,lineages);
[~,WWP1_colors,~] = colorArrays(embeddingStructWWP1,celltypes,lineages);

%%

all_colors = [WT_colors;
              APX1_colors;
              MOM1_colors;
              PIE1_colors;
              SKN1_colors;
              WWP1_colors];

          
%obtain the indices of each emb
%%
counter = 0;

WT_indices = [counter+1:counter+size(normed_percellfeature,1)];
counter = counter+size(normed_percellfeature,1);

APX1_indices = [counter+1:counter+size(normed_percellfeatureAPX1,1)];
counter = counter+size(normed_percellfeatureAPX1,1);

MOM1_indices = [counter+1:counter+size(normed_percellfeatureMOM1,1)];
counter = counter+size(normed_percellfeatureMOM1,1);

PIE1_indices = [counter+1:counter+size(normed_percellfeaturePIE1_edited,1)];
counter = counter+size(normed_percellfeaturePIE1_edited,1);

SKN1_indices = [counter+1:counter+size(normed_percellfeatureSKN1,1)];
counter = counter+size(normed_percellfeatureSKN1,1);

WWP1_indices = [counter+1:counter+size(normed_percellfeatureWWP1,1)];
counter = counter+size(normed_percellfeatureWWP1,1);

%%          
figure

subplot(2,4,1)
scatter(reduced_all(WT_indices,1),reduced_all(WT_indices,2),20,lineagecolorsWT)
xlabel('UMAP-1')
ylabel('UMAP-2')
title('WT | lineages: ABa (r), ABp (g), C (b), D (f), E (o), M (c), P (p)')

subplot(2,4,5)
scatter(reduced_all(WT_indices,1),reduced_all(WT_indices,2),20,timecolorsWT)
xlabel('UMAP-1')
ylabel('UMAP-2')
title('WT | time: green -> red (21 cells to 350 cells)')

subplot(2,4,2)
scatter(reduced_all(:,1),reduced_all(:,2),20,all_colors(:,:))
xlabel('UMAP-1')
ylabel('UMAP-2')
title('all embryos together')

subplot(2,4,3)
hold on
scatter(reduced_all(WT_indices,1),reduced_all(WT_indices,2),20,all_colors(WT_indices,:))
scatter(reduced_all(APX1_indices,1),reduced_all(APX1_indices,2),20,all_colors(APX1_indices,:))
hold off
xlabel('UMAP-1')
ylabel('UMAP-2')
title('WT (gray) and APX-1')

subplot(2,4,4)
hold on
scatter(reduced_all(WT_indices,1),reduced_all(WT_indices,2),20,all_colors(WT_indices,:))
scatter(reduced_all(MOM1_indices,1),reduced_all(MOM1_indices,2),20,all_colors(MOM1_indices,:))
hold off
xlabel('UMAP-1')
ylabel('UMAP-2')
title('WT (gray) and MOM-2')

subplot(2,4,6)
hold on
scatter(reduced_all(WT_indices,1),reduced_all(WT_indices,2),20,all_colors(WT_indices,:))
scatter(reduced_all(PIE1_indices,1),reduced_all(PIE1_indices,2),20,all_colors(PIE1_indices,:))
hold off
xlabel('UMAP-1')
ylabel('UMAP-2')
title('WT (gray) and PIE-1')

subplot(2,4,7)
hold on
scatter(reduced_all(WT_indices,1),reduced_all(WT_indices,2),20,all_colors(WT_indices,:))
scatter(reduced_all(SKN1_indices,1),reduced_all(SKN1_indices,2),20,all_colors(SKN1_indices,:))
hold off
xlabel('UMAP-1')
ylabel('UMAP-2')
title('WT (gray) and SKN-1')

subplot(2,4,8)
hold on
scatter(reduced_all(WT_indices,1),reduced_all(WT_indices,2),20,all_colors(WT_indices,:))
scatter(reduced_all(WWP1_indices,1),reduced_all(WWP1_indices,2),20,all_colors(WWP1_indices,:))
hold off
xlabel('UMAP-1')
ylabel('UMAP-2')
title('WT (gray) and WWP-1')

sgtitle('codecomposition of WT and 5 homeotic mutants')




%%
figure
scatter(reduced_all(:,1),reduced_all(:,2),20,all_colors(:,:))
xlabel('UMAP-1')
ylabel('UMAP-2')
title('all embryos: WT (gray), APX-1 (red), MOM-2 (green), PIE-1 (purple), SKN-1 (pink), WWP-1 (blue)')

%%
%doing above but without time

        
        
normed_percellfeature_all_NOTIME = normed_percellfeature_all(:,1:53);
        
reduced_all_NOTIME = run_umap([normed_percellfeature_all_NOTIME],'n_neighbors',neighborparam);



figure

subplot(2,4,1)
scatter(reduced_all_NOTIME(WT_indices,1),reduced_all_NOTIME(WT_indices,2),20,lineagecolorsWT)
xlabel('UMAP-1')
ylabel('UMAP-2')
title('lineages: ABa (r), ABp (g), C (b), D (f), E (o), M (c), P (p)')

subplot(2,4,5)
scatter(reduced_all_NOTIME(WT_indices,1),reduced_all_NOTIME(WT_indices,2),20,timecolorsWT)
xlabel('UMAP-1')
ylabel('UMAP-2')
title('time: green -> red (21 cells to 350 cells)')

subplot(2,4,2)
scatter(reduced_all_NOTIME(:,1),reduced_all_NOTIME(:,2),20,all_colors(:,:))
xlabel('UMAP-1')
ylabel('UMAP-2')
title('all embryos together')

subplot(2,4,3)
hold on
scatter(reduced_all_NOTIME(WT_indices,1),reduced_all_NOTIME(WT_indices,2),20,all_colors(WT_indices,:))
scatter(reduced_all_NOTIME(APX1_indices,1),reduced_all_NOTIME(APX1_indices,2),20,all_colors(APX1_indices,:))
hold off
xlabel('UMAP-1')
ylabel('UMAP-2')
title('WT (gray) and APX-1')

subplot(2,4,4)
hold on
scatter(reduced_all_NOTIME(WT_indices,1),reduced_all_NOTIME(WT_indices,2),20,all_colors(WT_indices,:))
scatter(reduced_all_NOTIME(MOM1_indices,1),reduced_all_NOTIME(MOM1_indices,2),20,all_colors(MOM1_indices,:))
hold off
xlabel('UMAP-1')
ylabel('UMAP-2')
title('WT (gray) and MOM-2')

subplot(2,4,6)
hold on
scatter(reduced_all_NOTIME(WT_indices,1),reduced_all_NOTIME(WT_indices,2),20,all_colors(WT_indices,:))
scatter(reduced_all_NOTIME(PIE1_indices,1),reduced_all_NOTIME(PIE1_indices,2),20,all_colors(PIE1_indices,:))
hold off
xlabel('UMAP-1')
ylabel('UMAP-2')
title('WT (gray) and PIE-1')

subplot(2,4,7)
hold on
scatter(reduced_all_NOTIME(WT_indices,1),reduced_all_NOTIME(WT_indices,2),20,all_colors(WT_indices,:))
scatter(reduced_all_NOTIME(SKN1_indices,1),reduced_all_NOTIME(SKN1_indices,2),20,all_colors(SKN1_indices,:))
hold off
xlabel('UMAP-1')
ylabel('UMAP-2')
title('WT (gray) and SKN-1')

subplot(2,4,8)
hold on
scatter(reduced_all_NOTIME(WT_indices,1),reduced_all_NOTIME(WT_indices,2),20,all_colors(WT_indices,:))
scatter(reduced_all_NOTIME(WWP1_indices,1),reduced_all_NOTIME(WWP1_indices,2),20,all_colors(WWP1_indices,:))
hold off
xlabel('UMAP-1')
ylabel('UMAP-2')
title('WT (gray) and WWP-1')

sgtitle('codecomposition of WT and 5 homeotic mutants with NO EXPLICIT TIME FEATURES')

%%
figure
scatter(reduced_all_NOTIME(:,1),reduced_all_NOTIME(:,2),20,all_colors(:,:))
xlabel('UMAP-1')
ylabel('UMAP-2')
title('all embryos: WT (gray), APX-1 (red), MOM-2 (green), PIE-1 (purple), SKN-1 (pink), WWP-1 (blue) NO EXPLICIT TIME FEATURES')



function [reduced1,reduced2] = codecompose(feature1,feature2,neighborparam)

reduced_together = run_umap([feature1;feature2],'n_neighbors',neighborparam);

reduced1 = reduced_together(1:size(feature1,1),:);
reduced2 = reduced_together(size(feature1,1)+1:end,:);

end
function [celltypecolors,lineagecolors,timecolors] = colorArrays(embeddingStruct,celltypes,lineages)

numPointsInEmbedding = size(embeddingStruct,2);
endtime = embeddingStruct(end).captureFrame;


%color scheme

%temporal evolution of embryo
timecolors=zeros(numPointsInEmbedding,3);
for i = 1:numPointsInEmbedding
    hsvcolor = [embeddingStruct(i).captureFrame/endtime,1,1];
    timecolors(i,:)=hsv2rgb(hsvcolor);
end

%cell type 

%%%*** remove this if statement instead make the option whether
%or not you want to have a 3d plot too

celltypecategories = [[0,1,0];[0,0,1];[1,0,0];[1,0,1];[0,1,1]];

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
%note: must be same length as lineages above
lineagecategories = [[1,0,0];[0,1,0];[0,0,1];[1,0,1];[1,.5,0];[0,1,1];[.5,0,1]];

%default color gray
lineagecolors = repmat([.25,.25,.25],numPointsInEmbedding,1);

for i = 1:numPointsInEmbedding
    for j = 1:length(lineages)
        if embeddingStruct(i).lineage == lineages(j)
            lineagecolors(i,:) = lineagecategories(j,:);
        end
    end
end

end
function [percellfeature,normed_percellfeature,embeddingStruct] = pipeline(emb,input_xy_res,anisotropy,tstart,tend)

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

starting_num_cells = 21; %number of cells at which to start computing features--must be one greater than k


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

%%% compute featured %%% 
times = [starttime:sampling:endtime];

params_to_use = lineage_optimized_param_values;

[percellfeature,normed_percellfeature,~] = computeFeatures4(embdat_stabilized,times,feats,params_to_use,normalize_bool,input_xy_res);
disp('feature extraction worked')



end















