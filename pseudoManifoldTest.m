%%

%may 14--testing the behavior of scoring the feature space by using a UMAP
%projection as a proxy for true manifold distance

%results will give feel for how to use it in a formal optimization


%from optimization_WT_mar22_workspace.mat
percellfeature_toUse = lineage_percellfeature_normed;
embeddingStruct_toUse = embeddingStruct;
embryo = embdat_stabilized;

%trying per cell dataset May 25
[percellfeature_toUse,embeddingStruct_toUse] = condenseToOneEntryPerCell(percellfeature_toUse,embeddingStruct_toUse,'median_timepoint_entry');


%number of neighbors looked at to compute scores
number_neighbors = 20;

%UMAP parameter unrelated to scoring--using default 30
neighborparam = 30;

%dimensions to embed space into
dimensions = [55,45,35,25,15,10,5,3,2];


%%%
%compute the scores
scores = zeros(length(dimensions),2); %1st column is fate scores, 2nd is lineage scores
for i = 1:length(dimensions)
    
    d = dimensions(i);
    
    if size(percellfeature_toUse,2) == d
        %if want to look at same-d space as original, do nothing
        lowD_percellfeature = percellfeature_toUse;
    else
        %if want to look at smaller-d space, compute it
        lowD_percellfeature = run_umap(percellfeature_toUse,'n_components',d,'n_neighbors',neighborparam);
    end
    
    %score, log
    [fate_overall,lineage_overall] = computeScore(lowD_percellfeature,embeddingStruct_toUse,embryo,number_neighbors);
    scores(i,:) = [fate_overall,lineage_overall];

end

scores_buffer = scores;

%%

%may 25--testing the behavior of scoring the feature space by using a UMAP
%projection as a proxy for true manifold distance

%results will give feel for how to use it in a formal optimization



%USING PER CELL DATASET, COMPUTESCORE2


%from optimization_WT_mar22_workspace.mat
percellfeature_toUse = lineage_percellfeature_normed;
embeddingStruct_toUse = embeddingStruct;
embryo = embdat_stabilized;

%trying per cell dataset May 25
[percellfeature_toUse,embeddingStruct_toUse] = condenseToOneEntryPerCell(percellfeature_toUse,embeddingStruct_toUse,'median_timepoint_entry');


%number of neighbors looked at to compute scores
number_neighbors = 5;

%UMAP parameter unrelated to scoring--using default 30
neighborparam = 30;

%dimensions to embed space into
dimensions = [55,45,35,25,15,10,5,3,2];

%%%
%compute the scores
scores = zeros(length(dimensions),2); %1st column is fate scores, 2nd is lineage scores
for i = 1:length(dimensions)
    
    d = dimensions(i);
    
    if size(percellfeature_toUse,2) == d
        %if want to look at same-d space as original, do nothing
        lowD_percellfeature = percellfeature_toUse;
    else
        %if want to look at smaller-d space, compute it
        lowD_percellfeature = run_umap(percellfeature_toUse,'n_components',d,'n_neighbors',neighborparam);
    end
    
    %score, log
    [fate_overall,lineage_overall] = computeScore2(lowD_percellfeature,embeddingStruct_toUse,embeddingStruct,embryo,number_neighbors);
    scores(i,:) = [fate_overall,lineage_overall];

end

scores_buffer = scores;

%%
%may 29
%redoing effect of dim reduction using per cell dataset for terminal fate
%and using a time window for lineage optimization

%note: basically same as above except diff scoring function used--computes 
%   the condensed dataset WITHIN it,computes lineage score using a timewindow 
%   (corresponding to the samplerate from the original dataset)


%{

load in the per cell per time dataset

define num neighbors, umap parameter, dimensions

new score function that
    computes the condensed dataset WITHIN it
    computes lineage score using a timewindow (corresponding to the sample
    rate from the original dataset)

%}

%from optimization_WT_mar22_workspace.mat
percellfeature_toUse = lineage_percellfeature_normed;
embeddingStruct_toUse = embeddingStruct;
embryo = embdat_stabilized;

%number of neighbors looked at to compute scores
number_neighbors = 20;

%timewindow forward and back (in minutes) for lineage scores
time_window = 7.5;

%UMAP parameter unrelated to scoring--using default 30
neighborparam = 30;

%dimensions to embed space into
dimensions = [55,45,35,25,15,10,5,4,3,2];



%compute the scores
scores = zeros(length(dimensions),2); %1st column is fate scores, 2nd is lineage scores
for i = 1:length(dimensions)
    
    d = dimensions(i);
    
    if size(percellfeature_toUse,2) == d
        %if want to look at same-d space as original, do nothing
        lowD_percellfeature = percellfeature_toUse;
    else
        %if want to look at smaller-d space, compute it
        lowD_percellfeature = run_umap(percellfeature_toUse,'n_components',d,'n_neighbors',neighborparam);
    end
    
    %make condensed version for fate scoring
    [condensed_lowD_percellfeature,condensed_embeddingStruct] = condenseToOneEntryPerCell(lowD_percellfeature,embeddingStruct_toUse,'median_timepoint_entry');


    
    %score, log
    [fate_overall,lineage_overall]  = computeScore3(lowD_percellfeature,embeddingStruct_toUse,embryo,number_neighbors,time_window,condensed_lowD_percellfeature,condensed_embeddingStruct,starttime,sampling,endtime);
    scores(i,:) = [fate_overall,lineage_overall];

end

scores_buffer = scores;




%other code


%%


%make labels for graph
names_chars = {'55 (original)';
               '45';
               '35';
               '25';
               '15';
               '10';
               '5';
               '4';
               '3';
               '2'};

%names_chars = {'55','45','35','25','15','10','5'}';

names_cat = categorical(names_chars);
names_cat = reordercats(names_cat,names_chars);


b = bar(names_cat,scores);

xlabel('number of dimensions')

tstr = ['percent of ',num2str(number_neighbors),' nns of same fate (b), same lineage (r) as measured in an X-d UMAP reduction of the original feature space (using lineage optimized paramters)'];
tstr = [tstr newline 'Lineage calculated with time window = ',num2str(time_window),' mins on per cell per time dataset, fate calculated on per cell dataset'];
title(tstr)




%%

%may 25 embedding space flow between WT, mutant (per cell space, not per
%cell per time)

neighborparam = 30;

%load wt, mutant, make condensed versions
%loaded from percellfeatureWT_mar25.mat and percellfeatureAPX1_mar25.mat,
%etc...

wt_percellfeature = percellfeature;%RAW
wt_embeddingStruct = embeddingStruct;%RAW

%APX-1
%{
mut_percellfeature = percellfeatureAPX1;%RAW
mut_embeddingStruct = embeddingStructAPX1;%RAW
%}

%MOM-2
%{
mut_percellfeature = percellfeatureMOM1;%RAW
mut_embeddingStruct = embeddingStructMOM1;%RAW
%}

%PIE-1
%{
mut_percellfeature = percellfeaturePIE1;%RAW
mut_embeddingStruct = embeddingStructPIE1;%RAW
%}

%SKN-1
%{
mut_percellfeature = percellfeatureSKN1;%RAW
mut_embeddingStruct = embeddingStructSKN1;%RAW
%}

%WWP-1
%
mut_percellfeature = percellfeatureWWP1;%RAW
mut_embeddingStruct = embeddingStructWWP1;%RAW
%}



%truncate wt dataset to to same length as mutant

indices_list_55 = load('indices_list_55.mat');
indices_list_55 = indices_list_55.indices_list_55;

max_num_cells = mut_embeddingStruct(end).numCellsAtCapture;
for i = 1:length(wt_embeddingStruct)
    if wt_embeddingStruct(i).numCellsAtCapture > max_num_cells
        break
    end
end
wt_truncated = wt_percellfeature(1:i,:);
wt_embeddingStruct_truncated = wt_embeddingStruct(1:i);

%normalize
wt_normed = quick_normalize(wt_truncated,indices_list_55);
mut_normed = quick_normalize(mut_percellfeature,indices_list_55);

[wt_percellfeature_condensed,wt_embeddingStruct_condensed] = condenseToOneEntryPerCell(wt_normed,wt_embeddingStruct_truncated,'median_timepoint_entry');
[mut_percellfeature_condensed,mut_embeddingStruct_condensed] = condenseToOneEntryPerCell(mut_normed,mut_embeddingStruct,'median_timepoint_entry');


%reduce -- AFTER CONDENSING
both = [wt_percellfeature_condensed;mut_percellfeature_condensed];
both_reduced  = run_umap(both,'n_neighbors',neighborparam);

wt_percellfeature_reduced = both_reduced(1:size(wt_percellfeature_condensed,1),:);
mut_percellfeature_reduced = both_reduced(size(wt_percellfeature_condensed,1)+1:end,:);


%make quivers

names_list_wt = ["seed"];
for i = 1:length(wt_embeddingStruct_condensed)
    names_list_wt(i) = convertCharsToStrings(wt_embeddingStruct_condensed(i).name);
end
names_list_mut = ["seed"];
for i = 1:length(mut_embeddingStruct_condensed)
    names_list_mut(i) = convertCharsToStrings(mut_embeddingStruct_condensed(i).name);
end

quivs = [];
quiv_points = [];
wt_inds_corresponding_to_quivers = [];
c = 1;
for i = 1:length(wt_embeddingStruct_condensed)
    
    dex = find(names_list_mut == wt_embeddingStruct_condensed(i).name);
    if ~isempty(dex)
        
        mut_point = [mut_percellfeature_reduced(dex,1),mut_percellfeature_reduced(dex,2)];
        wt_point = [wt_percellfeature_reduced(i,1),wt_percellfeature_reduced(i,2)];
        quivs(c,:) = [mut_point(1)-wt_point(1),mut_point(2)-wt_point(2)];
        quiv_points(c,:) = [wt_point(1),wt_point(2)];
        
        %log which cell in the wt has which quiver (i.e. quiver(1)
        %originates from cell x in the wt
        wt_inds_corresponding_to_quivers = [wt_inds_corresponding_to_quivers;i];
        c = c+1;

    end
    
end


%%
%plotting


%single plot version 
%{
%colors

celltypes = ["pharynx","neuron","bodywallmuscle","hypoderm","gut"];
lineages = ["ABa","ABp","C","D","E","MS","P"];

[wt_celltypecolors,wt_lineagecolors,wt_timecolors] = colorArrays(wt_embeddingStruct_condensed,celltypes,lineages);
[mut_celltypecolors,mut_lineagecolors,mut_timecolors] = colorArrays(mut_embeddingStruct_condensed,celltypes,lineages);



figure
hold on
scatter(wt_percellfeature_reduced(:,1),wt_percellfeature_reduced(:,2),20,wt_lineagecolors)
scatter(mut_percellfeature_reduced(:,1),mut_percellfeature_reduced(:,2),10,mut_lineagecolors,'filled')

quiver(quiv_points(:,1),quiv_points(:,2),quivs(:,1),quivs(:,2),'k');
title(['per cell wt, mut, reduction *after* removing repeated cell entries' newline 'ABa (r), ABp (g), C (b), D (f), E (o), M (c), P (p)' newline 'WT (circles), APX1 (filled)'])
%}



%multi plot fromm above


%define celltypes, lineages
celltypes = ["pharynx","neuron","bodywallmuscle","hypoderm","gut"];
lineages = ["ABa","ABp","C","D","E","MS","P"];

%generate colors
[wt_celltypecolors,wt_lineagecolors,wt_timecolors] = colorArrays(wt_embeddingStruct_condensed,celltypes,lineages);
[mut_celltypecolors,mut_lineagecolors,mut_timecolors] = colorArrays(mut_embeddingStruct_condensed,celltypes,lineages);

figure

%1: all lineages colored, all arrows shown scaled
subplot(2,4,1)
hold on
scatter(wt_percellfeature_reduced(:,1),wt_percellfeature_reduced(:,2),20,wt_lineagecolors)
scatter(mut_percellfeature_reduced(:,1),mut_percellfeature_reduced(:,2),10,mut_lineagecolors,'filled')
quiver(quiv_points(:,1),quiv_points(:,2),quivs(:,1),quivs(:,2),'k');
hold off

title('ABa (r), ABp (g), C (b), D (f), E (o), M (c), P (p)')
xlabel('UMAP-1')
ylabel('UMAP-2')

%5 colored by time, no arrows
subplot(2,4,5)
hold on
scatter(wt_percellfeature_reduced(:,1),wt_percellfeature_reduced(:,2),20,wt_timecolors)
scatter(mut_percellfeature_reduced(:,1),mut_percellfeature_reduced(:,2),10,mut_timecolors,'filled')
hold off

title('time (green->red)')
xlabel('UMAP-1')
ylabel('UMAP-2')

%define subplots on which to put 6 lineages that are not P
plots_inds = [2,3,4,6,7,8];

%do the lineages
for i = 1:6
    
    %create new colors where all cells not a particular lineage are gray
    wt_colors_new = wt_lineagecolors;
    wt_inds_correct_lineage = [];%for selecting relevant quivers
    for j = 1:length(wt_embeddingStruct_condensed)
        if wt_embeddingStruct_condensed(j).lineage == lineages(i)
            wt_inds_correct_lineage = [wt_inds_correct_lineage;j];
        else
            wt_colors_new(j,:) = [.25,.25,.25];
        end
    end
    mut_colors_new = mut_lineagecolors;
    for j = 1:length(mut_embeddingStruct_condensed)
        if mut_embeddingStruct_condensed(j).lineage ~= lineages(i)
            mut_colors_new(j,:) = [.25,.25,.25];
        end
    end
    
    %get relevant quivers
    relevant_quiv_inds = [];
    for j = 1:size(quivs,1)
        if ismember(wt_inds_corresponding_to_quivers(j),wt_inds_correct_lineage)
            relevant_quiv_inds = [relevant_quiv_inds;j];   
        end
    end
    
    
    %June 1 just keeping all colors
    wt_colors_new = wt_lineagecolors;
    mut_colors_new = mut_lineagecolors;
    
    
    subplot(2,4,plots_inds(i))
    hold on
    scatter(wt_percellfeature_reduced(:,1),wt_percellfeature_reduced(:,2),20,wt_colors_new)
    scatter(mut_percellfeature_reduced(:,1),mut_percellfeature_reduced(:,2),10,mut_colors_new,'filled')
    quiver(quiv_points(relevant_quiv_inds,1),quiv_points(relevant_quiv_inds,2),quivs(relevant_quiv_inds,1),quivs(relevant_quiv_inds,2),'k','AutoScale','off','ShowArrowHead','off');  
    hold off
    
    title(lineages(i))
    xlabel('UMAP-1')
    ylabel('UMAP-2')
    
end

%{
%APX-1
sgtitle(['lineage flows between WT and APX-1 (ABp->ABa)' newline 'one point per cell' newline 'lines connect WT (circles) cell to corresponding cell in mutant (dots)'])
savefig('APX-1_flows_June01.fig')
%}


%{
%MOM-2
sgtitle(['lineage flows between WT and MOM-2 (ABar->ABal,E->MS)' newline 'one point per cell' newline 'lines connect WT (circles) cell to corresponding cell in mutant (dots)'])
savefig('MOM-2_flows_June01.fig')
%}

%{
%PIE-1
sgtitle(['lineage flows between WT and PIE-1 (ABp->ABa,P2->EMS)' newline 'one point per cell' newline 'lines connect WT (circles) cell to corresponding cell in mutant (dots)'])
savefig('PIE-1_flows_June01.fig')
%}

%{
%SKN-1
sgtitle(['lineage flows between WT and SKN-1 (EMS->C, ABpra->ABpla,ABara->ABala,ABalp->ABarp )' newline 'one point per cell' newline 'lines connect WT (circles) cell to corresponding cell in mutant (dots)'])
savefig('SKN-1_flows_June01.fig')
%}

%
%WWP-1
sgtitle(['lineage flows between WT and WWP-1 (ABa->ABp,ABpla->ABpra)' newline 'one point per cell' newline 'lines connect WT (circles) cell to corresponding cell in mutant (dots)'])
savefig('WWP-1_flows_June01.fig')
%}



%%
function [fate_overall, lineage_overall] = computeScore3(adjusted_percellfeature,adjusted_embeddingStruct,embryo,number_neighbors,time_window,condensed_percellfeature,condensed_embeddingStruct,starttime,sampling,endtime)

    %compute scores where lineage has a time window (on per cell per time) and fate is on the
    %condensed dataset
    
    
    %compute time window for lineage (taken from computefeatures4--see for
    %more details)
    
    %value is number of minutes forward and back to look

    %   parameter value t becomes fraction of time from 4 cell to twitch ~~ 400 mins
    %   twitch may not be measurable in dataset, so approx using 350cell tpt

    %obtain first time at which emb is 4 cells
    fourcelltime = 1;
    while size(embryo(fourcelltime).finalpoints,1) < 4
        fourcelltime = fourcelltime + 1;
    end

    %obtain first time at which emb is 350 cells
    threefiftycelltime = 1;
    while size(embryo(threefiftycelltime).finalpoints,1) < 350
        threefiftycelltime = threefiftycelltime + 1;
    end

    %approximate twitch timepoint
    %   p = approximate ratio of time from 4cell to twitch over 4cell to 350 cell = (475-75)/(345-75)  (
    %   from https://www.wormatlas.org/hermaphrodite/introduction/mainframe.html)
    p = (475-75)/(345-75);
    approx_twitch_time = p*(threefiftycelltime - fourcelltime) + fourcelltime;

    %convert minutes to corresponding numbers of timesteps in embryo
    %dataset
    time_window = ceil(time_window/400 * (approx_twitch_time-fourcelltime));
   
   
    %LINEAGE CONTINUITY
    lineage_percellscores = [];
    
    for j = 1:length(adjusted_percellfeature)


        %score is  percentage of nns that are same cell, progenitor, or
        %daughter

        %find neighbors in feature space
        [~,inds] = pdist2(adjusted_percellfeature,adjusted_percellfeature(j,:),'euclidean','Smallest',number_neighbors+1);
        inds = inds(2:end);
        
        
        %find the embryo coords of cells that are in continuous lineage
        %w/ cell j
        
        continuous_cells = []; %embryo coords, first column embryo time, second column cell index
        
        
        start_time = adjusted_embeddingStruct(j).captureFrame;
        start_ind = adjusted_embeddingStruct(j).captureIndex;

        
        %obtain past cells  
        [~,~,~,~,all_prev_cell_times,all_prev_cell_inds] = backTraverse_keepall(start_time,start_ind,'timesteps',time_window,embryo);%change to an amount of time
        continuous_cells = [continuous_cells;all_prev_cell_times,all_prev_cell_inds];
        
        
        %same process; make a new function that returns all future cell
        %times, inds, map them to coords in feature space, tally
        
        %obtain future cells
        [~,~,all_later_cell_times,all_later_cell_inds] = forwardTraverse_keepall(start_time,start_ind,'timesteps',time_window,embryo);%change to an amount of time
        continuous_cells = [continuous_cells;all_later_cell_times,all_later_cell_inds];
        
        
        %convert included cells to embeddingStruct (featurespace) indices
        
        [embdatToEmbeddingStruct,~] = mappingMaker(embryo,starttime,sampling,endtime);
        
        continuous_cells_featureCoords = zeros(size(continuous_cells,1),1);
        for l = 1:size(continuous_cells,1)
            
            continuous_cells_featureCoords(l) = embdatToEmbeddingStruct(continuous_cells(l,1)).map(continuous_cells(l,2));

        end
        
        
        %tally how many neighbors are 'continuous' (in included cells)
        num_continuous = 0;
        for k = 1:length(inds)
   
            if ismember(inds(k),continuous_cells_featureCoords)
                num_continuous = num_continuous+1;
            end

        end

        %find percent continuous, add to lineage_percellscores
        percent_continuous = num_continuous/number_neighbors;
        lineage_percellscores = [lineage_percellscores;percent_continuous];

    end
    
    
    
    
    %TERMINAL FATE SIMILARITY
    
    fate_percellscores = [];
    for j = 1:length(condensed_percellfeature)
        
        %find neighbors in feature space
        [~,condensed_inds] = pdist2(condensed_percellfeature,condensed_percellfeature(j,:),'euclidean','Smallest',number_neighbors+1);
        condensed_inds = condensed_inds(2:end);

        
        if condensed_embeddingStruct(j).celltype ~= "OTHER"
            %score is  percentage of nns that are same type as cell in
            %question

            %tally how many neighbors are of same type
            num_same_type = 0;
            for k = 1:number_neighbors

                if condensed_embeddingStruct(condensed_inds(k)).celltype == condensed_embeddingStruct(j).celltype
                    num_same_type = num_same_type+1;
                end

            end
            %percent of neighbors that are the same type
            percent_same_type = num_same_type/number_neighbors;

            %log score
            fate_percellscores = [fate_percellscores;percent_same_type];

        end
        
    end
    
    
    %average (overall) scores
    fate_overall = mean(fate_percellscores);
    lineage_overall = mean(lineage_percellscores);

end

function [fate_overall, lineage_overall] = computeScore2(adjusted_percellfeature,adjusted_embeddingStruct,refStruct,embryo,number_neighbors)

    %clone of computescore that is able to compute lineage scores on per cell dataset
    
    lineage_percellscores = [];
    fate_percellscores = [];

    
    %make array of names showing lineage continuity relationships (using
    %a full embeddingStruct as a reference)
    %   1st column is names of cells
    %   2nd column is names of progenitors
    %   3rd,4th columns are names of daughters
    %   "X" indicates that progenitor/daughter doesn't exist 
    
    
    lineage_array = ["seed","seed","seed","seed"];
    c = 1;

    for i = 1:length(refStruct)
        
        %if current cell hasn't already been catalogued 
        cells_name = convertCharsToStrings(refStruct(i).name);
        if ~ismember(cells_name,lineage_array(:,1))
            
            %add cell's own name
            lineage_array(c,1) = cells_name;
            
            %obtain name of progenitor, add
            start_time = refStruct(i).captureFrame;
            start_ind = refStruct(i).captureIndex;

            [prog_time_pt,prog_ind_pt,~,~] = backTraverse(start_time,start_ind,'generations',1,embryo);
            prog_name = convertCharsToStrings(embryo(prog_time_pt).cellnames(prog_ind_pt));
            
            if prog_name ~= cells_name
                lineage_array(c,2) = prog_name;
            else
                lineage_array(c,2) = "X";
            end
            
            %obtain names of daughter(s), add

            [daughter_inds,daughter_time_pt] = forwardTraverse(start_time,start_ind,'generations',1,embryo);
            daughter_names = ["seed"];
            d = 1;
            for l = 1:length(daughter_inds)
                daughter_names(d) = convertCharsToStrings(embryo(daughter_time_pt).cellnames{daughter_inds(l)});
                d = d + 1;
            end
            
            daughter_names = setdiff(daughter_names,"seed");
            switch length(daughter_names)
                case 0
                    lineage_array(c,3:4) = ["X","X"];
                case 1
                    if daughter_names(1) ~= cells_name
                        lineage_array(c,3:4) = [daughter_names(1),"X"];
                    else
                        lineage_array(c,3:4) = ["X","X"];
                    end
                case 2
                    lineage_array(c,3:4) = [daughter_names(1),daughter_names(2)];
                    
            end
            c = c+1;
            
        end
        
    end

    
    for j = 1:length(adjusted_percellfeature)

        %LINEAGE CONTINUITY

        %score is  percentage of nns that are same cell, progenitor, or
        %daughter

        %find neighbors
        [~,inds] = pdist2(adjusted_percellfeature,adjusted_percellfeature(j,:),'euclidean','Smallest',number_neighbors+1);
        inds = inds(2:end);

        %tally how many neighbors are 'continuous'
        num_continuous = 0;
        correct_names = lineage_array(find(lineage_array(:,1)==adjusted_embeddingStruct(j).name),:);
        
        for k = 1:number_neighbors

            if ~isempty(find(correct_names == adjusted_embeddingStruct(inds(k)).name))
                num_continuous = num_continuous+1;
            end

        end

        %find percent continuous, add to lineage_percellscores
        percent_continuous = num_continuous/number_neighbors;
        lineage_percellscores = [lineage_percellscores;percent_continuous];


        %TERMINAL FATE SIMILARITY
        if adjusted_embeddingStruct(j).celltype ~= "OTHER"
            %score is  percentage of nns that are same type as cell in
            %question

            %tally how many neighbors are of same type
            num_same_type = 0;
            for k = 1:number_neighbors

                if adjusted_embeddingStruct(inds(k)).celltype == adjusted_embeddingStruct(j).celltype
                    num_same_type = num_same_type+1;
                end

            end
            %percent of neighbors that are the same type
            percent_same_type = num_same_type/number_neighbors;

            %log score
            fate_percellscores = [fate_percellscores;percent_same_type];

        end

    end

    %average (overall) scores
    fate_overall = mean(fate_percellscores);
    lineage_overall = mean(lineage_percellscores);

end

function [fate_overall, lineage_overall] = computeScore(adjusted_percellfeature,adjusted_embeddingStruct,embryo,number_neighbors)

    %compute scores [taken directly from hyperparameter tuning]
    
    lineage_percellscores = [];
    fate_percellscores = [];

    for j = 1:length(adjusted_percellfeature)

        %LINEAGE CONTINUITY

        %score is  percentage of nns that are same cell, progenitor, or
        %daughter

        %find neighbors
        [~,inds] = pdist2(adjusted_percellfeature,adjusted_percellfeature(j,:),'euclidean','Smallest',number_neighbors+1);
        inds = inds(2:end);

        %create list of names that qualify as same cell, progenitor, or
        %daughter

        %add same cell
        correct_names = [convertCharsToStrings(adjusted_embeddingStruct(j).name)];

        start_time = adjusted_embeddingStruct(j).captureFrame;
        start_ind = adjusted_embeddingStruct(j).captureIndex;

        %obtain name of progenitor

        [prog_time_pt,prog_ind_pt,~,~] = backTraverse(start_time,start_ind,'generations',1,embryo);
        prog_name = convertCharsToStrings(embryo(prog_time_pt).cellnames(prog_ind_pt));

        correct_names = [correct_names;prog_name];

        %obtain names of daughters
        [daughter_inds,daughter_time_pt] = forwardTraverse(start_time,start_ind,'generations',1,embryo);

        for l = 1:length(daughter_inds)
            correct_names = [correct_names;convertCharsToStrings(embryo(daughter_time_pt).cellnames{daughter_inds(l)})];
        end


        %tally how many neighbors are 'continuous'
        num_continuous = 0;
        for k = 1:number_neighbors

            if ~isempty(find(correct_names == adjusted_embeddingStruct(inds(k)).name))
                num_continuous = num_continuous+1;
            end

        end

        %find percent continuous, add to lineage_percellscores
        percent_continuous = num_continuous/number_neighbors;
        lineage_percellscores = [lineage_percellscores;percent_continuous];


        %TERMINAL FATE SIMILARITY
        if adjusted_embeddingStruct(j).celltype ~= "OTHER"
            %score is  percentage of nns that are same type as cell in
            %question

            %tally how many neighbors are of same type
            num_same_type = 0;
            for k = 1:number_neighbors

                if adjusted_embeddingStruct(inds(k)).celltype == adjusted_embeddingStruct(j).celltype
                    num_same_type = num_same_type+1;
                end

            end
            %percent of neighbors that are the same type
            percent_same_type = num_same_type/number_neighbors;

            %log score
            fate_percellscores = [fate_percellscores;percent_same_type];

        end

    end

    %average (overall) scores
    fate_overall = mean(fate_percellscores);
    lineage_overall = mean(lineage_percellscores);

end

function [condensed_feature,condensed_embeddingStruct] = condenseToOneEntryPerCell(feature,embeddingStruct,mode)
    %take in a feature space/embeddingStruct where each cell is
    %represented and multiple times and return versions of both where each
    %cell is represented only once
    
    
    condensed_feature = [];
	condensed_embeddingStruct = embeddingStruct(1);%seed with first entry to avoid nonlike struct error
    c = 1;
    names_list = ["seed"];
    for i = 1:size(feature,1)
        names_list(i) = convertCharsToStrings(embeddingStruct(i).name);
    end
    queue = [1:size(feature,1)];
    
    while ~isempty(queue)
        %get name of first cell in there
        top_name = embeddingStruct(queue(1)).name;
        
        
        %get inds of all entries of that cell
        same_cell_inds = find(names_list == top_name);
        
        %make embStruct with all entries for that particular cell
        same_cellStruct = embeddingStruct(same_cell_inds);
        
        switch mode
            case 'average_of_all_entries' %conceptually doesn't make sense lol
                
                %update condensed_feature and embeddingStruct
                feat_entry = mean(feature(same_cell_inds,:),1);
                condensed_feature(c,:) = feat_entry;
                
                condensed_embeddingStruct(c) = embeddingStruct(queue(1));
                c = c+1;
                
            case 'median_timepoint_entry'
                times = [embeddingStruct(same_cell_inds).captureFrame];%correspond to same_cell_inds
                
                %find the 'median' timept entry
                sorted_times = sort(times);
                med_tpt = sorted_times(ceil(length(sorted_times)/2));%ceil ensures that it works even w/one time
                entry_to_use = same_cell_inds(find(times==med_tpt));
                
                %update condensed_feature and embeddingStruct
                feat_entry = feature(entry_to_use,:);
                condensed_feature(c,:) = feat_entry;
                condensed_embeddingStruct(c) = embeddingStruct(entry_to_use);
                c = c+1;
        end
        %remove indices from queue
        queue = setdiff(queue,same_cell_inds);
    end
    
    %remove now meaningless fields from embeddingStruct
    
    
    
    %condensed_embeddingStruct = rmfield(condensed_embeddingStruct,{'captureFrame','captureIndex'});
    
    
    
    %note: numCellsAtCapture and realSpaceCoord arent fully meaningful or
    %accurate if 'average_of_all_entries' option is used, but may be a
    %rough reference
    

    
    %{
    queue of indices to visit = all indices
    
    while queue not empty
        see what cell name it is
        find all entries with same cell name
        organize them by time
        either take avg median entry
        append to vars condensed feature, condensed_embeddingStruct
        remove indices corresponding to that cell name from queue
        
    
        reason through whether indices (from find) and queue will work
            
    %}
  
end 

function [normed_featurevector] = quick_normalize(featurevector,indices_list)
    %* scales time features from min to max within the dataset, not using
    %some external landmark like 4 cells stage--this means any two feature
    %vectors being compared using this function need to span the same
    %developmental window (in terms of num cells in embryo)


    %*assumes that all feats are present and ordered as in:
    %feature list
    %{

    1     "beginning_of_current_cell_cycle"
    2     "end_of_current_cell_cycle"
    3     "duration_of_current_cell_cycle"
    4     "time_until_next_division"
    5     "time_since_last_division"
    6     "fraction_of_current_cycle_elapsed"
    7     "similarity_to_ancestors_cell_cycle"
    8     "similarity_to_cousins_cell_cycle"
    9     "first_deriv_x"
    10    "first_deriv_y"
    11    "first_deriv_z"
    12    "laplacian"
    13    "local_normal_vector"
    14    "residual"
    15    "distances_to_k_nns"
    16    "signed_distances"
    17    "vec_to_position_x_mins_ago"
    18    "vec_to_position_x_mins_in_future"
    19    "vec_to_site_of_last_division"
    20    "vec_to_site_of_next_division"
    21    "displacement_vec_over_cell_cycle"
    22    "proportion_of_present_nns_that_are_descendants_of_past_nns"
    23    "filtered_nuc_size"
    24    "number_of_cell_in_embryo"
    25    "developmental_time"

    %}


    %normalizations corresponding to feats
    ordered_norms = ["scale_4cell_to_end"; 
                  "scale_4cell_to_end";       
                  "scale_4cell_to_end";  
                  "scale_4cell_to_end";        
                  "scale_4cell_to_end";        
                  "nothing";
                  "nothing";          
                  "nothing";            
                  "column_zscore";                    
                  "column_zscore";                    
                  "column_zscore";                    
                  "column_zscore";                        
                  "divide_avg_mag";              
                  "column_zscore";                         
                  "all_zscore"                  
                  "divide_avg_mag"                  
                  "divide_avg_mag";                             
                  "divide_avg_mag";                             
                  "divide_avg_mag";                            
                  "divide_avg_mag";                             
                  "divide_avg_mag";                             
                  "nothing"; 
                  "column_zscore";                
                  "scale_0_to_max_num_cells";         
                  "scale_4cell_to_end"];   

    %make copy
    normed_featurevector = featurevector;

    for i = 1:size(indices_list,1)

        %get piece of normed_featurevector
        entries = normed_featurevector(:,indices_list(i,1):indices_list(i,2));

        if ordered_norms(i) == "scale_4cell_to_end"
            %linearly scale all entries on [0,1] (using range of values in dataset, not external 4 cell stage)
            entries = rescale(entries);

        elseif ordered_norms(i) == "scale_0_to_max_num_cells"
            %linearly scale all entries on [0,1] (using range of values in dataset, not external 4 cell stage)
            entries = rescale(entries);

        elseif ordered_norms(i) == "column_zscore"
            %take zscore of column ( indices_list(i,1)==indices_list(i,2) )
            entries = zscore(entries);

        elseif ordered_norms(i) == "all_zscore"
            %take zscore of all elements
            entries = zscore(entries,0,'all');

        elseif ordered_norms(i) == "divide_avg_mag"
            %divide array elements by avg magitude of a vector in the array
            entries = avgMagnitude(entries);
        end

        %put newly normalized entries back in normed_featurevector
        normed_featurevector(:,indices_list(i,1):indices_list(i,2)) = entries;

    end

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

function [embdatToEmbeddingStruct,embeddingStructToEmbdat] = mappingMaker(embdat,starttime,sampling,endtime)
%Makes struct mappings between embdat format data and embeddingStruct-data 

embdatToEmbeddingStruct = struct();
embeddingStructToEmbdat = struct();


times = [starttime:sampling:endtime];

c = 1;
for time_index  = 1:length(times)
    
    t = times(time_index);
    
   	map_temp = zeros(size(embdat(t).finalpoints,1),1);
    for i = 1:size(embdat(t).finalpoints,1)
        
        map_temp(i) = c;
        
        %embeddingStruct->embdat
        embeddingStructToEmbdat(c).embryoTime = t;
        embeddingStructToEmbdat(c).embryoInd = i;
 
        %increment
        c = c+1;
    end
    %embdat->embeddingStruct
    embdatToEmbeddingStruct(t).map = map_temp;
    
    
end

%populate empty embdatToEmbeddingStruct entries with arrays so that they
%can be indexed into even if they arent present in the embeddingStruct

times_notIncluded = setdiff([1:length(embdat)],times);

for i = 1:length(times_notIncluded)
    
    t = times_notIncluded(i);
    embdatToEmbeddingStruct(t).map = zeros(size(embdat(t).finalpoints,1),1);
    
end

end












