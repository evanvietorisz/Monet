% analyzing clustering results

%%
%prep

%raw data, parameters
result_struct = importdata('Nov6_louvain_dim55_clus.txt');
clustermatrix = result_struct.data;

%load struct representation
%   1st index is the clustering experiment
%   2nd index is the cluster within that experiment
%   contains the names of the cells in that cluster
%other file: load('clustering_results_sep30.mat')

%dataset
emb = load('dispim_emb2_full.mat');
featurespace = emb.normed_percellfeature;
embStruct = emb.embeddingStruct;

%minimum length of a continuous sublineage (# generations)
min_length = 3;


%%
%aux

%compile the names
all_names = cell(1);
for i = 1:length(emb.embeddingStruct)
    
    all_names{i} = char(emb.embeddingStruct(i).name);

end
all_names = all_names';
%%
%get the umap
reduced = run_umap(featurespace,'n_components',15,'n_neighbors',30);


%%

%get terminal fate composition of all cells in embryo

fates_to_look_for = ["pharynx";
                     "neuron";
                     "hypoderm";
                     "bodywallmuscle";
                     "gut";
                     "OTHER"];

[terminal_fates] = fateOfTerminalChildren2(emb,fates_to_look_for);

%normalize to percent composition of each term fate
terminal_composition = terminal_fates;
for i = 1:size(terminal_composition,1)
    if sum(terminal_composition(i,:)) ~= 0 %if all zeros (programmed cell death), leave as all zeros
        terminal_composition(i,:) = terminal_composition(i,:)/sum(terminal_composition(i,:));
    end
end





%%

for i = 1:length(seven_inds)
    
    if startsWith(embStruct(seven_inds(i)).name,'E')
        
        
        if embStruct(seven_inds(i)).celltype ~= "gut"
        
            disp(seven_inds(i))
            disp(embStruct(seven_inds(i)).name)
            disp(embStruct(seven_inds(i)).celltype)

        end
    end
end


%%

%%% pie charts %%%

%generate the report

for experiment_num = 1:size(clustermatrix,2)
    
    
    fignum = 1;
    figure(fignum)
    
    for cluster_num = 1:max(clustermatrix(:,experiment_num))
        
        %isolate cells contained in the cluster
        cell_inds = clustermatrix(:,experiment_num) == cluster_num;
        %cell_inds = find(cell_inds == 1);
        
        %%% content for pie charts %%%
        
        %tissue
        tissues = [embStruct(cell_inds).celltype];
        [tissue_counts, tissue_groups] = groupcounts(tissues');
        
        %top-level lineage
        lins = [embStruct(cell_inds).lineage];
        [lin_counts, lin_groups] = groupcounts(lins');
        
        
        %progenitor type
        
        comp_terminal_children = terminal_composition(cell_inds,:);

        

        
        
        %%% content for scatter %%%
        
        %time the cells come from
        num_cells_at_capture = [embStruct(cell_inds).numCellsAtCapture]';
        
        %bounds
        lower_bound = 1;
        upper_bound = max([embStruct.numCellsAtCapture]);
        
        %assign top level lineage to values
        scatter_yvals_lineage = zeros(length(find(cell_inds==1)),1);
        colors = repmat([.5,.5,.5],length(find(cell_inds==1)),1);
        %colors = zeros(length(find(cell_inds==1)),3);
        for i = 1:length(find(cell_inds==1))
            
            switch lins(i)
                case 'ABa'
                    scatter_yvals_lineage(i) = 2;
                case 'ABp'
                    scatter_yvals_lineage(i) = 3;
                case 'MS'
                    scatter_yvals_lineage(i) = 4;
                case 'C'
                    scatter_yvals_lineage(i) = 5;
                case 'D'
                    scatter_yvals_lineage(i) = 6;
                case 'E'
                    scatter_yvals_lineage(i) = 7;
                otherwise
                    scatter_yvals_lineage(i) = 1;
            end
            
            %assign colors to committed cells
            
            if sum(comp_terminal_children(i,:)~=0) == 1 %if cell's terminal children are all one type
                
                switch find(comp_terminal_children(i,:)~=0)
                    case 1 %pharynx
                        colors(i,:) = [0,1,0];
                    case 2 %neuron
                        colors(i,:) = [0,0,1];
                    case 3 %hypoderm
                        colors(i,:) = [1,0,1];
                    case 4 %bodywallmuscle
                        colors(i,:) = [1,0,0];
                    case 5 %gut
                        colors(i,:) = [0,1,1];
                    %note: 'OTHER' will just stay as gray. uncommitted
                    %cells will remain gray   
                end
                        
            end
            

            
        end
        


        %%% plot %%%
        
        switch mod(cluster_num,4)
            case 1
                subplot_row = 1;
            case 2
                subplot_row = 2;
            case 3
                subplot_row = 3;
            case 0
                subplot_row = 4;
        end
        
        
        base_subplot_num = (subplot_row-1)*5;
        
        % 1. fate pie chart
        
        vals = tissue_counts';
        labs = cellstr(tissue_groups');
        titl = 'fate identity of cells';
        subplot_column_num = 1;
        %ylabel(['cluster ',num2str(cluster_num)])

        
        %repeated
        subplot(4,5,base_subplot_num + subplot_column_num)
        p = pie(vals,labs);
        
        %{
        h = findobj(p,'Type','Text'); % text handles
        isSmall = startsWith({h.String}, '<');
        set(h(isSmall),'String', '')
        %}
        title(titl)
        
        tt = text(-1.75,.3,['cluster ',num2str(cluster_num)],'HorizontalAlignment','center','VerticalAlignment','middle');
        set(tt,'Rotation',90)
        set(tt,'FontSize',14)
        %
        
        % 2. lineage pie chart
        
        vals = lin_counts';
        labs = cellstr(lin_groups');
        titl = 'lineage identity of cells';
        subplot_column_num = 2;
        
        %
        subplot(4,5,base_subplot_num + subplot_column_num)
        pie(vals,labs)
        title(titl)
        %
        
        

        
        % 3. fate of terminal children pie chart
        
        fate_counts = sum(comp_terminal_children,1)';
        vals = fate_counts(fate_counts~=0);
        labs = fates_to_look_for(fate_counts~=0);
        
        titl = 'fate of terminal children';
        subplot_column_num = 3;
        
        
        %
        subplot(4,5,base_subplot_num + subplot_column_num)
        pie(vals,labs)
        title(titl)
        %

        
        
        % 4. scatter
        
        subplot_column_num = 4;
        subplot(4,5,base_subplot_num + subplot_column_num)
        
        scatter(num_cells_at_capture,scatter_yvals_lineage,15,colors)
        xlim([lower_bound,upper_bound])
        ylim([0,8])
        yticks([1,2,3,4,5,6,7])
        yticklabels({'else','ABa','ABp','MS','C','D','E'})
        xlabel('capture frame')
        ylabel('top-level lineage')
        %title(['top-level lineage vs capture frame | color indicates fate commitment' newline 'pharynx (g), neuron (b), hypoderm (p), bodymuscle (r), gut (c)'])
        title('timeline | color indicates fate commitment')

        
        
        % 5. umap
        
        umap_colors = zeros(size(cell_inds,1),3);
        for i = 1:size(cell_inds,1)
            if cell_inds(i) == 1
                umap_colors(i,:) = [0,1,0];
            else
                umap_colors(i,:) = [.5,.5,.5];
            end
        end

        subplot_column_num = 5;
        subplot(4,5,base_subplot_num + subplot_column_num)
        
        
        scatter(reduced(:,1),reduced(:,2),12,umap_colors)
        xlabel('UMAP-1')
        ylabel('UMAP-2')
        title(['cluster ',num2str(cluster_num),' | global context'])
        tt = text(.15,.9,[' # cells = ',num2str(length(find(cell_inds==1)))],'Units','normalized');
        
        
        
        sgtitle(['experiment ',num2str(experiment_num),', ',num2str(max(clustermatrix(:,experiment_num))),' clusters total'])
    
        
        
        
        
        %%% organization %%%
        if mod(cluster_num,4) == 0 || cluster_num == max(clustermatrix(:,experiment_num))
            
            %save current figure
            savename = ['Oct14clusterAnalysis_clustering',num2str(experiment_num),'_figure',num2str(fignum),'.fig'];
            savefig(savename)
            
            %make a new one
            fignum = fignum+1;
            
        end
        
        if cluster_num == max(clustermatrix(:,experiment_num))
            close all
        end
        


        
    end


end

%%

%statistics
tot = 0;
for experiment_num = 1:10

    max(clustermatrix(:,experiment_num))
    tot = tot + max(clustermatrix(:,experiment_num));
end

%%

%looking for continuous sublineages

min_length = 3; %number of generations that counts as a legitimate sublineage
threshold = .1; %percent of the cells in the cluster that have to come from the sublineage

cluster_bools = [];
c = 1;

for experiment_num = 1:10
    
    for cluster_num = 1:max(clustermatrix(:,experiment_num))
        
        
        cell_inds = find(clustermatrix(:,experiment_num)==cluster_num);
        names = all_names(cell_inds);
        
        [beginnings, endings] = detectContinuousSublineages(names,min_length);
        
        %
        for i = 1:size(beginnings,1)
            
            cells_that_meet_criteria = 0;
            
            len_of_end_cell_name = length(endings{i});
            
            num_cells_in_cluster = length(names);
            
            %count up
            for j = 1:num_cells_in_cluster
                
                if startsWith(names{j},beginnings{i}) && length(names{j}) < len_of_end_cell_name
                    
                    cells_that_meet_criteria = cells_that_meet_criteria + 1;
                    
                end
                
                
            end
            
            if cells_that_meet_criteria/num_cells_in_cluster >= threshold
                cluster_bools(c) = 1;
                c
            else
                cluster_bools(c) = 0;
                c
            end
            
        end
        %}
        
      c = c+1;  

    end
    
    
end

%%


%progenitor fates

fates_to_look_for = ["pharynx";
                     "neuron";
                     "hypoderm";
                     "bodywallmuscle";
                     "gut";
                     "OTHER"];
                 
                 
[fates] = fateOfTerminalChildren3(emb,fates_to_look_for);

c = 1;
tmp = ["seed"];
for i = 1:length(emb.embeddingStruct)
    
    if sum(fates(i,:)) == 0
        tmp(c) = string(emb.embeddingStruct(i).name);
        
        disp(emb.embeddingStruct(i).name)
        
        c = c+1;
    end

end


%%

%%% Oct 15 coarse-level sublineage pie charts %%%

%seeing how clusters map to ~20 categories


all_names = cell(1);
for i = 1:length(emb.embeddingStruct)
    
    all_names{i} = char(emb.embeddingStruct(i).name);

end
all_names = all_names';




categories_and_mainfates = ["ABalaa","neuron"; %for coloring the text labels
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
          
categories = categories_and_mainfates(:,1);
mainfates = categories_and_mainfates(:,2);
          
for experiment_num = 1:size(clustermatrix,2)
    
    figure
    sgtitle(['coarse-level sublineage composition of 10 Louvain clustering experiments, experiment ',num2str(experiment_num) newline 'label color indicates dominant progenitor type: pharynx (g), neuron (b), hypoderm (pk), bodywallmuscle (r), gut (c),  germline (pr), early/else (bk)'])
    
    %find optimal subplot dimensions
    
    num_clusters = max(clustermatrix(:,experiment_num));
    
    dim = ceil(num_clusters^.5);
    
    total_number_wedges = 0;
    total_number_distinct_progenitortypes = 0;
    
    for cluster_num = 1:num_clusters
        
        cell_inds = find(clustermatrix(:,experiment_num)==cluster_num);
        
        cluster_fates = zeros(length(cell_inds),length(categories));
        
        names = all_names(cell_inds);

        %count up what kind of cell each cell is
        for i = 1:length(cell_inds)
            
            for j = 1:length(categories)
                
                if startsWith(all_names(cell_inds(i)),categories(j))
                    
                    cluster_fates(i,j) = 1;
                    break

                end
                
            end
            
            %fill in  final column if cell doesn't fit into other categories
            if sum(cluster_fates(i,1:end-1)) == 0
                cluster_fates(i,end) = 1;
            end
            
        end
        
        %tally
        
        fate_counts = sum(cluster_fates,1)';
        vals = fate_counts(fate_counts~=0);
        labs = categories(fate_counts~=0);
        
        subplot(dim,dim,cluster_num)
        p = pie(vals,labs);
        title(['cluster ',num2str(cluster_num),' | # cells = ',num2str(length(cell_inds))])
        
        %color text labels by/count different progenitor identities in cluster
        
        progenitortypes_counts = zeros(7);%6 tissues below + other (includes early cells/
        
        for i = 2:2:length(p)%access the text objects for each wedge
            
            text_obj = p(i);
            
            wedge_fate = mainfates(find(categories == text_obj.String));

            switch wedge_fate
                case "pharynx"
                    text_obj.Color = [0,1,0];
                    progenitortypes_counts(1) = 1;
                case "neuron"
                    text_obj.Color = [0,0,1];
                    progenitortypes_counts(2) = 1;
                case "hypoderm"
                    text_obj.Color = [1,0,1];
                    progenitortypes_counts(3) = 1;
                case "bodywallmuscle"
                    text_obj.Color = [1,0,0];
                    progenitortypes_counts(4) = 1;
                case "gut"
                    text_obj.Color = [0,1,1];
                    progenitortypes_counts(5) = 1;
                case "germline"
                    text_obj.Color = [.5,0,1];
                    progenitortypes_counts(6) = 1;
                otherwise
                    %type is early/other; color stays black
                    progenitortypes_counts(1) = 7;
            end
            
        end
        
        %number wedges per cluster
        
        total_number_wedges = total_number_wedges + length(p)/2;
        
        distinct_progenitortypes = length(find(progenitortypes_counts == 1));
        total_number_distinct_progenitortypes = total_number_distinct_progenitortypes + distinct_progenitortypes;
        
    end
    
    avg_wedges_per_cluster = total_number_wedges/num_clusters;
    avg_number_distinct_progenitortypes_per_cluster = total_number_distinct_progenitortypes/num_clusters;
    disp(['Experiment ',num2str(experiment_num),': # clusters = ',num2str(num_clusters),', avg # wedges per cluster = ',num2str(avg_wedges_per_cluster),', avg # distinct progenitor types per cluster = ',num2str(avg_number_distinct_progenitortypes_per_cluster)])
    
    savefig(['Nov16_coarselevel_cluster_composition_experiment',num2str(experiment_num),'.fig'])
    
    close all
    
end




%%

%10/26

embS = emb.embeddingStruct;

%get limits for plots
lower = 0;
upper = embS(end).captureFrame;

starting_num_cells = embS(1).numCellsAtCapture;
ending_num_cells = embS(end).numCellsAtCapture;


for experiment_num = 1:size(clustermatrix,2)
    
    %for each cluster, figure out where its timepoints come from (from
    %embeddingStruct
    
    num_clusters = max(clustermatrix(:,experiment_num));
    
    subplot(2,5,experiment_num)
    
    
    data_to_plot = [];
    
    
    for cluster_num = 1:num_clusters
        
        cell_inds = find(clustermatrix(:,experiment_num)==cluster_num);
        capture_times = [embS(cell_inds).captureFrame]';

        data_to_plot = [data_to_plot;capture_times,repmat([cluster_num],length(capture_times),1)];

    end

    %put all times on a y value corresponding to the cluster
    scatter(data_to_plot(:,1),data_to_plot(:,2),15)
    ylim([0,num_clusters+1])
    ylabel('clusters')
    xlim([lower,upper])
    xlabel('capture frame cells in cluster')
    title(['experiment ',num2str(experiment_num)])
    
    
end

sgtitle(['where clusters fall in time, starting from ',num2str(starting_num_cells),' to ',num2str(ending_num_cells),' cells'])

%%

%space

%{
make movies of the embryo going through time with cells colored based on
what cluster they're in

%}

%get appropriate viewing window

embryo = emb.embdat_stabilized;


%make a master list of colors
load('fiftycolors.mat')
colorStruct = struct();
for experiment_num = 1:size(clustermatrix,2)
        
    for i = 1:size(clustermatrix,1)
        
        colorStruct(experiment_num).cols(i,:) = fiftycolors.colors(clustermatrix(i,experiment_num),:);
        
    end

end








%%


%initialize
xmin = min(embryo(1).finalpoints(:,1));
xmax = max(embryo(1).finalpoints(:,1));
ymin = min(embryo(1).finalpoints(:,2));
ymax = max(embryo(1).finalpoints(:,2));
zmin = min(embryo(1).finalpoints(:,3));
zmax = max(embryo(1).finalpoints(:,3));
%find values
for t = [emb.starttime:emb.sampling:emb.endtime]
    
    temp = min(embryo(t).finalpoints(:,1));
    if temp < xmin
        xmin = temp;
    end
    temp = min(embryo(t).finalpoints(:,2));
    if temp < ymin
        ymin = temp;
    end
    temp = min(embryo(t).finalpoints(:,3));
    if temp < zmin
        zmin = temp;
    end
    
    temp = max(embryo(t).finalpoints(:,1));
    if temp > xmax
        xmax = temp;
    end
    temp = max(embryo(t).finalpoints(:,2));
    if temp > ymax
        ymax = temp;
    end
    temp = max(embryo(t).finalpoints(:,3));
    if temp > zmax
        zmax = temp;
    end
end






%%
%render movie

clear movie

video_name = 'cluster10ONLY_video_test';

pause_length = .01;%seconds between frames

counter = 1; %to keep track of what colors to use

f = figure;
%f.Position = [100, 100, 1400, 1000];
f.Position = [616, 598, 1400, 1050];
movie_count = 1;
for t = [emb.starttime:emb.sampling:emb.endtime]
    
    %plot on blank figure
    clf
    
    points = embryo(t).finalpoints;
    num_points_in_frame = size(points,1);
    
    x = points(:,1);
    y = points(:,2);
    z = points(:,3);
    
    cells_being_viewed = [counter:counter + num_points_in_frame-1];%indices in embeddingStruct/colors
    
    
    %plot for all experiments
    
    for experiment_num = 10%1:size(clustermatrix,2)
        
        num_clusters = max(clustermatrix(:,experiment_num));
        
        colors = colorStruct(experiment_num).cols(cells_being_viewed,:);
        
        %subplot(2,5,experiment_num)
        scatter3(x,y,z,30,colors,'filled')%'MarkerFaceAlpha',.3)
        %decorate
        grid on
        xlabel('x (A->P)')
        ylabel('y (D->V)')
        zlabel('z (L->R)')
        title(['experiment ',num2str(experiment_num),' | ',num2str(num_clusters),' clusters'])
        pbaspect([1.5,1,1])
        
        %consistent viewing window
        xlim([xmin, xmax])
        ylim([ymin, ymax])
        zlim([zmin, zmax])
        
    end
    
    sgtitle(['clusters over time Oct 26 | 10 cell to 570 cell | t = ',num2str(t)])
    

    %force Matlab to plot and play slowly enough
    %pause(pause_length)
    %movie(t) = getframe(f,[110,110,1400-20,1050-20]);
    %drawnow
    %movie(t) = getframe(f,[10,10,1400-20,1000-20]);
    movie(movie_count) = getframe(f,[10,10,1400-20,1050-20]);
    
    %update counter
    counter = counter + num_points_in_frame;
    movie_count = movie_count + 1;
end

%write to video

videofile = VideoWriter(video_name,'MPEG-4');
videofile.FrameRate = 15;


open(videofile);
writeVideo(videofile,movie);
close(videofile);


%}

%%

%make color key

X = 1:50;
Y = repmat([1],1,50);

subplot(3,1,1)
scatter(X,Y,50,fiftycolors.colors,'filled','MarkerFaceAlpha',.3)
xlim([0,51])
for i = 1:50
    text(i-.2,1.2,num2str(i))
end




%%

%%% Video of development colored by clusters %%% 

%get appropriate viewing window

embryo = emb.embdat_stabilized;
experiment_num = 21;


%make a master list of colors
thirtycolors = load('thirtycontrastingcolors.mat');

video_name = 'Nov6_clustering_21_withsublineagesandfates_clusters91to120';
colors_map_to = [91:120];%numbers of the THIRTY clusters you want to map the colors to; all else will be light gray
sg_titl = ['clusters over time nov6 clustering 21 | 10 cell to twitch | clusters shown: 91:120 of 120'];

allcolors = zeros(length(emb.embeddingStruct),3);
for i = 1:size(allcolors,1)
    
    cluster_of_cell = clustermatrix(i,experiment_num);
    ind_of_col = find(colors_map_to == cluster_of_cell);
    
    if ~isempty(ind_of_col)
        allcolors(i,:) = thirtycolors.colors(ind_of_col,:);
    else
        allcolors(i,:) = [.5,.5,.5];
    end
    
end
celltypes = ["pharynx","neuron","bodywallmuscle","hypoderm","gut"];
lineages = ["ABa","ABp","C","D","E","MS","P"];
coarselevelcolors = colorByCoarseLevelSublineage_special(emb.embeddingStruct);
[celltypecolors,~,~] = colorArrays_Nov9(emb.embeddingStruct,celltypes,lineages);

%{
colorStruct = struct();
for experiment_num = 1:size(clustermatrix,2)
        
    for i = 1:size(clustermatrix,1)
        
        colorStruct(experiment_num).cols(i,:) = fiftycolors.colors(clustermatrix(i,experiment_num),:);
        
    end

end
%}


%initialize
xmin = min(embryo(1).finalpoints(:,1));
xmax = max(embryo(1).finalpoints(:,1));
ymin = min(embryo(1).finalpoints(:,2));
ymax = max(embryo(1).finalpoints(:,2));
zmin = min(embryo(1).finalpoints(:,3));
zmax = max(embryo(1).finalpoints(:,3));
%find values
for t = [emb.starttime:emb.sampling:emb.endtime]
    
    temp = min(embryo(t).finalpoints(:,1));
    if temp < xmin
        xmin = temp;
    end
    temp = min(embryo(t).finalpoints(:,2));
    if temp < ymin
        ymin = temp;
    end
    temp = min(embryo(t).finalpoints(:,3));
    if temp < zmin
        zmin = temp;
    end
    
    temp = max(embryo(t).finalpoints(:,1));
    if temp > xmax
        xmax = temp;
    end
    temp = max(embryo(t).finalpoints(:,2));
    if temp > ymax
        ymax = temp;
    end
    temp = max(embryo(t).finalpoints(:,3));
    if temp > zmax
        zmax = temp;
    end
end


%render movie

clear movie


pause_length = .01;%seconds between frames

counter = 1; %to keep track of what colors to use

f = figure;
%f.Position = [100, 100, 1400, 1000];
f.Position = [616, 598, 1400, 1050];
movie_count = 1;
for t = [emb.starttime:emb.sampling:emb.endtime]
    
    %plot on blank figure
    clf
    
    points = embryo(t).finalpoints;
    num_points_in_frame = size(points,1);
    
    x = points(:,1);
    y = points(:,2);
    z = points(:,3);
    
    cells_being_viewed = [counter:counter + num_points_in_frame-1];%indices in embeddingStruct/colors
    
    
    %plot for all experiments
    
    for experiment_num = 10%1:size(clustermatrix,2)
        
        num_clusters = max(clustermatrix(:,experiment_num));
        
        colors = allcolors(cells_being_viewed,:);
        
        coarsecolors = coarselevelcolors(cells_being_viewed,:);
        
        fatecolors = celltypecolors(cells_being_viewed,:);
        
        subplot(1,3,2)
        view(3)
        hold on
        %subplot(2,5,experiment_num)
        scatter3(x,y,z,30,colors,'filled')%'MarkerFaceAlpha',.3)
        scatter3(x,y,z,100,colors,'filled','MarkerFaceAlpha',.3)%'MarkerFaceAlpha',.3)
        %decorate
        grid on
        xlabel('x (A->P)')
        ylabel('y (D->V)')
        zlabel('z (L->R)')
        title('clustering')
        pbaspect([1.5,1,1])
        
        %consistent viewing window
        xlim([xmin, xmax])
        ylim([ymin, ymax])
        zlim([zmin, zmax])
        
        subplot(1,3,1)
        view(3)
        hold on
        %subplot(2,5,experiment_num)
        scatter3(x,y,z,30,coarsecolors,'filled')%'MarkerFaceAlpha',.3)
        scatter3(x,y,z,100,coarsecolors,'filled','MarkerFaceAlpha',.3)%'MarkerFaceAlpha',.3)
        %decorate
        grid on
        xlabel('x (A->P)')
        ylabel('y (D->V)')
        zlabel('z (L->R)')
        title('coarse-level sublineages')
        pbaspect([1.5,1,1])
        
        %consistent viewing window
        xlim([xmin, xmax])
        ylim([ymin, ymax])
        zlim([zmin, zmax])
        
        subplot(1,3,3)
        view(3)
        hold on
        %subplot(2,5,experiment_num)
        scatter3(x,y,z,30,fatecolors,'filled')%'MarkerFaceAlpha',.3)
        scatter3(x,y,z,100,fatecolors,'filled','MarkerFaceAlpha',.3)%'MarkerFaceAlpha',.3)
        %decorate
        grid on
        xlabel('x (A->P)')
        ylabel('y (D->V)')
        zlabel('z (L->R)')
        title('fates of terminal cells')
        pbaspect([1.5,1,1])

        %consistent viewing window
        xlim([xmin, xmax])
        ylim([ymin, ymax])
        zlim([zmin, zmax])
        
        
    end
    
    sgtitle([sg_titl newline 'frame = ',num2str(t)])
    
    movie(movie_count) = getframe(f,[10,10,1400-20,1050-20]);
    
    %update counter
    counter = counter + num_points_in_frame;
    movie_count = movie_count + 1;
end

%write to video

videofile = VideoWriter(video_name,'MPEG-4');
videofile.FrameRate = 15;


open(videofile);
writeVideo(videofile,movie);
close(videofile);



%%

%make color key

X = 1:30;
Y = repmat([1],1,30);
figure
scatter(X,Y,50,thirtycolors.colors,'filled')
xlim([0,31])
for i = 1:30
    text(i-.2,1.2,num2str(i))
end
%%
%make color key for sublineages


categories_and_mainfates = ["ABalaa","neuron"; %for coloring the text labels
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

X = 1:29;
Y = repmat([1],1,29);

cols = load('thirtycontrastingcolors.mat');
cols.colors(end-1,:) = [];
cols.colors(end,:) = [.5,.5,.5];

figure
scatter(X,Y,200,cols.colors,'filled')
xlim([0,31])
for i = 1:29
    text(i-.2,1.2,categories_and_mainfates(i))
end

title('sublineage color key')


%%

%%% look for apical/basal clusters %%%

%plot single timepoints

experiment_num = 10; %which clstering experiment do you want to look at


time_to_view = embeddingStruct(end).captureFrame;

%figure out what indices to plot
inds_to_plot = find([embeddingStruct.captureFrame] == time_to_view);

%get points
pts = [];
for i = 1:length(inds_to_plot)
    
    pts = [pts;embeddingStruct(inds_to_plot(i)).realSpaceCoord];
        
end
cls = colorStruct(experiment_num).cols(inds_to_plot,:);

figure
view(3)
hold on 
scatter3(pts(:,1),pts(:,2),pts(:,3),15,cls,'filled')
scatter3(pts(:,1),pts(:,2),pts(:,3),150,cls,'filled','MarkerFaceAlpha',.3)
xlabel('x (A->P)')
ylabel('y (D->V)')
zlabel('z (L->R)')
title('dispim emb 2 last timepoint, colored by clustering 10')


%%

%compsitopn of cluseters in last time
%   uses archietcture above

%plot single timepoints

experiment_num = 10; %which clstering experiment do you want to look at

cluster_to_view = 3;%cluster to color in second subplot


time_to_view = embeddingStruct(end).captureFrame;

%figure out what indices to plot
inds_to_plot = find([embeddingStruct.captureFrame] == time_to_view);

%get points
pts = [];
for i = 1:length(inds_to_plot)
    
    pts = [pts;embeddingStruct(inds_to_plot(i)).realSpaceCoord];
        
end
cls = colorStruct(experiment_num).cols(inds_to_plot,:);

specific_colors = repmat([.5,.5,.5],length(embeddingStruct),1);
specific_colors(clustermatrix(:,experiment_num)==cluster_to_view,:) = repmat([1,0,0],length(find(clustermatrix(:,experiment_num)==cluster_to_view)==1),1);
specific_colors = specific_colors(inds_to_plot,:);

celltypes = ["pharynx","neuron","bodywallmuscle","hypoderm","gut"];
lineages = ["ABa","ABp","C","D","E","MS","P"];

[celltypecolors,~,~] = colorArrays_Nov9(embeddingStruct,celltypes,lineages);

celltypecolors = celltypecolors(inds_to_plot,:);


subplot(1,3,1)
view(3)
hold on 

scatter3(pts(:,1),pts(:,2),pts(:,3),15,cls,'filled')
scatter3(pts(:,1),pts(:,2),pts(:,3),150,cls,'filled','MarkerFaceAlpha',.3)
xlabel('x (A->P)')
ylabel('y (D->V)')
zlabel('z (L->R)')
title('colored by clustering 10')

subplot(1,3,2)
view(3)
hold on 
scatter3(pts(:,1),pts(:,2),pts(:,3),15,specific_colors,'filled')
scatter3(pts(:,1),pts(:,2),pts(:,3),150,specific_colors,'filled','MarkerFaceAlpha',.3)
xlabel('x (A->P)')
ylabel('y (D->V)')
zlabel('z (L->R)')
title(['cluster',num2str(cluster_to_view)])


subplot(1,3,3)
view(3)
hold on 
scatter3(pts(:,1),pts(:,2),pts(:,3),15,celltypecolors,'filled')
scatter3(pts(:,1),pts(:,2),pts(:,3),150,celltypecolors,'filled','MarkerFaceAlpha',.3)
xlabel('x (A->P)')
ylabel('y (D->V)')
zlabel('z (L->R)')
title('colored by terminal fate')


sgtitle('dispim emb 2 (through comma) last timepoint (570 cells)')





%%

%10/26

% how many clusters does a given cell (ID) fall into?

number_of_experiments = size(clustermatrix,2);

ordered_names_unique = unique(all_names);
number_of_unique_names = length(ordered_names_unique);

%make string form of all_names
all_names_string = ["seed"];
for i = 1:length(all_names)
    all_names_string(i) = string(all_names{i});
end

%matrix to store numbers 
%   rows are *unique* cell names 
%   columns are experiments 1:10
%   entries are numbers of distinct clusters that cell ID falls into

vals_matrix = zeros(number_of_unique_names,number_of_experiments);
median_time_matrix = zeros(number_of_unique_names,1);
percent_of_timepoints_that_are_transitions = zeros(number_of_unique_names,1);
when_the_transitions_occur = [];

for unique_cell_id = 1:number_of_unique_names
    
    for experiment_num = 10%1:number_of_experiments
        
        cell_inds_to_look_at = all_names_string == ordered_names_unique{unique_cell_id};
        
        %log what the median time of the existence of that cell was
        median_time_matrix(unique_cell_id) = median([emb.embeddingStruct(cell_inds_to_look_at).captureFrame]);
        
        clusters_the_cell_is_in = clustermatrix(cell_inds_to_look_at,experiment_num);
        
        %percent of timepts tht are transitions
        num_transitions = 0;
        for i = 2:length(clusters_the_cell_is_in)
            if clusters_the_cell_is_in(i) ~= clusters_the_cell_is_in(i-1)
                when_the_transitions_occur = [when_the_transitions_occur;i/length(clusters_the_cell_is_in)];
                num_transitions = num_transitions + 1;
            end
        end
        percent_transitions = num_transitions/length(clusters_the_cell_is_in);
        percent_of_timepoints_that_are_transitions(unique_cell_id) = percent_transitions;
        
        
        %count up # clusters a given cell id can be found in for experiment
        number_of_distinct_clusters = length(unique(clusters_the_cell_is_in));
        
        %record
        vals_matrix(unique_cell_id,experiment_num) = number_of_distinct_clusters;
        
        if number_of_distinct_clusters > 1
            
            %these are messing things up and dont display the right
            %information
            %all_numbers_in_line = str2double(sprintf('%d',clusters_the_cell_is_in')); %from internet 
            %sprintf('%f',all_numbers_in_line) %also from internet

            disp(unique_cell_id)
            
            disp(clusters_the_cell_is_in)

            
        end
        
        
        
    end
    
end

%%
%figuring out how the assortment into clusters over time works
%only on exp 10. (experiment_num = 10 in for loop in prev section)

num_clusters_each_cell = vals_matrix(:,experiment_num);

%mean number of clusters separated into
mean_num_clusters_per_cellID = mean(num_clusters_each_cell);

%histogram
histogram(num_clusters_each_cell);

%number of clusters each cell falls into (ordered by lineage, time)
    %already done, in google doc
    
%number of clusters each cell falls into (ordered by time only)
scatter(median_time_matrix,num_clusters_each_cell);


%what percent of a cell's timepoints are a transition (cluster assignment is different from the previous timepoint

percent_transitions_transitioningcellsonly = percent_of_timepoints_that_are_transitions(percent_of_timepoints_that_are_transitions~=0);

histogram(percent_transitions_transitioningcellsonly,20)


%plot the above

subplot(2,2,1)
histogram(num_clusters_each_cell);
title('number of cluters each cell is assigned to')

subplot(2,2,2)
scatter(median_time_matrix,num_clusters_each_cell);
title(['number of clusters each cell is assigned to (ordered chronologically)' newline '(cells placed at the median capture frame at which they exist)'])
xlabel('capture frame')

subplot(2,2,3)
histogram(percent_transitions_transitioningcellsonly,20)
title(['percent of a cells timepoints assigned to a different cluster than the timepoint before them' newline '(shown only for cells assigned to more than 1 cluster)'])

subplot(2,2,4)
histogram(when_the_transitions_occur,20)
title(['when in a cells lifetime its timepoints are assigned to a different cluster than the timepoint before them' newline '(shown only for cells assigned to more than 1 cluster)'])
xlabel("relative position through a cell's lifetime scaled [0,1]")

sgtitle('distribution of time-cells into clusters Oct 28')




%%
%define the specific set of cells you label
%   discard all past: AB16, MS4, D2, E

labels = ordered_names_unique;
for i = 1:number_of_unique_names
    
    nm = ordered_names_unique{i};
    l = length(nm);
    
    if startsWith(nm,'AB') && l ~= 6
        labels{i} = '';
    elseif startsWith(nm,'MS') && l ~= 4
        labels{i} = '';
    elseif startsWith(nm,'C') && l ~= 3
        labels{i} = '';
    elseif startsWith(nm,'D') && l ~= 2
        labels{i} = '';
    elseif startsWith(nm,'E') && l ~= 1
        labels{i} = '';
    elseif startsWith(nm,'P') && nm ~= string('P2') %eliminate close similar names
        labels{i} = '';
    elseif startsWith(nm,'Z') && nm ~= string('Z2') %eliminate close similar names
        labels{i} = '';
    end
    
end

%make ticks

ticks = [];
for i = 1:length(labels)
    if string(labels{i})~=''
        ticks = [ticks;i];
    end
end
labels = labels(ticks); %should discard all the empty ones


%make scatters 

X = 1:number_of_unique_names;

for i = 1:number_of_experiments
    
    subplot(10,1,i)
    Y = vals_matrix(:,i);
    scatter(X,Y,15)
    xticks(ticks)
    xticklabels(labels)
    ylim([0,max(Y)+1])
    xlim([0,number_of_unique_names])
    
    
    num_clusters = max((clustermatrix(:,i)));
    title(['experiment ',num2str(i),' | contains ',num2str(num_clusters),' clusters'])
    
    
end

meanss = "1.2974, 1.4493, 1.5147, 1.5302, 1.5752, 1.6307, 1.7083, 1.7149, 1.6961, 1.7753";
sgtitle(['number of clusters each cell ID gets placed into (in early October Louvian clustering) -- Oct 26' newline 'experiment means: ',meanss])

%probability that all cell instances are in same cluster
disp('probabilities that all cell instances are in same cluster, by cluster')
for i = 1:10
    
    %num cell IDs in one cluster
    num_cell_IDs_in_one_cluster = length(vals_matrix(vals_matrix(:,i)==1,1));
    
    prob = num_cell_IDs_in_one_cluster/number_of_unique_names;
    
    disp(prob)
    
    
    
    
end

%%

%Oct 26

%number of time-cells and number of distinct cell IDs per cluster

for i = 1:number_of_experiments
    
    num_clusters = max(clustermatrix(:,i));
    
    vals_holder = zeros(num_clusters,2); %first column is number of time-cells, second is number of distinct cell IDs
    
    for j = 1:num_clusters
        
        num_cells_in_cluster = length(clustermatrix(clustermatrix(:,i)==j,i));
        num_unique_cell_IDs_in_cluster = length(unique(all_names(clustermatrix(:,i)==j)));
        
        vals_holder(j,:) = [num_cells_in_cluster,num_unique_cell_IDs_in_cluster];
        
    end
    dexes = [1:5,11:15];
    subplot_index = dexes(i);
    
    subplot(4,5,subplot_index)
    
    bar(vals_holder(:,1))
    %set(gca, 'YScale', 'log')
    title(['experiment ',num2str(i)])
    
    subplot(4,5,subplot_index+5)
    
    bar(vals_holder(:,2),'Facecolor',[0.8500, 0.3250, 0.0980])
    %set(gca, 'YScale', 'log')
    %title(['experiment ',num2str(i)])

    
end

sgtitle('number of time-cells (b) and unique cell IDs (o) in each cluster (early October Louvain clustering) Oct 26')



%%

num_clusters = max(clustermatrix(:,experiment_num));
for experiment_num = 1:size(clustermatrix,2)
    
    
    num_clusters = max(clustermatrix(:,experiment_num));
    
    for cluster_num = 1:num_clusters
        
        
        clus_names=all_names(clustermatrix(:,experiment_num)==cluster_num);
        
        if ~isempty(find(clus_names == "ABal"))
            
            experiment_num
            cluster_num
            
            
            break
            
            
        end
    
    end
end

%%

for i = 1:size(fates,1)
    
    if sum(fates(i,:)) == 0
        disp(i)
    end
    
end
%%

%pie text properties

dataa = [1 1 1 1]';
labelll = ["a";"b";"c";"d"];
figure
pp = pie(dataa,labelll);




%%

%trying to read in parts list

C = readcell('partslist_web.txt');

%%

%find cells that are not classifed

fates_to_look_for = ["pharynx";
                     "neuron";
                     "hypoderm";
                     "bodywallmuscle";
                     "gut";
                     "OTHER"];

[terminal_fates] = fateOfTerminalChildren3(emb,fates_to_look_for);

%find the cells that are left uncategorized
unclassified_inds = find(sum(terminal_fates,2)==0);

unclassified_cells = all_names(unclassified_inds);

unclassified_cells_unique = unique(unclassified_cells);



%{
%%
%normalize to percent composition of each term fate
terminal_composition = terminal_fates;
for i = 1:size(terminal_composition,1)
    if sum(terminal_composition(i,:)) ~= 0 %if all zeros (programmed cell death), leave as all zeros
        terminal_composition(i,:) = terminal_composition(i,:)/sum(terminal_composition(i,:));
    end
end
%}


%%

%{

there are cases to this: there are
    - unclassified cells that are not present in the partslist (ex.
    ABalaaaalar)
        - in some of these, the sibling  is present (ex. ABalaaaarla is not
        there, but ABalaaaarlap is)
    - unclassified cells that are a type of cells not considered in the
    other lists (ex. ABalaaaalp is inner labial sheath)
    

%}

%looking to see how easy it is to find an unclassified cell in the master
%partslist

unclassified_cells_unique_string = ["seed"];
for i = 1:length(unclassified_cells_unique)
    unclassified_cells_unique_string(i) = string(unclassified_cells_unique{i});
end

second_column_partslist = ["seed"];
for i = 1:size(partslist,1)
    second_column_partslist(i) = string(partslist{i,2});
end

still_unclassifiable = [];
become_classifiable = [];
for i = 1:1:length(unclassified_cells_unique)
    
    if ~isempty(find(second_column_partslist == unclassified_cells_unique{i}))
        become_classifiable = [become_classifiable;i];
    else
        still_unclassifiable = [still_unclassifiable;i];
    end
    
end

%of the ones that are still unclassifiable


%question: are they the *children* of recognized tissues types?

alltermcells_string = ["seed"];%just from the tissue lists

for i = 1:length(alltermcells)
    alltermcells_string(i) = string(alltermcells{i});
end

for i = 1:length(still_unclassifiable)
    
    cell_name = unclassified_cells_unique(still_unclassifiable(i));
    

    if sum(startsWith(alltermcells,'cell_name')) > 0
        disp(cell_name)
    end
    
end
%answer-->no.

%so. there are cells in the embryo not in the tissue lists and not in the
%master partslist.(still_unclassifiable). number exceeds number of
%programmed cell deaths in c elegans (131?).




%%

function [fates] = fateOfTerminalChildren2(emb,fates_to_look_for)




embdat = emb.embdat_stabilized;
embeddingStruct = emb.embeddingStruct;


embStruct_times = [emb.starttime:emb.sampling:emb.endtime];
last_time = embStruct_times(end);


%get fates of terminal cells (last frame)

last_frame_l = size(embdat(last_time).cellnames,1);
big_l = length(embeddingStruct);

last_time_fates = ["seed"];
c = 1;
for i = big_l - last_frame_l + 1:big_l
    
    last_time_fates(c) = embeddingStruct(i).celltype;
    
    c = c+1;
    
end

fates = zeros(length(embeddingStruct),6);
c = 1;


for t = emb.starttime:emb.sampling:emb.endtime %may throw an error for a differnt emb--no time to deal with now. has to do with how last entry in embdat (which points to nothing) is handled
    for i = 1:size(embdat(t).cellnames,1)
        
        %t
        %i
        [terminal_inds] = terminalChildren(t,i,emb);
        
        term_fates = last_time_fates(terminal_inds);
        
        for j = 1:6
            
            fates(c,j) = length(find(term_fates==fates_to_look_for(j)));
            
        end
        
        
        %{
        debugging 
        if sum(fates(c,:)) == 0
            t
            i
            embdat(t).cellnames(i)
            term_fates

        end
        %}
        
        
        c = c+1;
        
        
        

    end
end


end

function [cells_to_look_at] = terminalChildren(start_time,start_ind,emb)

%replace forward traverse, find the indices in the embryo struct of the
%terminal children of the idicated cell


cells_to_look_at = start_ind;
time_pt = start_time;

%find point to end search
times_in_embeddingStruct = [emb.starttime:emb.sampling:emb.endtime];
cutoff = times_in_embeddingStruct(end);


while true
    
    
    if time_pt == cutoff || isempty(cells_to_look_at)
        break
    end
    %note: if you input the cutoff time, it'll just return the start_ind
    
    
    cells_to_look_at_next = [];
    
    for i = 1:length(cells_to_look_at)
        
        %time_pt
        %i
        cells_to_look_at_next = [cells_to_look_at_next,emb.embdat_stabilized(time_pt).suc(cells_to_look_at(i),:)];
    end
    
    %remove all -1s (dead/nonexistent cells
    
    cells_to_look_at_next = cells_to_look_at_next(cells_to_look_at_next~=-1);
    
    %update
    cells_to_look_at = cells_to_look_at_next;
    time_pt = time_pt+1;

end

time_pt


end

function [fates] = fateOfTerminalChildren3(emb,fates_to_look_for)

embdat = emb.embdat_stabilized;
embeddingStruct = emb.embeddingStruct;

%get fates of terminal cells
load('allpharynx.mat','allpharynx');
load('allneuron.mat','allneuron');
load('allbodywallmuscle.mat','allbodywallmuscle');
load('allhypoderm.mat','allhypoderm');
load('allgut.mat','allgut');


alltermcells = [allpharynx;
                allneuron;
                allhypoderm;
                allbodywallmuscle;
                allgut];
    
alltermfates = [repmat({'pharynx'},size(allpharynx,1),1);
                repmat({'neuron'},size(allneuron,1),1);
                repmat({'hypoderm'},size(allhypoderm,1),1);
                repmat({'bodywallmuscle'},size(allbodywallmuscle,1),1);
                repmat({'gut'},size(allgut,1),1)];
            
%basic case

fates = zeros(length(embeddingStruct),6);
for i = 1:length(emb.embeddingStruct)
    
    cell_name = emb.embeddingStruct(i).name;
    
    term_fates = alltermfates(startsWith(alltermcells,cell_name));
    
    for j = 1:6    
        fates(i,j) = length(find(term_fates==fates_to_look_for(j)));    
    end
    
end


%special cases


%//

end

function [celltypecolors,lineagecolors,timecolors] = colorArrays_Nov9(embeddingStruct,celltypes,lineages)

%note: timecolors edited to exclude the low hue value red from the HSV
%spectrum. Scales to [.3,1], no [0,1], as do other versions of this
%function

numPointsInEmbedding = size(embeddingStruct,2);
endtime = embeddingStruct(end).captureFrame;


%color scheme

%temporal evolution of embryo
timecolors=zeros(numPointsInEmbedding,3);
min_frame = embeddingStruct(1).captureFrame;
max_frame = embeddingStruct(end).captureFrame;
for i = 1:numPointsInEmbedding
    
    sandwich = [min_frame,embeddingStruct(i).captureFrame,max_frame];
    sandwich = rescale(sandwich,.3,1);%.3 = green, 1 = red
    hsvcolor = [sandwich(2),1,1];    
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

function [coarselevelcolors] = colorByCoarseLevelSublineage_special(embeddingStruct)

%color by these sublineages

cols = load('thirtycontrastingcolors.mat');
cols.colors(end,:) = [.5,.5,.5];

categories_and_mainfates = ["ABalaa","neuron"; %for coloring the text labels
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
          
sublineages = categories_and_mainfates(:,1);
          
coarselevelcolors = zeros(length(embeddingStruct),3);

for i = 1:length(embeddingStruct)
    
    
    for j = 1:length(sublineages)-1%exclude the catch-all early/other category
        
        if startsWith(embeddingStruct(i).name,sublineages(j))
            
            coarselevelcolors(i,:) = cols.colors(j,:);
            
        end

    end
    
    %if still not categorized, give it the 30th color
    
    if coarselevelcolors(i,:) == [0,0,0]
        
        coarselevelcolors(i,:) = cols.colors(end,:);
        
    end

end

end



