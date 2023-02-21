%dump pile for pieces of code discarded from the original pipeline


%{
%remove rotation during development
[embdat_stabilized] = internallyAlignNamedEmbryo(embdat,tstart,tend);
[cells_stabilized] = parseCellsFromEmb(embdat_stabilized,tend);
%}






%coloring points in embedding by terminal fate

%old implementation
%{
%cell type
if coloring == "celltype"
    %load cells from file
    %**should switch it back to 'allneurons.mat?'
    allpharynx=load('allpharynx.mat');
    allpharynx=allpharynx.allpharynx;
    allneuron=load('allneuron.mat');
    allneuron=allneuron.allneuron;
    allbodywallmuscle=load('allbodywallmuscle.mat');
    allbodywallmuscle=allbodywallmuscle.allbodywallmuscle;
    allhypoderm=load('allhypoderm.mat');
    allhypoderm=allhypoderm.allhypoderm;
    allgut=load('allgut.mat');
    allgut=allgut.allgut;

    %for objective function
    celltypes = zeros(size(percellfeature,1),1);
    %hard coding bc this is useless anyway and dont want to learn key val pairs
    %right now:
    %pharynx 1
    %neuron 2
    %bodywallmuscle 3
    %hypoderm 4
    %gut 5
    organ_names = ["pharynx: 1","neuron: 2","body wall muscle: 3","hypoderm: 4","gut: 5"];

    %color scheme
    c=1;
    colors=zeros(size(percellfeature,1),3);
    for t=starttime:sampling:endtime %length(embdat)
        for j=1:length(embdat(t).finalpoints)
            name=embdat(t).names{j};
            ispharynx=false;
            isneuron=false;
            isbodywallmuscle=false;
            ishypoderm=false;
            isgut=false;

            %categorize cell, log in celltypes{}
            if ~isempty(find(strcmp(name,allpharynx), 1))
                ispharynx=true;
                celltypes(c) = 1;
            end
            if ~isempty(find(strcmp(name,allneuron), 1))
                isneuron=true;
                celltypes(c) = 2;
            end
            if ~isempty(find(strcmp(name,allbodywallmuscle), 1))
                isbodywallmuscle=true;
                celltypes(c) = 3;
            end
            if ~isempty(find(strcmp(name,allhypoderm), 1))
                ishypoderm=true;
                celltypes(c) = 4;
            end
            if ~isempty(find(strcmp(name,allgut), 1))
                isgut=true;
                celltypes(c) = 5;
            end

            %add color
            if (ispharynx)
                colors(c,:)=[0,1,0]; %pharynx is green
            elseif (isneuron)
                colors(c,:)=[0,0,1]; %neuron is blue
            elseif (isbodywallmuscle)
                colors(c,:)=[1,0,0]; %bodywallmuscle is red
            elseif (ishypoderm)
                colors(c,:)=[1,0,1]; %hypoderm is pink (fuschia in RGB speak)
            elseif (isgut)
                colors(c,:)=[0,1,1]; %gut is cyan (aqua in RGB speak)
            else
                colors(c,:)=[.25,.25,.25]; %else is gray 
            end

            c=c+1;
        end
    end
    disp('cell type annotation worked')
end
%}