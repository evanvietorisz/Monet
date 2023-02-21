function [embryo_stabilized]=internallyAlignNamedEmbryo(embryo,tstart,tend)
%removes global rotations of the embryo during development. since alignment
%is not an issue while embryo is small, uses frame at which embryo contains
%four cells to start alignment (necessary because before then cell division
%does not follow the naming rule)

%embryo is embryo struct
%tstart is first frame in the dataset
%tend is last (relevant) frame in dataset--should crop out times after
%which the embryo starts moving around

%%%

%find first frame containing 4 cells
startframe = tstart;
num_cells = size(embryo(1).finalpoints,1);
while(num_cells < 4)
    startframe = startframe + 1;
    num_cells = size(embryo(startframe).finalpoints,1); 
end

%copy over struct structure to stabilized version
embryo_stabilized=embryo;

%stabilize embryo affinely over time period provided
for t=startframe:tend-1

    names_1=embryo(t).names;
    names_2=embryo(t+1).names;

    %aligning t+1 to t so t needs to be aready stabilized
    pos_1=embryo_stabilized(t).finalpoints;

    pos_2=embryo(t+1).finalpoints;

    %c=1;
    matchpoint1=[];
    matchpoint2=[];
    
    for i=1:length(names_1)
        
        %generate arrays of positions of correspoinding points
        for j=1:length(names_2)
            if strcmp(names_1{i},names_2{j}) && (isempty(strfind(names_1{i},'Nuc')))
                %indmatch1=j;
                matchpoint1=[matchpoint1;pos_1(i,:)]; % 3 x number cells i at time time
                matchpoint2=[matchpoint2;pos_2(j,:)]; % 3 x number cells j at time time
            end
        end
    end
    
    %rename
    lmpositions1=matchpoint1;
    lmpositions2=matchpoint2;
    
    %compute transformation mapping points at t+1 to t
        % / operator finds least squares solution X to X*A = B
        %add in extra spatial dimension for affine transformation
    transform2to1=[lmpositions1,ones(size(lmpositions1,1),1)]'/[lmpositions2,ones(size(lmpositions2,1),1)]'; % 4 x 4
    
%{
    if t == 9 || t == 9 || t == 10
        disp(t)
        disp([lmpositions1,ones(size(lmpositions1,1),1)]')
        disp([lmpositions2,ones(size(lmpositions2,1),1)]')
        disp(transform2to1)
    end
%}
    
    %apply transformation on points at t+1
    embryo_stabilized(t+1).finalpoints=(transform2to1*[embryo(t+1).finalpoints,ones(size(embryo(t+1).finalpoints,1),1)]')';
    %remove extra dimension used for transformation
    embryo_stabilized(t+1).finalpoints=embryo_stabilized(t+1).finalpoints(:,1:3);
end