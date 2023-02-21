%visualize normals
%   time is time at which to view embryo.
%   knorm is number of nearest neighbors used to calculate normals.
%   c is initial grouping factor representing the number of connections
%   in the graph that determines the propagation trajectory over the cell.
%   note: c <= knorm

%parameters
time = 52; 
knorm = 5;
c = 2;

%load and prepare data
tstart=1;
tend=305;
anisotropy=1;
emb='/Users/evanvietorisz/Documents/MATLAB/Decon_emb1_MGedits.zip';
templocation='temp_unzip/';

%unzip zipfile to temp file
if ~exist(emb,'file')
    errors.zipfilemissing=1;
    return
end
try
    unzip(emb,templocation);
catch exception
    errors.zipfilecorrupted=1;
    return
end

[cells,embdat] = loadcells_unnamed(templocation,tend,4,false );
rmdir(templocation,'s');

[embdat_stabilized] = internallyAlignNamedEmbryo(embdat,tstart,tend,anisotropy);
[cells_stabilized] = parseCellsFromEmb( embdat_stabilized,tend );

cells_backup=cells;
cells=cells_stabilized;

points = embdat_stabilized(time).finalpoints;

%calc normals
[untreated_norms,treated_norms] = cohereNormals(points,knorm,c);

%plot
figure

subplot(1,2,1)
quiver3(points(:,1),points(:,2),points(:,3),untreated_norms(:,1),untreated_norms(:,2),untreated_norms(:,3),.4)
xlabel('x')
ylabel('y')
zlabel('z')
title('Untreated Normals')
pbaspect([1 1 1])

subplot(1,2,2)
quiver3(points(:,1),points(:,2),points(:,3),treated_norms(:,1),treated_norms(:,2),treated_norms(:,3),.4)
xlabel('x')
ylabel('y')
zlabel('z')
title('Treated Normals')
pbaspect([1 1 1])

titlestr = strcat('Time=',num2str(time),', knorm=',num2str(knorm),', c=',num2str(c));
sgtitle(titlestr)




