function [] = showEmbMovie(emb,tstart,tend,pause_length)
%plots the point cloud over the development of the embryo for checking
%stabilization

%emb is stabilized or unstabilized embryo
%tstart,tend are first and last frame plotted
%pause_length is the duration between frames in seconds

%get appropriate viewing window

%initialize
xmin = min(emb(1).finalpoints(:,1));
xmax = max(emb(1).finalpoints(:,1));
ymin = min(emb(1).finalpoints(:,2));
ymax = max(emb(1).finalpoints(:,2));
zmin = min(emb(1).finalpoints(:,3));
zmax = max(emb(1).finalpoints(:,3));

%find values
for t = tstart:tend
    temp = min(emb(t).finalpoints(:,1));
    if temp < xmin
        xmin = temp;
    end
    temp = min(emb(t).finalpoints(:,2));
    if temp < ymin
        ymin = temp;
    end
    temp = min(emb(t).finalpoints(:,3));
    if temp < zmin
        zmin = temp;
    end
    
    temp = max(emb(t).finalpoints(:,1));
    if temp > xmax
        xmax = temp;
    end
    temp = max(emb(t).finalpoints(:,2));
    if temp > ymax
        ymax = temp;
    end
    temp = max(emb(t).finalpoints(:,3));
    if temp > zmax
        zmax = temp;
    end
end

%render movie
for t = tstart:tend
    %plot on blank figure
    clf
    
    points = emb(t).finalpoints;
    x = points(:,1);
    y = points(:,2);
    z = points(:,3);
    
    scatter3(x,y,z)
    
    %decorate
    grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title(['t = ',num2str(t)])
    
    %consistent viewing window
    xlim([xmin, xmax])
    ylim([ymin, ymax])
    zlim([zmin, zmax])
    
    %force Matlab to plot and play slowly enough
    pause(pause_length)
end
end

