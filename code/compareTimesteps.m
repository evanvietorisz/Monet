function [] = compareTimesteps(emb,times)
%conveniently compare different timesteps in embryo's development
%emb is embryo struct
%times is array of times to compare separately

%get appropriate viewing window

%initialize
xmin = min(emb(times(1)).finalpoints(:,1));
xmax = max(emb(times(1)).finalpoints(:,1));
ymin = min(emb(times(1)).finalpoints(:,2));
ymax = max(emb(times(1)).finalpoints(:,2));
zmin = min(emb(times(1)).finalpoints(:,3));
zmax = max(emb(times(1)).finalpoints(:,3));

%find values
for t = times(1):times(end)
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

%plot for desired times
for i = 1:length(times)
    t = times(i);
    
    figure
    scatter3(emb(t).finalpoints(:,1),emb(t).finalpoints(:,2),emb(t).finalpoints(:,3))
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title(['t = ',num2str(t)])
    xlim([xmin, xmax])
    ylim([ymin, ymax])
    zlim([zmin, zmax])
end

end

