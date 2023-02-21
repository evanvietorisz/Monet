function [transformedemb]=coalignNamedEmbryosPerTime2(embref,tstartref,tendref,xyscaleref,emb,tstart,tend,anisotropy,xyscale,centering_time)

%given a template embryo and list of n embryos to be aligned to it
%along with start and end time for computing linear temporal alignment 
%for end frame compute transform for each emb 
%do this anew at each frame

%align an input embryo to a reference embryo at each timestep
%embref is reference embryo struct
%tstartref is start time for reference
%tendref is end time for reference
%anisotropy_ref is anisotropy of reference
%xyscaleref is scale of x-y coordinates in microns for reference
%emb is input embryo struct to be aligned
%tstart is start time for input
%tend is end time for input
%anisotropy is anisotropy of input
%xyscale is scale of x-y coordinates in microns of input

%end times for both reference and input must be determined manually. are
%times when the embryo starts moving, at which point the point clouds are
%no longer meaningful.




%landmarks
%seam cell landmarks
landmarks={'ABplaaappa';'ABplaaappp';'ABarppaaap';'ABarppapaa';'ABarppapap';'ABplappapa';'ABarppappa';'ABplapapaa';'ABarppappp';'ABarpapppp';'ABarpppaap';'ABarppppaa';'ABarppppap';'ABprappapa';'ABarpppppa';'ABprapapaa';'ABarpppppp'};
%gut
landmarks2={'Ealaad';'Earaad';'Ealaav';'Earaav';'Earpa';'Ealpa';'Earap';'Ealap';'Epraa';'Eplaa';'Earpp';'Ealpp';'Eprap';'Eplap';'Eprpa';'Eplpa';'Eprppa';'Eplppa';'Eplppp';'Eprppp'};
%headlm exc hyp6 7,4 6x3 4x2 hyp7x2
landmarks3={'ABplpappaap','ABplaaaapp','ABarpaapap','ABarpapapa','ABplaaaapa','ABarpaapaa','ABarpapapp','ABplaappaa','ABpraappaa','ABplaapppp','ABpraapppp'};
%tail lm %pvqr pvql p11,p12
landmarks4={'ABprapppaaa','ABplapppaaa','ABplapappa','ABprapappa','Cappppv','Cpppppv'};

%don't use landmarks2
lmtargetnames={landmarks{:},landmarks3{:},landmarks4{:}};

%prepare input emb

%note: last frame in struct generally contains empty point array
for frame  = 1:length(emb)-1
    %remove anisotropy
    emb(frame).finalpoints(:,3) = emb(frame).finalpoints(:,3).*anisotropy;
    %scale to units of reference
    emb(frame).finalpoints = (emb(frame).finalpoints).*(xyscale/xyscaleref);
end

%center input 
cent = [mean(emb(centering_time).finalpoints(:,1)),mean(emb(centering_time).finalpoints(:,2)),mean(emb(centering_time).finalpoints(:,3))];
for i = 1:tend
    emb(i).finalpoints = emb(i).finalpoints - cent;
end


%compute times when both embryos first have 4 cells for alignment

%input
t4cells = tstart;
numcells = size(emb(1).finalpoints,1);
while(numcells < 4)
    t4cells = t4cells + 1;
    numcells = size(emb(t4cells).finalpoints,1); 
end

%reference
t4cellsref = tstartref;
numcellsref = size(embref(1).finalpoints,1);
while(numcellsref < 4)
    t4cellsref = t4cellsref + 1;
    numcellsref = size(embref(t4cellsref).finalpoints,1); 
end

%compute transformations
for frame=t4cells:tend
    
    %keep track of landmark positions and names
    lmpositions1=[];
    lmpositions2=[];
    lmnames={};

    %compute corresponding frame in reference embryo
    a = (frame-t4cells)/(tend-t4cells);
    corresponding_frame=round(t4cellsref+(tendref-t4cellsref)*a);
    %keeping things in bounds (should delete)
    if corresponding_frame == 0
        corresponding_frame = 1;
    end

    %grab snapshots
    names_ref=embref(corresponding_frame).names;
    pos_ref=embref(corresponding_frame).finalpoints;

    names=emb(frame).names;
    pos=emb(frame).finalpoints;

    %get positions of landmark cells in both embryos
    if(length(names_ref)>8)
        
        %for keeping track of landmark names
        c=1;
        
        for i=1:length(lmtargetnames)
            
            lmtarget = lmtargetnames{i};
            
            matchpoint1=[];
            matchpoint2=[];
            
            %find landmark in reference
            for j=1:length(names_ref)
                %if same or if current cell is ancestor of target
                if strcmp(lmtarget,names_ref{j})||~isempty(strfind(lmtarget,names_ref{j}))
                    matchpoint1=pos_ref(j,:);
                end
            end
            %find landmark in input
            for j=1:length(names)
                %if same or if current cell is ancestor of target
                if strcmp(lmtarget,names{j})||~isempty(strfind(lmtarget,names{j}))
                    matchpoint2=pos(j,:);
                end
            end
            
            %if landmark is shared between both embs, collect
            if (~isempty(matchpoint1)&&~isempty(matchpoint2))
                lmnames{c}=lmtarget;
                c=c+1;
                lmpositions1=[lmpositions1;matchpoint1];
                lmpositions2=[lmpositions2;matchpoint2];
            end
        end

        %compute affine transformation -- if no matches, transformation is empty
        
        %https://igl.ethz.ch/projects/ARAP/svd_rot.pdf (Kabsch Algorithm)
        
        %P is input emb; Q is reference where points are column vectors
        P = [lmpositions2]';
        Q = [lmpositions1]';
        centroidP = [mean(P(1,:));mean(P(2,:));mean(P(3,:))];
        centroidQ = [mean(Q(1,:));mean(Q(2,:));mean(Q(3,:))];
        
        X = P;
        Y = Q;
        
        %if wanting to use translation, uncomment
        %X = P - repmat(centroidP,1,size(P,2));
        %Y = Q - repmat(centroidQ,1,size(Q,2));
        
        
        [U,S,V] = svd(X * Y');
        
        %eliminate stretching
        S = eye(size(S));
        %eliminates reflection--pure rotation
        S(end,end) = det(V*U');
        
        %rotation
        R = V * S * U';
        
        T = [0;0;0];
        
        %if wanting to use translation, uncomment
        %T = centroidQ - R * centroidP;
        
        transform2to1 = [R,T];
        %add in extra dimension (so it doesnt break)
        transform2to1 = [transform2to1;0,0,0,1];
        %
        
        %log transform
        alltransforms{frame}=transform2to1;
    end
end

%apply transforms to all timepoints

for frame=1:tend
    transcur=alltransforms{frame};
    %walk forward till you find a transformation
    c=1;
    while (isempty(transcur))
        transcur=alltransforms{frame+c};
        c=c+1;
    end
    %apply affine transformation
    emb(frame).finalpoints=(transcur*[emb(frame).finalpoints,ones(size(emb(frame).finalpoints,1),1)]')';
    emb(frame).finalpoints=emb(frame).finalpoints(:,1:3);
end

transformedemb = emb;
end