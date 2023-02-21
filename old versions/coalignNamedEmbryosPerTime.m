function [transformedembs]=coalignNamedEmbryosPerTime(emb_e1,temb1,t4emb1,anisotropy_1,embryostoalign,tsembstoalign,t4embstoalign,anisotropy_2)

%given a template embryo and list of n embryos to be aligned to it
%along with start and end time for computing linear temporal alignment 
%for end frame compute transform for each emb 
%do this anew at each frame


%landmarks
%seam cell landmarks
landmarks={'ABplaaappa';'ABplaaappp';'ABarppaaap';'ABarppapaa';'ABarppapap';'ABplappapa';'ABarppappa';'ABplapapaa';'ABarppappp';'ABarpapppp';'ABarpppaap';'ABarppppaa';'ABarppppap';'ABprappapa';'ABarpppppa';'ABprapapaa';'ABarpppppp'};
%gut
landmarks2={'Ealaad';'Earaad';'Ealaav';'Earaav';'Earpa';'Ealpa';'Earap';'Ealap';'Epraa';'Eplaa';'Earpp';'Ealpp';'Eprap';'Eplap';'Eprpa';'Eplpa';'Eprppa';'Eplppa';'Eplppp';'Eprppp'};
%headlm exc hyp6 7,4 6x3 4x2 hyp7x2
landmarks3={'ABplpappaap','ABplaaaapp','ABarpaapap','ABarpapapa','ABplaaaapa','ABarpaapaa','ABarpapapp','ABplaappaa','ABpraappaa','ABplaapppp','ABpraapppp'};
%tail lm %pvqr pvql p11,p12
landmarks4={'ABprapppaaa','ABplapppaaa','ABplapappa','ABprapappa','Cappppv','Cpppppv'};

lmtargetnames={landmarks{:},landmarks3{:},landmarks4{:}};

transformedembs={};


for e=1:length(embryostoalign) %e is index of each embryo
    for frame=1:length(embryostoalign{e})-1
       
        lmpositions1=[];
        lmpositions2=[];
        lmnames={};
        
        %compute corresponding frame in reference embryo
        
        %t4 is a manual anchor point between two embs, so subtracting it
        %aligns the 'zeros' of both embryos' timesteps. ts is the total
        %length of the runtime (max frame value)? alpha1 is the fraction
        %of the total lifetime of the emb to be aligned that elapses
        %between the anchor point and the 'frame' in question. temb1match
        %is the anchor point in the reference embryo plus the total
        %distance from the anchor point to the end of its capture
        %multiplied by the fraction alpha1, which gives the reference emb
        %anchor point plus the raw movement in frames corresponding to the
        %same fraction of its total lifetime as the difference between the
        %'frame' and anchor point in the emb to be aligned. thus,
        %temb1match is the raw frame that corresponds to the 'frame' in
        %question in the emb to be aligned. it is the same percentage of
        %the reference's total lifetime past the reference's anchor point
        %as the 'frame' is past the emb to be aligned's anchor point. 

        alpha1=(frame-t4embstoalign{e})/(tsembstoalign{e}-t4embstoalign{e});
        temb1match=round(t4emb1+(temb1-t4emb1)*alpha1);
        
        %grab snapshots
        names_1=emb_e1(temb1match).names;
        pos_1=emb_e1(temb1match).finalpoints;
        
        names_2=embryostoalign{e}(frame).names;
        pos_2=embryostoalign{e}(frame).finalpoints;
        
        %get positions of landmark cells in both embryos
        if(length(names_1)>8) %why only good if greater than 8? cell naming?
            c=1;
            for i=1:length(lmtargetnames)
                matchpoint1=[];
                matchpoint2=[];
                
                for j=1:length(names_1)
                    %if same or if current cell is ancestor of target
                    if strcmp(lmtargetnames{i},names_1{j})||~isempty(strfind(lmtargetnames{i},names_1{j}))
                        indmatch1=j;
                        matchpoint1=pos_1(j,:);
                    end
                end
                for j=1:length(names_2)
                    %if same or if current cell is ancestor of target
                    if strcmp(lmtargetnames{i},names_2{j})||~isempty(strfind(lmtargetnames{i},names_2{j}))
                        indmatch2=j;
                        matchpoint2=pos_2(j,:);
                    end
                end
                
                if (~isempty(matchpoint1)&&~isempty(matchpoint2))
                    lmnames{c}=lmtargetnames{i};
                    c=c+1;
                    lmpositions1=[lmpositions1;matchpoint1];
                    lmpositions2=[lmpositions2;matchpoint2];
                end
            end
            
            %scale out anisotropy
            lmpositions1(:,3)=lmpositions1(:,3).*anisotropy_1; %here anisotropy is xyscale/zscale -> 1:1:1
            lmpositions2(:,3)=lmpositions2(:,3).*anisotropy_2;
            
            %compute,log transformation
            transform2to1=[lmpositions1,ones(size(lmpositions1,1),1)]'/[lmpositions2,ones(size(lmpositions2,1),1)]';
            alltransforms{e,frame}=transform2to1;
        end
    end
end

%apply transforms to all timepoints
for i=1:length(embryostoalign)
    e=embryostoalign{i}; % e is now embryo struct
    for t=1:length(e)-1
        transcur=alltransforms{i,t};
        %walk forward till you find a 
        c=1;
        while (isempty(transcur))
            transcur=alltransforms{i,t+c};
            c=c+1;
        end
        %apply affine transformation (using an extra dimension)
        e(t).finalpoints=(transcur*[e(t).finalpoints,ones(length(e(t).finalpoints),1)]')';
        e(t).finalpoints=e(t).finalpoints(:,1:3);
    end
    %log transformed embryo
    transformedembs{i}=e;
end