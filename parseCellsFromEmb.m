function [ cells ] = parseCellsFromEmb( embinfouned,endtime )
%reparse an embryo into a cell data structure Detailed explanation goes here

cells={};
cellcount=1;
for t=1:endtime
   %  if (t<length(embinfouned))&& ~isempty(embinfouned(t))
    for i=1:size(embinfouned(t).finalpoints,1)
        %parse new cell if from div or start
        if(cellcount==366)
        'odd';
        end
        if(embinfouned(t).pred(i)==-1||embinfouned(t-1).suc(embinfouned(t).pred(i),2)~=-1)
            continueparse=true;
            count=[];
            count.exp=[];
            count.con=[];
            count.pos=[];
            count.matchnames={};
            count.divides=false;
            count.birthposition=embinfouned(t).finalpoints(i,:);
            count.name=embinfouned(t).names{i};
            currenti=i;
            currentt=t;
            while(continueparse)
                if( currentt~=endtime&& embinfouned(currentt).suc(currenti,2)~=-1&& embinfouned(currentt).suc(currenti,1)~=-1)
                    count.divides=true;
                    count.sucnames={embinfouned(currentt+1).names{ embinfouned(currentt).suc(currenti,1)},...
                        embinfouned(currentt+1).names{embinfouned(currentt).suc(currenti,2)}};
                    continueparse=false;
                    count.endtime=currentt;
                else if (currentt==endtime||embinfouned(currentt).suc(currenti,1)==-1)
                        continueparse=false;
                        count.endtime=currentt;
                    end
                end
                count.pos=[count.pos;embinfouned(currentt).finalpoints(currenti,:)];
                count.exp=[count.exp;embinfouned(currentt).expression(currenti)];
                count.con=[count.con;embinfouned(currentt).confidence(currenti)];
                count.matchnames={count.matchnames{:},' '};
                currenti=embinfouned(currentt).suc(currenti,1);
                currentt=currentt+1;
            end
            cells{cellcount}=count;
            cellcount=cellcount+1;
        end
   % end
      end
end

end

