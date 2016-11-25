%this function is used to test the time consumed with the dimension increasing.
%The first experiment in section Comparison to KD-tree for retrieving neighbors
%is done by this function.

function testKDTree2(times,radius)
    dimRecorder=zeros(8,1);
    timeRecorder=zeros(8,1);
    timeRecorder1=zeros(8,1);
    j=1;
    distMetric='euclidean';
%     distMetric='chebychev';
    for (i=10:10:80)
        j=j+1;
        dimRecorder(j)=i;
        points=rand(100000,i);
        [x,dim]=size(points);
%         IdxMdl=KDTreeSearcher(points,'Distance','chebychev');
       IdxMdl=KDTreeSearcher(points,'Distance',distMetric);   
        tic;
        t1=clock;
        for (i=1:times)
             rangesearch(IdxMdl,points(10,:), radius ); 
        end
        toc;
        timeRecorder(j)=etime(clock,t1);
        
        
        sortedValues= zeros(x,dim);
        %global sortedLocs;
        sortedLocs=zeros(x,dim);
        for i=1:dim
            [sortedV,sortedLoc]= sort(points(:,i));
            sortedValues(:,i)= sortedV;
            sortedLocs(:,i)=sortedLoc;
        end

        innerLowRecorder=zeros(x,dim);
        innerHighRecorder=zeros(x,dim);
        innerMinDimRecorder=zeros(x,dim);

        possibleNeibors=0;
        pos=0;
        
        CBisearch1(points,sortedValues, radius ,innerLowRecorder,innerHighRecorder,innerMinDimRecorder);
        tic;
        t1=clock;
        for (i=1:times)
            [possibleNeibors,pos]=findPossibleNeibsByRecorder(points,10,dim,radius,sortedLocs,innerLowRecorder,innerHighRecorder,innerMinDimRecorder,possibleNeibors,distMetric);
        end
        toc;
        timeRecorder1(j)=etime(clock,t1);
    
    end
    

    plot(dimRecorder, timeRecorder,'-k','linewidth',2);
    hold on
    plot(dimRecorder, timeRecorder1,'--r','linewidth',2);
    xlabel('dimension','fontsize',15);
    ylabel('time consumed (second)','fontsize',15);
    legend('KD-tree based algorithm','Algorithm 3');
    %titleCaption=['Time comparison under scanning radius=', num2str(radius),' with dimension varing'];
    %title(['\epsilon=' num2str(radius)]);
    axis([10 80 0 1]);

end


function [possibleNeibors,num]=findPossibleNeibsByRecorder(points,i,dim,radius,sortedLocs,LowRecorder,HighRecorder,MinDimRecorder,possibleNeibors,distMetric)
    %把mindim维度上的所有可能点做为possibleNeibors
    %begin:找出在所有维度上距离小于radius的点
    %end:找出在所有维度上距离小于radius的点
    mindim=MinDimRecorder(i,1);
    minnum=MinDimRecorder(i,2);
    elowLoc=LowRecorder(i,mindim);
    ehighLoc=HighRecorder(i,mindim);
    possNeibs1=sortedLocs(elowLoc:ehighLoc,mindim);
    
%     dists=pdist2(points(i,:),points(possNeibs1,:),'chebychev');
%     possNeibs2=find (dists<=radius);
%     possNeibs1=possNeibs1(possNeibs2);
    for (j=1:dim)
        if (j~=mindim) 
            cv=points(i,j);
            possNeibs2=find ((cv-radius<=points(possNeibs1,j))& (cv+radius>=points(possNeibs1,j)));
            possNeibs1=possNeibs1(possNeibs2);
        end
    end
     if (strcmp(distMetric,'Euclidean'))
        dists=pdist2(points(i,:), points(possNeibs1,:));
        poss=find(dists<=EPS);
        possNeibs1=possNeibs1(poss);
    end
    num=length(possNeibs1);
    possibleNeibors(1:num)=possNeibs1;
end