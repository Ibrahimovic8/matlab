%this function is used to test the time consumed with the radius varying.
%The first experiment in section Comparison to KD-tree for retrieving neighbors
%is done by this function.
function testKDTree(points,index,times,radius)
    [x,dim]=size(points);
    count=10;
    indent=1;
    distMetric='euclidean';
%     distMetric='chebychev';
    IdxMdl=KDTreeSearcher(points,'Distance',distMetric);
%     IdxMdl=KDTreeSearcher(points,'Distance','euclidean');
    radiusV= zeros(count/indent,1);
    timeRecord=zeros(count/indent,1);
    k=0;
    neiNumbers=zeros(count/indent,1);
    for (j=1:indent:count)
        k=k+1;
        tic;
        t1=clock;
        for (i=1:times)
             [tmp]=rangesearch(IdxMdl,points(index,:), j*radius ); 
             radiusV(k)= j*radius;
        end
        toc;
        timeRecord(k)=etime(clock,t1);
        neiNumbers(k)=length(tmp{1});
    end
    


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
    innerMinDimRecorder=zeros(x,2);

    possibleNeibors=0;
    pos=0;
    radiusV1= zeros(count/indent,1);
    timeRecord1=zeros(count/indent,1);
    neiNumbers1=zeros(count/indent,1);
    k=0;
    for (j=1:indent: count)
        CBisearch1(points,sortedValues, j*radius ,innerLowRecorder,innerHighRecorder,innerMinDimRecorder);
        k=k+1;
        tic;
        t1=clock;
%         minnum=0;
%         mindim=1;
%         for i=1:dim
%            low=bisearch(sortedValues,points(index,i)-j*radius);
%            high=bisearch(sortedValues,points(index,i)+j*radius);
%            innerLowRecorder(index,i)=low;
%            innerHighRecorder(index,i)=high;
%            num=high-low;
%            if num<minnum
%                minnum=num;
%                mindim=i;
%            end
%         end;
%         innerMinDimRecorder(index,2)=minnum;
%         innerMinDimRecorder(index,1)=mindim;
        for (i=1:times)
            [possibleNeibors,pos]=findPossibleNeibsByRecorder(points,index,dim,j*radius,sortedLocs,innerLowRecorder,innerHighRecorder,innerMinDimRecorder,possibleNeibors,distMetric);
             radiusV1(k)= j*radius;
        end
        toc;
        neiNumbers1(k)=pos;
        timeRecord1(k)=etime(clock,t1);
    end
    
    plot(radiusV, timeRecord,'-k','linewidth',2);
    hold on
    plot(radiusV1, timeRecord1,'--r','linewidth',2);
    xlabel('Scanning radius','fontsize',15);
    ylabel('Time consumed (second)','fontsize',15);
    legend('KD-tree based algorithm','Algorithm 3');
    titleCaption=['Time comparison under fixed dimension=', num2str(dim),' with scanning radius varing'];
    %title(titleCaption,'fontsize',13);
    axis([radius radius*count 0 1.5* max(max(timeRecord),max(timeRecord1))]);
    
    figure(2);
    hold on;
    plot(radiusV, neiNumbers,'-k','linewidth',2);
    plot(radiusV, neiNumbers1,'--r','linewidth',2);
    xlabel('Scanning radius','fontsize',15);
    ylabel('The numbers of neighbors','fontsize',15);
    titleCaption=['Movement of the total number of neighbors of 10^{th} point with scanning radius expanding'];
    %title(titleCaption,'fontsize',15);
    
    figure(3);
    hold on;
    plot(neiNumbers1, timeRecord,'-k','linewidth',2);
    plot(neiNumbers1, timeRecord1,'--r','linewidth',2);
    xlabel('The numbers of neighbors','fontsize',15);
    ylabel('time consumed (second)','fontsize',15);
    legend('KD-tree based algorithm','Algorithm 3');
    titleCaption=['Time comparison with neighbors increasing'];
    %title(titleCaption,'fontsize',15);

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

%binary searching for a sorted array
function [loc]=bisearch(array, key)
    len=length(array);
    front=1; back=len;    
    while (front<back)
        mid=fix((front+back)/2);             %fix(5.7)=5
        if array(mid)==key;
            loc=mid;
            return;
        end
        if (array(mid)<key)
            front = mid+1;
        else
            back = mid-1;
        end
    end
    loc=front;
end

% int binary_search(double a[], int low, int high, double key) {
%     //double key= target[0];
%     int result=-1;
%     int mid=-1;
%     while (low <high){
%         mid= (low + high)/2;
%         if (a[mid]==key)
%             return mid;
%         else{
%            if (a[mid]>key)
%                high=mid-1;
%            else
%                low=mid+1;
%         }
%     }
%     return low;
% }