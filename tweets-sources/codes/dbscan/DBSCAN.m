function Clust = DBSCAN(DistMat,Eps,MinPts)
%A simple DBSCAN implementation of the original paper:
%"A Density-Based Algorithm for Discovering Clusters in Large Spatial
%Databases with Noise" -- Martin Ester et.al.
%Since no spatial access method is implemented, the run time complexity
%will be N^2 rather than N*logN
%**************************************************************************
%Input: DistMat, Eps, MinPts
%DistMat: A N*N distance matrix, the (i,j) element contains the distance
%from point-i to point-j.

%Eps:     A scalar value for Epsilon-neighborhood threshold.

%MinPts:  A scalar value for minimum points in Eps-neighborhood that holds
%the core-point condition.
%**************************************************************************
%Output: Clust
%Clust:  A N*1 vector describes the cluster membership for each point. 0 is
%reserved for NOISE.
%**************************************************************************
%Written by Tianxiao Jiang, jtxinnocence@gmail.com
%Nov-4-2015
%**************************************************************************

%Initialize Cluster membership as -1, which means UNCLASSIFIED
Clust=zeros(size(DistMat,1),1)-1;
ClusterId=1;

%randomly choose the visiting order
VisitSequence=randperm(length(Clust));

for i=1:length(Clust)
    % For each point, check if it is not visited yet (unclassified)
    pt=VisitSequence(i);
    if Clust(pt)==-1
        %Iteratively expand the cluster through density-reachability
        [Clust,isnoise]=ExpandCluster(DistMat,pt,ClusterId,Eps,MinPts,Clust);
        if ~isnoise
            ClusterId=ClusterId+1;
        end
    end
end

end

function [Clust,isnoise]=ExpandCluster(DistMat,pt,ClusterId,Eps,MinPts,Clust)

%region query
seeds=find(DistMat(:,pt)<=Eps);

if length(seeds)<MinPts
    Clust(pt)=0; % 0 reserved for noise
    isnoise=true;
    return
else
    Clust(seeds)=ClusterId;
    %delete the core point
    seeds=setxor(seeds,pt);
    while ~isempty(seeds)
        currentP=seeds(1);
        %region query
        result=find(DistMat(:,currentP)<=Eps);
        if length(result)>=MinPts
            for i=1:length(result)
                resultP=result(i);
                if Clust(resultP)==-1||Clust(resultP)==0 % unclassified or noise
                    if Clust(resultP)==-1 %unclassified
                        seeds=[seeds(:);resultP];
                    end
                    Clust(resultP)=ClusterId;
                end
                
            end
        end
        seeds=setxor(seeds,currentP);
    end
    isnoise=false;
    return 
end
end

