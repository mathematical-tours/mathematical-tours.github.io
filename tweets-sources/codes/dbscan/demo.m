%A demo to cluster a synthetic two spirals dataset 
%%Written by Tianxiao Jiang, jtxinnocence@gmail.com
%Nov-4-2015

%generate two spirals synthetic data
%data=twospirals(200,360,50,1.5,15);

data = twospirals(200, 360, 0, 1.5);
% data=corners(500);

% figure
% scatter(data(:,1),data(:,2),'filled','markerfacecolor',[0.8,0.8,0.8]);


%%
%using euclidean distance measurement
distmat=zeros(size(data,1),size(data,1));

for i=1:size(data,1)
    for j=i:size(data,1)
        distmat(i,j)=sqrt((data(i,1:2)-data(j,1:2))*(data(i,1:2)-data(j,1:2))');
    end
end

for i=1:size(data,1)
    for j=i:size(data,1)
        distmat(j,i)=distmat(i,j);
    end
end

% k_dist=zeros(size(data,1),1);
% figure
% for k=3:5
%     for i=1:size(data,1)
%         tmp=sort(distmat(:,i),'ascend');
%         k_dist(i)=tmp(k);
%     end
%     hold on
%     plot(1:size(data,1),k_dist);
% end
%%
Eps=0.5;
MinPts=4;
Clust = DBSCAN(distmat,Eps,MinPts);

x=data(:,1);
y=data(:,2);
%plot different clusters
figure
for i=1:max(Clust)
    hold on
    scatter(x(Clust==i),y(Clust==i),'filled');
end

%plot the noise
hold on
scatter(x(Clust==0),y(Clust==0),'*')