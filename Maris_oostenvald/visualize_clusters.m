function []=visualize_clusters(data_a,data_b,times,clusters,ps,alpha,chanlocs)

% this function visualized significant clusters found using test_clusters_2D.
% INPUT: data like in test_clusters_2D, times of samples, p values of
% clusters, alpha (type I err thresh)
%(chanlocs) (see topoplot for chanlocs format).
% if no chanlocs are given, the cluster topo maps  are not shown




pixmap=zeros(size(data_a,1),size(data_a,2));






q=sum(ps<alpha);
colors=[1:-(1/q):(1/q)];

cnt=0;
while 1
    cnt=cnt+1;
    p=ps(cnt);
    if p>alpha
        go=0;
        break;
        
    end
    elects=clusters{cnt}.elects;
    for i=1:length(elects)
        tp1=clusters{cnt}.borders(1,i);tp2=clusters{cnt}.borders(2,i);
        pixmap(i,tp1:tp2)=colors(cnt);
    end
    
    clustmap{cnt}.t=mean(mean(clusters{cnt}.borders));
    clustmap{cnt}.map=zeros(64,1);
    clustmap{cnt}.map(elects)=1;
end

figure;subplot(3,1,1);
mean_a=mean(data_a,3);
mean_b=mean(data_b,3);
maxval=max(max(cat(2,mean_a,mean_b)));
minval=min(min(cat(2,mean_a,mean_b)));

imagesc(times,1:size(mean_a,1),mean_a);a=gca;a.CLim=[minval maxval];
subplot(3,1,2);
imagesc(times,1:size(mean_b,1),mean_b);a=gca;a.CLim=[minval maxval];
subplot(3,1,3);
imagesc(times,1:size(pixmap,1),pixmap)
if nargin>6
for i=1:length(clustmap)
    figure;topoplot(clustmap{i}.map,chanlocs(1:64),'conv','off')
    title(['cluster ' num2str(i) 'centered at' num2str(times(round(clustmap{i}.t))) 'ms'])
end
end