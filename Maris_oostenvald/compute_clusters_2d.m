%% in  single electrode
function [elec_tsums_all,isclust,t,all_c_edges,all_clusts,all_powers ]=compute_clusters_2d(data_a,data_b,t_crit,tolerance_samps,nmat,varargin)

% data_a=Rs_i_fa;
% data_b=Rs_i_ra;

if nargin>5
    dim_subjs=varargin{1}
    dim_elecs=varargin{2}
    dim_time=varargin{3}
else

dim_subjs=3;
dim_elecs=1;
dim_time=2;
end

data_a=permute(data_a,[dim_elecs,dim_time,dim_subjs]);
data_b=permute(data_b,[dim_elecs,dim_time,dim_subjs]);
%size(data_a)
t=produce_tmap(data_a,data_b);
if sum(sum(abs(t)>t_crit))==0
    %disp('Failed, no significant points!')
    elec_tsums_all=0;isclust=[];all_c_edges=[];all_clusts=[];all_powers=0;
    return
end
% t_crit=2;
% tolerance_time=5;
nloops=size(data_a,1)*length(t);
%disp(repmat('-',1,20))
%t=t(:,1:1024);
elec_tsums_all={};
totcnt=0;
isclust=zeros(size(data_a,1),length(t));
for elec=1:size(data_a,1)
    clust_tsums=[];
    elec_c_edges=[];
for timepoint=1:length(t)
    totcnt=totcnt+1;
%     if mod(totcnt,round(nloops/20))==0
%         fprintf('-')
%     end
    if ~ isclust(elec,timepoint)
        bads=0;
        if abs(t(elec,timepoint))>t_crit
            
            isclust(elec,timepoint)=1;
            go=1;cnt=0;
            while go
                cnt=cnt+1;
                if (timepoint+cnt)>length(t)
                    if isclust(elec,timepoint+cnt-1)
                        clust_tsum=sum(abs(t(elec,timepoint:timepoint+cnt-1)));
                        clust_tsums=[clust_tsums clust_tsum];
                        clust_edges=[timepoint , timepoint+cnt-1];
                        elec_c_edges=[elec_c_edges  clust_edges'];
                    end
                    break
                end
                if abs(t(elec,timepoint+cnt))>t_crit
                    isclust(elec,timepoint+cnt)=1;
                    bads=0;
                else
                    bads=bads+1;
                    isclust(elec,timepoint+cnt)=1;
                    if bads==tolerance_samps
                        isclust(elec,(timepoint+cnt-tolerance_samps+1):(timepoint+cnt))=0;
                        clust_tsum=sum(abs(t(elec,timepoint:timepoint+cnt-tolerance_samps)));
                        clust_tsums=[clust_tsums clust_tsum];
                        clust_edges=[timepoint , timepoint+cnt-tolerance_samps];
                        elec_c_edges=[elec_c_edges  clust_edges'];
                        break
                    end
                    
                    
                    
                end
            end
        end
    end
end
elec_tsums_all{elec}=clust_tsums;
all_c_edges{elec}=elec_c_edges;
end             

isclust_new=zeros(size(isclust));
clusternum=0;
for elec=1:size(data_a,1)
    
     for clustnum=1:size(all_c_edges{elec},2)
         
        
         
         clust_power=elec_tsums_all{elec}(clustnum);
         clust_elects=[elec];
         clust_borders=all_c_edges{elec}(:,clustnum);
         if ~any(isclust_new(elec,clust_borders(1):clust_borders(2))==1);
             clusternum=clusternum+1;
             isclust_new(elec,clust_borders(1):clust_borders(2))=1;
         [clust_elects,clust_power,clust_borders]=nighbour_clusters(elec,clust_borders(1),clust_borders(2),clust_elects,clust_power,clust_borders);
         all_clusts{clusternum}.elects=clust_elects;
         all_clusts{clusternum}.power=clust_power;
         all_clusts{clusternum}.borders=clust_borders;
         all_powers(clusternum)=clust_power;
         
         
         end
     end
end

    function [clust_elects,clust_power,clust_borders]=nighbour_clusters(elec,tp1,tp2,clust_elects,clust_power,clust_borders)
        
        vecinos=find(nmat(elec,:));
        % disp(elec)
        for k=1:length(vecinos)
            % disp(vecinos)
            if (~ismember (vecinos(k),clust_elects)) && any(isclust(vecinos(k),tp1:tp2)) && ~ any(isclust_new(vecinos(k),tp1:tp2)==1)
                clust_elects=[clust_elects,vecinos(k)];
                
                for q=1:size(all_c_edges{vecinos(k)},2)
                    inds=all_c_edges{vecinos(k)}(1,q):all_c_edges{vecinos(k)}(2,q);
                    if any(ismember(inds,tp1:tp2)) | any(ismember(tp1:tp2,inds))
                        channel_clust_ind=q;
                        break
                    end
                end
                
                new_channel_borders=all_c_edges{vecinos(k)}(:,channel_clust_ind);
                 isclust_new(vecinos(k),new_channel_borders(1):new_channel_borders(2))=1;
                clust_tsumi=sum(abs(t(vecinos(k),new_channel_borders(1):new_channel_borders(2))));
                clust_power=clust_power+clust_tsumi;
                clust_borders=[clust_borders new_channel_borders];
                [clust_elects,clust_power,clust_borders]=nighbour_clusters(vecinos(k),new_channel_borders(1),new_channel_borders(2),clust_elects,clust_power,clust_borders);
            end
        end
        
    end

assert(all(all(isclust==isclust_new)));
[all_powers,ind]=sort(all_powers,'descend');
all_clusts=all_clusts(ind);
end

%%