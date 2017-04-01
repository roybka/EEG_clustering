%% in  single electrode
function [elec_tsums_all,isclust,t,all_c_edges ]=compute_clusters(data_a,data_b,t_crit,tolerance_time,varargin)

% data_a=Rs_i_fa;
% data_b=Rs_i_ra;

if nargin>4
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
                    break
                end
                if abs(t(elec,timepoint+cnt))>t_crit
                    isclust(elec,timepoint+cnt)=1;
                    bads=0;
                else
                    bads=bads+1;
                    isclust(elec,timepoint+cnt)=1;
                    if bads==tolerance_time
                        isclust(elec,(timepoint+cnt-tolerance_time+1):(timepoint+cnt))=0;
                        clust_tsum=sum(abs(t(elec,timepoint:timepoint+cnt-tolerance_time)));
                        clust_tsums=[clust_tsums clust_tsum];
                        clust_edges=[timepoint , timepoint+cnt-tolerance_time];
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

