function [ps,clusters]=test_clusters_3D(data_a,data_b,t_crit,tolerance_samps)
% This functions implements Maris & Oostenvald's 2007 method of clustering, for 3D data.
% INPUT: arrays of data (data_a,data_b) which are electrodes*timepoints*(trials or subjects)
% critical t-val; tolerance_samps is number of consecutive nonsignificant
% samples in time which can be tolerated within a cluster. 
% OUTPUT: p values for all clusters; all the clusters, ordered by sum of
% t's (=power). 

% Roy Amit, spring 2017


 
load('neighboursBIO64.mat')
nperms=1000;

data_c=cat(3,data_a,data_b);
real_conds=[zeros(size(data_a,3),1);ones(size(data_b,3),1)];
a=GetSecs;
[elec_tsums_all,isclust,t,edges ,all_clusts,allpowers]=compute_clusters_2d(data_a,data_b,t_crit,tolerance_samps,nmat);
timetaken=GetSecs-a;
q=nmat;
disp('calculating permutations') 
disp(['Taking ' num2str(timetaken*nperms/4) 'seconds'])
disp(repmat('_',1,20))
parfor i=1:nperms
%         if mod(i,round(nperms/20))==0
%             fprintf('-')
%         end
    new_conds=real_conds(randperm(size(data_c,3)));
      data_d=data_c(:,:,new_conds==0);
    data_e=data_c(:,:,new_conds==1);
    
[~,~,~,~ ,~,allpows]=compute_clusters_2d(data_d,data_e,t_crit,tolerance_samps,q);
permpowers(i)=allpows(1);
end

for k=1:length(all_clusts)
    ps(k)=sum(allpowers(k)<permpowers)/nperms;
end
clusters=all_clusts;
