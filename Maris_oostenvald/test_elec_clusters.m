

nperms=1000;
data_c=cat(3,data_a,data_b);
real_conds=[zeros(size(data_c,3)/2,1);ones(size(data_c,3)/2,1)];
t_crit=1.5;
tolerance_samps=10;



[elec_tsums_all,isclust,t,edges ]=compute_clusters(data_a,data_b,t_crit,tolerance_samps);
disp(repmat('_',1,20))

for i=1:nperms
        if mod(i,round(nperms/20))==0
            fprintf('-')
        end
    new_conds=real_conds(randperm(size(data_c,3)));
    data_d=data_c(:,:,new_conds==0);
    data_e=data_c(:,:,new_conds==1);
    [elec_tsums_all_perms{i},~,~ ,~]=compute_clusters(data_d,data_e,t_crit,tolerance_samps);
end

t_maxi_real=zeros(1,size(data_a,1));
for elec=1:size(data_a,1)
    if ~isempty(elec_tsums_all{elec})
     [t_maxi_real(elec),ind]=max(elec_tsums_all{elec});
     largest_c_edges(:,elec)=edges{elec}(:,ind);
    end
    
    for i=1:nperms
        var=elec_tsums_all_perms{i}(elec);
        if ~isempty(var{1})
            
            t_maxi_perm(elec,i)=max(var{1});
        end
    end
end

for elec=1:size(data_a,1)
    p(elec)=sum(t_maxi_perm(elec,:)>t_maxi_real(elec))/nperms;
end
