
%% produce neighbours array for biosemi (or any other EEG data) from EEGLAB.
% usage - run once an EEGLAB dataset is open in EEGLAB (and has electrode
% locations). You must have fieldtrip installed as well. 
% If you increase the "neighbouring_distance" parameter, you'll get more
% neighbours for each electrode. 
%
% Roy Amit, spring 2017

neighbouring_distance_crit=.5;
save_result=1;

if ndims(EEG.data)==3
    data = eeglab2fieldtrip( EEG, 'timelockanalysis', 'none' );
else
    data = eeglab2fieldtrip( EEG, 'preprocessing', 'none' );
end

cfg=data.cfg;
cfg.method='distance'; % alternatively, you can define neighbours by putting 'triangulation' or 'template' , and then "neighbouring_distance_crit" is irrelevant
cfg.channel=1:64;
cfg.elec=data.elec;
cfg.neighbourdist=2;
cfg.channel=1:size(EEG.data,1);
cfg.neighbourdist=neighbouring_distance_crit;
%ft_neighbourplot(cfg)
w=ft_neighbourplot(cfg);
q=w.neighbours


neighbournames={};
neighbournums={};
for elec=1:64
neighbournames{elec}=[];
neighbournums{elec}=[];
for k=1:length(q(elec).neighblabel)
    num= elec_number_by_name(q(elec).neighblabel(k));
    name=q(elec).neighblabel(k);
    if num<65           % if you want the external electrodes, remove this. 
        neighbournames{elec}=[neighbournames{elec} name];
        neighbournums{elec}=[neighbournums{elec} num];
    end
end
end

if save_result
    save('neighboursBIO64','neighbournames','neighbournums')
end