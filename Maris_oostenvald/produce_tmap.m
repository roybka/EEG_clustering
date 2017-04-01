function [ t ] = produce_tmap( data_a,data_b )
% compares 2 datasets along the 3rd dimention using a paired-sample ttest.

 [~,~,~,st]=ttest(data_a,data_b,'dim',3);
 t=st.tstat;
end

