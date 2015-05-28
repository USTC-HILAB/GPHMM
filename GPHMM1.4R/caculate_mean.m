function caculate_mean(acn,w,a,o,depend_table, fn,listfile )
global data_gc_sep
% global data_pos_sep
global data_baf_sep
global data_lrr_sep
global data_genotype_sep2
global data_genotype_sep
% fp=fopen('./output.txt','a+');
% if fp==-1
%     error('can not open the output.txt !');
% end
% fprintf(fp,'%s:\n', fn );
% fprintf(fp,'a=%f\n', a );
% fprintf(fp,'o=%f\n', o );
% fprintf(fp,'w=%f\n', w );
% fprintf(fp,'acn=%f\n', acn );
% Nc = depend_table(:,3)'; %row vector of copy numbers of different entries
% Muc = depend_table(:,4)'; %row vector of baf means of different entries
% Num_US = 20; % the number of unique states regulated by global parameters
% tv = (depend_table(:,1)<=Num_US);
% ns = [0.01 1 2 2 3 3 4 4 4 5 5 5 6 6 6 6 7 7 7 7 ];   %copy number of stromal cells
% ns_multiply_mus=[0 1 1 2 2 3 3 2 4 4 3 5 5 4 3 6 6 5 4 7 ]; % added by lyn
% 
% Y=[];
% Z=[];
% 
% for i=1:length(ns)
%     Y = [Y [w*ns(i)+(1-w)*Nc(tv)]'];
%     Z = [Z [w*ns_multiply_mus(i)+(1-w)*Nc(tv).*Muc(tv)]'];
% end
% fprintf(fp,'\tCHANTS:\n');
% for n=[6 16]
%     obs_Y=[];
%     obs_Z=[];
%     for i=1:length(data_genotype_sep{n})
%         obs_Y =[obs_Y Y(data_genotype_sep2{n}(i) , data_genotype_sep{n}(i) )];
%         obs_Z =[obs_Z Z(data_genotype_sep2{n}(i) , data_genotype_sep{n}(i) )];
%     end
%     %     obs_gc = data_gc_sep{n};
%     %     temp=log10(obs_Y./2).*2+a.*obs_gc+o;
%     temp=log10(obs_Y./2).*2+o;
%     m=mean(temp);
%     fprintf(fp,'\t\tchr%d: lrr_mean=%.5f\n',n,m);
%     temp=obs_Z./obs_Y;
%     m=mean(temp);
%     fprintf(fp,'\t\tchr%d: baf_mean=%.5f\n',n,m);
% end
% 
% fprintf(fp,'\tASCAT:\n');
% for n=[6 16]
%     obs_Y=[];
%     obs_Z=[];
%     for i=1:length(data_genotype_sep{n})
%         obs_Y =[obs_Y Y(data_genotype_sep2{n}(i) , data_genotype_sep{n}(i) )];
%         obs_Z =[obs_Z Z(data_genotype_sep2{n}(i) , data_genotype_sep{n}(i) )];
%     end
%     %     nc_map=Nc(data_genotype_sep2{n});
%     %     ns_map=ns(data_genotype_sep2{n});
%     %     acn=(1-w)*mean(nc_map)+w*mean(ns_map);
%     temp=log2(obs_Y./acn)*0.55;
%     m=mean(temp);
%     fprintf(fp,'\t\tchr%d: lrr_mean=%.5f\n',n,m);
%     temp=obs_Z./obs_Y;
%     m=mean(temp);
%     fprintf(fp,'\t\tchr%d: baf_mean=%.5f\n',n,m);
% end
% 
% fprintf(fp,'\tGAP:\n');
% for n=[6 16]
%     if strcmp(fn,'CRL2324_tQN')==1
%         p=1;
%         q=0.45;
%     elseif strcmp(fn,'CRL2324_21pc_Tum_tQN')==1
%         p=0.14;
%         q=0.12;
%     elseif strcmp(fn,'CRL2324_23pc_Tum_tQN')==1
%         p=0.26;
%         q=0.14;
%     elseif strcmp(fn,'CRL2324_30pc_Tum_tQN')==1
%         p=0.25;
%         q=0.12;
%     elseif strcmp(fn,'CRL2324_34pc_Tum_tQN')==1
%         p=0.27;
%         q=0.16;
%     elseif strcmp(fn,'CRL2324_45pc_Tum_tQN')==1
%         p=0.35;
%         q=0.18;
%     elseif strcmp(fn,'CRL2324_47pc_Tum_tQN')==1
%         p=0.42;
%         q=0.24;
%     elseif strcmp(fn,'CRL2324_50pc_Tum_tQN')==1
%         p=0.42;
%         q=0.22;
%     elseif strcmp(fn,'CRL2324_79pc_Tum_tQN')==1
%         p=0.8;
%         q=0.35;
%     else
%         p=0.5;
%         q=0.2;
%         fprintf(fp,'Can not find the corresponding p and q!\n');
%     end
%     fprintf(fp,'\tp=%f\tq=%f\n',p,q);
%     obs_Y=[];
%     obs_Z=[];
%     Y=[];
%     Z=[];
%     for i=1:length(ns)
%         Y = [Y [(1-p)*ns(i)+p*Nc(tv)]'];
%         Z = [Z [(1-p)*ns_multiply_mus(i)+p*Nc(tv).*Muc(tv)]'];
%     end
%     for i=1:length(data_genotype_sep{n})
%         obs_Y =[obs_Y Y(data_genotype_sep2{n}(i) , data_genotype_sep{n}(i) )];
%         obs_Z =[obs_Z Z(data_genotype_sep2{n}(i) , data_genotype_sep{n}(i) )];
%     end
%     obs_gc = data_gc_sep{n};
%     temp=q*log2(obs_Y./2);
%     m=mean(temp);
%     fprintf(fp,'\t\tchr%d: lrr_mean=%.5f\n',n,m);
%     temp=obs_Z./obs_Y;
%     m=mean(temp);
%     fprintf(fp,'\t\tchr%d: baf_mean=%.5f\n',n,m);
% end
% 
% fclose(fp);
% 
% fp=fopen(fn,'w');
% if fp==-1
%     error('can not open output file!(#)');
% end
% fprintf(fp,'Chr\tBAF\tLRR\tGeno1\tGeno2\tGC\n');
% for n=1:22
%     for i=1:length(data_genotype_sep{n})
%         fprintf(fp, '%d\t%f\t%f\t%d\t%d\t%f\n',n,data_baf_sep{n}(i),data_lrr_sep{n}(i),data_genotype_sep{n}(i),data_genotype_sep2{n}(i),data_gc_sep{n}(i) );
%     end
% end
% fclose(fp);
fp=fopen(listfile,'a+');
if fp==-1
    error('can not open output config file!(&)');
end
fprintf(fp,'%s\t%s\tAbnormal Proportion=%f,GC coeff=%f,LRR baseline shift=%f,ACN=%f\n',datestr(clock, 'mmm-dd-yyyy HH:MM:SS'),fn,1-w,a,o,acn);
fclose(fp);
end