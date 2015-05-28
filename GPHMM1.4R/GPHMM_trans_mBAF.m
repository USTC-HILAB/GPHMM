function [data_baf,data_lrr] = GPHMM_trans_mBAF(data_baf,data_lrr,flag1)
%This function is used to transfer BAF signals to mirrored BAF signals.
%When the BAF signals are not symmetric,use imputation for correction.
%data_baf: can be 1xN or Nx1


if nargin==3 && flag1~=0
    %thresold to determine homo bands
    thres_homo = 0.95;
    w_len = 10; %length of one side of a sliding window
    N = length(data_baf);
    indx_all = 1:N;
    
    tv = data_baf>=thres_homo;
    indx_l = find(~tv);
    for i=1:length(indx_l)
        tv_li = indx_all>=max(1,indx_l(i)-w_len);%left indx
        tv_ri = indx_all<=min(indx_l(i)+w_len,N);%right indx
        data_lrr_n = data_lrr(tv&tv_li&tv_ri); %neigbors
        if ~isempty(data_lrr_n)
            %             data_lrr(indx_l(i)) = median(data_lrr_n);
            [temp,indx] = min(abs(data_lrr_n-data_lrr(indx_l(i))));
            data_lrr(indx_l(i)) = data_lrr_n(indx);
        end
    end

    %first transform low homo bands
    tv_l = (data_baf<0.5)&(data_baf<=(1-thres_homo));
    data_baf(tv_l) = 1-data_baf(tv_l);
    %then transform low het bands
    tv_l = (data_baf<0.5)&(data_baf>(1-thres_homo));
    tv_u = (data_baf>0.5)&(data_baf<thres_homo);
    indx_all = 1:N;
    indx_l = find(tv_l);
    for i=1:length(indx_l)
        tv_li = indx_all>=max(1,indx_l(i)-w_len);%left indx
        tv_ri = indx_all<=min(indx_l(i)+w_len,N);%right indx
        data_baf_n = data_baf(tv_u&tv_li&tv_ri); %neigbors
% % %         data_lrr_n = data_lrr(tv_u&tv_li&tv_ri); %neigbors
        if ~isempty(data_baf_n)
            %             data_baf(indx_l(i)) = median(data_baf_imp);
            [temp,indx] = min(abs(data_baf_n-(1-data_baf(indx_l(i)))));
            data_baf(indx_l(i)) = data_baf_n(indx);
% % %             data_lrr(indx_l(i)) = median(data_lrr_n);
        else
            data_baf(indx_l(i)) = 1-data_baf(indx_l(i));
        end
    end%for i=1:length(indx_l)

end %if 0

