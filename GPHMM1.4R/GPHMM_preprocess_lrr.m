function data = GPHMM_preprocess_lrr(data,num1,flag)
%this function is used to impute missing values in genotyping array data.
%Basically, all missing BAF scores will be assigned to 0, which will
%generate exactly the same probablity for every hidden state (so it won't
%cause any bias). For missing LRR scores, you can't just replace them with
%0 as sometimes it will generate significant biased distribution (see
%4240271535_A_tQN, chromosome 2 for example), thus a better way will be
%used median filter.
%flag: if 0 impute with num1;otherwise use median filter
%num1: if 0 value for imputation;otherwise sliding window length
%data: can be 1xN or Nx1

%first make sure the median of LRR signal is 0            
tv = isnan(data);
temp = median(data(~tv));
data(~tv) = data(~tv)-temp;

%then impute NaNs with local LRR signals
len = length(data);
half_w = max(ceil((num1-1)/2),1);
indx1 = find(tv);
%flag == 0, impute with num1
if flag == 0
    data(tv) = num1;
else %impute with median filter, num1 is the length of sliding window
    for i=1:length(indx1)
        indx_w = max(1,indx1(i)-half_w):min(len,indx1(i)+half_w+1);
        data_w = data(indx_w);
        tv_nan = isnan(data_w);
        if any(~tv_nan)
            data(indx1(i)) = median(data_w(~tv_nan));
        else
            data(indx1(i)) = 0; 
        end
    end
end
        
