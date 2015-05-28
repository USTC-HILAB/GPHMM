function [p_states,num_SNP,aCN,p_s,aCN_C1,aCN_C2,states_pre,state_seq,het_seq] = GPHMM_process_results_new(w,depend_table)
%-----------------------------------------------------
%GPHMM_process_results_new is based on GPHMM_process_results, but can
%automatically handle the results of either multiple clone model or single
%clone model

%------over-all information of the cancer sample------
%p_states: proportions of all hidden states
%num_SNp_all: total number of SNPs investigated
%aCN: include averaged copy numbers with three diffferent ways
%        aCN:
%       waCN:
%       maCN:
%p_s: proportions of chromosomal regions with different scenarios: 1)clone 1 2)clone 2 3)clone 1+2
%p_AACR: proportion of all abberated chromosomal regions identified in both cancer clones
%------sub-clonal information of the cancer sample------
%aCN_C1: averaged copy number for clone 1
%p_ACR_C1: proportion of abberated chromosomal regions identified in clone 1
%aCN_C2: averaged copy number for clone 2
%p_ACR_C2: proportion of abberated chromosomal regions identified in clone 2
%state_pre_C1: predicted MAP states for clone 1
%state_pre_C2: predicted MAP states for clone 2

global gamma_sep
global condi_probs_sep
global condi_probs_fluct_sep

state_seq=[];
het_seq=[];

Num_US = 20;
tv_S = depend_table(:,2)~=0;
US_indx = depend_table(tv_S,1);%uniqe state mapped to states in GPHMM
MS_indx = depend_table(tv_S,2);%multi-clonal status mapped to states in GPHMM
CN_mapping = depend_table(tv_S,3)'; %make sure it's 1xN

%initialize output parameters
p_states = [];
num_SNP = [];
aCN = [];
p_s = [];
if nargout > 4%
    aCN_C1 = [];
    aCN_C2 = [];
    states_pre = [];
end

%initialize intermediate variables
current_num_SNP = 0;
exp_num_states = [];
if nargout > 4
    current_states_pre = [];
end

if length(w)>1 %multiple clones
    w_all = [w(1) w(2) w(1)+w(2)];
    weight = [ w(1)/sum(w)  w(2)/sum(w) 1];
    CN_aCN = [];
    CN_waCN = [];
    CN_maCN = [];
    CN_aCN_C1 =[];
    CN_aCN_C2 = [];
    for j=1:length(w_all)
        tv = depend_table(:,2)==j;
        CN_aCN = [CN_aCN depend_table(tv,3)'];
        CN_waCN = [CN_waCN weight(j)*depend_table(tv,3)'];
        CN_maCN = [CN_maCN weight(j)*depend_table(tv,3)'+(1-weight(j))*2];
        if j~=2
            CN_aCN_C1 = [CN_aCN_C1 depend_table(tv,3)'];
        else
            CN_aCN_C1 = [CN_aCN_C1 2*ones(1,sum(tv))];
        end
        if j~=1
            CN_aCN_C2 = [CN_aCN_C2 depend_table(tv,3)'];
        else
            CN_aCN_C2 = [CN_aCN_C2 2*ones(1,sum(tv))];
        end
    end
end

for i=1:length(gamma_sep) %for the ith chromosome
    post_probs = gamma_sep{i};

    %---handle p_states and num_SNP---
    if isempty(exp_num_states) %initialization
        exp_num_states = zeros(size(post_probs,1),1);
    end
    exp_num_states = exp_num_states+sum(post_probs,2);
    current_num_SNP = current_num_SNP+size(post_probs,2);

    %---handle MAP states---
    if nargout > 4 %output predicted MAP states for both clones
        if length(w)>1 %multiple clones
            [temp,MAP_state] = max(post_probs,[],1);
            state_pre_seg = GPHMM_segment_results(MAP_state,0);
            for k = 1: size(state_pre_seg,1)
                if MS_indx(state_pre_seg(k,3))==1 %clone 1
                    current_states_pre = [current_states_pre;[i,state_pre_seg(k,1:2),...
                        US_indx(state_pre_seg(k,3)),3,1]];
                elseif MS_indx(state_pre_seg(k,3))==2 %clone 2
                    current_states_pre = [current_states_pre;[i,state_pre_seg(k,1:2),...
                        3,US_indx(state_pre_seg(k,3)),2]];
                else % clone 1 and 2,>2*S
                    current_states_pre = [current_states_pre;[i,state_pre_seg(k,1:2),...
                        US_indx(state_pre_seg(k,3)),US_indx(state_pre_seg(k,3)),3]];
                end
            end
        else %single clone
            [temp,MAP_state] = max(post_probs,[],1);
            state_seq=[state_seq {MAP_state}];
            condi_probs=zeros(size(MAP_state));
            condi_probs_fluct=zeros(size(MAP_state));
            for k=1:length(MAP_state)
                condi_probs(k)=condi_probs_sep{i}(MAP_state(k),k);
                condi_probs_fluct(k)=condi_probs_fluct_sep{i}(MAP_state(k),k);
            end
            het_seq=[het_seq {condi_probs>(1-condi_probs-condi_probs_fluct)}];
            state_pre_seg = GPHMM_segment_results(MAP_state,0);
            for k = 1: size(state_pre_seg,1)
                current_states_pre = [current_states_pre;[i,state_pre_seg(k,1:2),...
                    US_indx(state_pre_seg(k,3)),0,1]]; %here use 0 to indicate clone 2 doesn't exist
            end
        end %if length(w)>1
    end %if nargout > 4%output predicted MAP states for both clones

end % for i=1:length(gamma) %for the ith chromosome

%---handle p_states and num_SNP---
current_p_states = exp_num_states/current_num_SNP;
p_states = [p_states current_p_states];
num_SNP = [num_SNP current_num_SNP];

%---handle aCN---
if length(w)>1 %multiple clones
    aCN = [aCN [CN_aCN*current_p_states;CN_waCN*current_p_states;CN_maCN*current_p_states]];
else
    temp = CN_mapping*current_p_states;
    aCN = [aCN [temp;temp;temp]];
end

%---handle p_s---
if length(w)>1 %multiple clones
    normal_tv = (depend_table(tv_S,1)==3)&(depend_table(tv_S,2)==3); %normal state
    w_all = [w(1) w(2) w(1)+w(2)];
    temp = [];
    for j=1:length(w_all)
        tv = (depend_table(tv_S,2)==j)&(depend_table(tv_S,1)<=Num_US)&(~normal_tv);
        temp = [temp;sum(current_p_states(tv))];
    end
    p_s = [p_s [temp;sum(current_p_states(normal_tv));sum(current_p_states(depend_table(tv_S,1)>Num_US))]];%last state
else
    normal_tv = (depend_table(tv_S,1)==3)&(depend_table(tv_S,2)==1); %normal state
    p_s = [p_s [sum(current_p_states((depend_table(tv_S,1)<=Num_US)&(~normal_tv)));0;0;...
        sum(current_p_states(normal_tv));sum(current_p_states(depend_table(tv_S,1)>Num_US))]];
end

%---handle MAP states---
if nargout > 4%output predicted MAP states for both clones
    if length(w)>1 %multiple clones
        %---handle aCN_C1,aCN_C2---
        aCN_C1 = [aCN_C1 CN_aCN_C1*current_p_states];
        aCN_C2 = [aCN_C2 CN_aCN_C2*current_p_states];
    else %single clone
        %---handle aCN_C1,aCN_C2---
        aCN_C1 = [aCN_C1 CN_mapping*current_p_states];
        aCN_C2 = [aCN_C2 0];
    end
    states_pre = [states_pre {current_states_pre}];
end
