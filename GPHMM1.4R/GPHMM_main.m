function [LL_all,GPHMM_paras,best_indx,ds_info] = ...
    GPHMM_main(screening_method,init_GPHMM_paras,depend_table,depend_table_sc,stepsize_ds,thres_EM,verbose1,fn_nosuffix)
%--------------------------------------------------------------------%
%------------------>       version 0.0.1       <---------------------
%--------------------------------------------------------------------%
%05/25/2010 by Ao
%this is the first version of the GPHMM main function,basically everything
%is done here
%--------------------------- screening -------------------------
global clamp_thres
global mc_w
global current_version
global sin_or_mult
global NoSolutionFlag
clamp_thres = 1-1e-5;
mc_w = 0.8;
thres_EM = 1e-4;
ds_info = [];
verbose1 = 0;

if screening_method==1 %use w0=0.5, fastest, not always accurate, or any given parameters
elseif screening_method==2 %d-sampling screening->w0, screen
    %-------------------------------------------------------------------
    %               ---> d-sampling screening <---
    %-------------------------------------------------------------------
    %initialize parameters
    thres_del = 0.006;
    %     thres1_ratio = 0.003/2;

    init_GPHMM_paras = GPHMM_GPHMM_paras(init_GPHMM_paras,depend_table,depend_table_sc,0,sin_or_mult);
    [LL,GPHMM_paras,p_states,num_SNP,aCN,p_s] = GPHMM_screening...
        (stepsize_ds,init_GPHMM_paras,depend_table,depend_table_sc,thres_EM,30,verbose1);
    %-------------------------------------------------------------------
    %               ---> summarize ds results <---
    %-------------------------------------------------------------------
    %get detailed info for further ananlysis
    %CN           0 1 2 2 3 3 4 4 4 5 5 5 6 6 6 6 7 7 7 7 0
    N_genotype = [0 1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 4 4 4 4 0];
    tv_del = depend_table(depend_table(:,2)~=0,3)<1;
    US_indx = depend_table(depend_table(:,2)~=0,1);

    o_all = cell2mat(GPHMM_paras{3});
    a_all = cell2mat(GPHMM_paras{4});
    w_all = cell2mat(GPHMM_paras{5});
    temp = zeros(1,size(w_all,2));
    tv = w_all(2,:)>0;
    temp(tv) = sum(w_all(:,tv),1);
    w_all = [w_all;temp];
    sigmal_all = sqrt(cell2mat(GPHMM_paras{6}));
    sigmab_all = sqrt(cell2mat(GPHMM_paras{7}));
    entropy_all = GPHMM_entropy(p_states);
    p_total_del = sum(p_states(tv_del,:),1);
    DF_bias = (1-1./(N_genotype(US_indx)+1))*p_states;

    %-------------------------------------------------------------------
    %                   ---> find the best model <---
    %-------------------------------------------------------------------
    %porportion of chromosomal regions with total deletion should be less
    %than the given threshold. Also the average copy number can not be
    %larger than a pre-defined maximum for single clone & high normal
    %contamination(version 0.2.0)
    tv = (p_total_del<thres_del) & ((aCN(1,:)<=4.0 | w_all(2,:)>0 )| w_all(1,:)<0.2);
    NoSolutionFlag=false;
    if ~any(tv)
        warning('Can not find a feasible solution with pre-defined criteria!');
        tv=~tv;
        NoSolutionFlag=true;
    end
    [temp,indx] = max(LL.*tv);
    if w_all(2,indx)>0 %MC
        thres1_ratio = 0.003/2;
    else
        thres1_ratio = 0.0006;
    end
    best_indx = indx;
    %correction for possible signal noise
    %     ratio = (log(num_SNP)/2).*(DF_bias(best_indx)-DF_bias)./(LL(best_indx)-LL+eps);
    ratio = [];
    for i=1:length(LL)
        temp1 = LL(best_indx)-LL(i);
        temp2 =(log(num_SNP(i))/2).*(DF_bias(best_indx)-DF_bias(i));
        if abs(temp1)<=1 && abs(temp2)<=0.01 %identical solutions
            ratio = [ratio 0];
        else
          ratio = [ratio temp2/temp1];
        end
    end

    ds_candi1 = (ratio>=thres1_ratio) & tv;
    if any(ds_candi1)
        [temp,best_indx] = max(ratio.*ds_candi1);
    end

    ds_info = [{[p_total_del;LL;aCN;o_all;a_all;w_all;sigmal_all;sigmab_all;...
        p_s;entropy_all;num_SNP]},{[best_indx,indx]},{[mc_w,thres1_ratio,thres_del]},{[ratio;ds_candi1;DF_bias]},{current_version}];

    %-------------------------------------------------------------------
    init_GPHMM_paras = GPHMM_GPHMM_paras(GPHMM_paras,depend_table,depend_table_sc,best_indx,sin_or_mult);
    [LL_all,GPHMM_paras] = GPHMM_screening(1,init_GPHMM_paras,depend_table,depend_table_sc,5*thres_EM,20,verbose1);
    best_indx = 1;
    %----------------------------------------------------------------------
elseif screening_method==3 %screening w/o d-sampling, best performance, most
else
    error('screening method is not recognized!');
end


function GPHMM_paras = GPHMM_GPHMM_paras(init_GPHMM_paras,depend_table,depend_table_sc,best_indx,multi_clone)
%this function is used to initialize/process parameters for GPHMM training,
%best_indx is used to indicate which parameter configuration (ususally
%multiple generated in previous screening procedure) are selected. If
%best_indx equals 0, parameters will be initialized
global clamp_thres
global tumor_range
% % S = sum(depend_table(:,2)~=0); %number of states in the table
% % S_sc = sum(depend_table_sc(:,2)~=0);

GPHMM_paras = cell(1,9);
%parameter initialization
if best_indx == 0
    %---w---, ###variable in grid searching###
    if isempty(init_GPHMM_paras{5}) %B w
        if multi_clone
            w1_all = [0.10 0.15 0.25 0.40 0.30 0.45 0.50 0.50 0.65 0.75];
            w2_all = [0.00 0.00 0.00 0.00 0.20 0.30 0.20 0.35 0.30 0.20];
            o = [-0.6 0];
            w_0 = [];
            o_all = [];
            for i=1:length(w1_all)
                if w1_all(i)+w2_all(i)>0.9
                    w_0 = [w_0 repmat([w1_all(i);w2_all(i)],1,length(o))];
                    o_all = [o_all o];
                else
                    w_0 = [w_0 repmat([w1_all(i);w2_all(i)],1,1)];
                    o_all = [o_all 0];
                end
            end

        else %single clone
            w_all = [0 0 0 0.2 0.2 0.2 0.4 0.6 0.8 0.85 0.9];
            o_all = [-0.6 -0.3 0 -0.6 -0.3 0 0 0 0 0 0];
            % The soultions candidate are chosen by the tumor proportion range
            tv=w_all>=(1-tumor_range(2))&w_all<=(1-tumor_range(1));
            if sum(tv)==0
                %if no solution candidate, find the closest one, maybe two
                tv=(abs(w_all-(1-tumor_range(2)))==min(abs(w_all-(1-tumor_range(2)))));
                tv=tv|(abs(w_all-(1-tumor_range(2)))==min(abs(w_all-(1-tumor_range(1)))));
            end
            w_all=[w_all(tv) 1-(tumor_range(1)+tumor_range(2))/2];
            o_all=[o_all(tv) 0];
            clear tv
            
            w_0 = w_all;
        end
    else %only one fixed (w,lrr) pair from previous searching
        w_0 = init_GPHMM_paras{5};
    end
    N = size(w_0,2);
    GPHMM_paras{5} = mat2cell(w_0,size(w_0,1),ones(1,size(w_0,2)));

    %---o---
    if isempty(init_GPHMM_paras{3}) %B o
        o_0 = o_all;
    elseif length(init_GPHMM_paras{3})==1%fixed lrr baseshift
        o_0 = repmat(init_GPHMM_paras{3},1,N);
    else % multiple (w,lrr) pairs
        o_0 = init_GPHMM_paras{3};
    end
    GPHMM_paras{3} = mat2cell(o_0,size(o_0,1),ones(1,size(o_0,2)));

    %---a---
    if isempty(init_GPHMM_paras{4}) %B a
        a_0 = 0.0;
    else
        a_0 = init_GPHMM_paras{4};
    end
    GPHMM_paras{4} = repmat({a_0},1,N);

    %---pi---
    if isempty(init_GPHMM_paras{1})
        if multi_clone
            for i=1:N
                if w_0(2,i)>0 %multiple clone
                    S = sum(depend_table(:,2)~=0);
                else %single clone
                    S = sum(depend_table_sc(:,2)~=0);
                end
                prior_0 = 1/(S)*ones(S,1); %Nx1
                GPHMM_paras{1} = [GPHMM_paras{1} {prior_0}];
            end %for i=1:N
        else %if multi_clone
            S = sum(depend_table_sc(:,2)~=0);
            prior_0 = 1/(S)*ones(S,1);
            GPHMM_paras{1} = repmat({prior_0},1,N);
        end
    else %if isempty(init_GPHMM_paras{1})
        prior_0 = init_GPHMM_paras{1};
        GPHMM_paras{1} = repmat({prior_0},1,N);
    end

    %---A---
    if isempty(init_GPHMM_paras{2})
        if multi_clone
            for i=1:N
                if w_0(2,i)>0 %multiple clone
                    S = sum(depend_table(:,2)~=0);
                else %single clone
                    S = sum(depend_table_sc(:,2)~=0);
                end
                transmat_0 = norm_trans(ones(S,S),clamp_thres);
                GPHMM_paras{2} = [GPHMM_paras{2} {transmat_0}];
            end %for i=1:N
        else %if multi_clone
            S = sum(depend_table_sc(:,2)~=0);
            transmat_0 = norm_trans(ones(S,S),clamp_thres);
            GPHMM_paras{2} = repmat({transmat_0},1,N);
        end
    else %if isempty(init_GPHMM_paras{2})
        transmat_0 = init_GPHMM_paras{2};
        GPHMM_paras{2} = repmat({transmat_0},1,N);
    end

    %---varl---
    if isempty(init_GPHMM_paras{6}) %B varl
        varl_0 = 0.2^2;
    else
        varl_0 = init_GPHMM_paras{6};
    end
    GPHMM_paras{6} = repmat({varl_0},1,N);

    %---varb---
    if isempty(init_GPHMM_paras{7}) %B varb
        varb_0 = [(0.03^2);(0.03^2)];
    else
        varb_0 = init_GPHMM_paras{7};
    end
    GPHMM_paras{7} = repmat({varb_0},1,N);

    %---indicator vector---
    if isempty(init_GPHMM_paras{8}) %indicator vector: '1' for update '0' fixed
        adj_all = ones(1,7);
    else
        adj_all = init_GPHMM_paras{8};
    end
    GPHMM_paras{8} = repmat({adj_all},1,N);

    if isempty(init_GPHMM_paras{9}) %parameters for observation function and EM algorithm
        other_paras = [{[1 1.7 1.7 1.7]} 0 0]; %normal prior(0,2,4,6 copies), LRR var0, BAF var0
    else
        other_paras = init_GPHMM_paras{9};
    end
    GPHMM_paras{9} = repmat({other_paras},1,N);

else %parse the results from previous screening
    for i = 1:length(best_indx)
        %--pi--
        GPHMM_paras{1} = [GPHMM_paras{1} init_GPHMM_paras{1}(best_indx(i))];
        % % %         S = length(init_GPHMM_paras{1}{best_indx(i)});
        % % %         GPHMM_paras{1} = [GPHMM_paras{1} {1/(S)*ones(S,1)}];
        %--A--
        GPHMM_paras{2} = [GPHMM_paras{2} init_GPHMM_paras{2}(best_indx(i))];
        % % %         S = length(init_GPHMM_paras{2}{best_indx(i)});
        % % %         GPHMM_paras{2} = [GPHMM_paras{2} {norm_trans(ones(S,S),clamp_thres)}];
        %--o--
        GPHMM_paras{3} = [GPHMM_paras{3} init_GPHMM_paras{3}(best_indx(i))];
        %--a--
        GPHMM_paras{4} = [GPHMM_paras{4} init_GPHMM_paras{4}(best_indx(i))];
        % % %         GPHMM_paras{4} = [GPHMM_paras{4} {0}];
        %--w--
        GPHMM_paras{5} = [GPHMM_paras{5} init_GPHMM_paras{5}(best_indx(i))];
        %--var_l--
        GPHMM_paras{6} = [GPHMM_paras{6} init_GPHMM_paras{6}(best_indx(i))];
        % % %         GPHMM_paras{6} = [GPHMM_paras{6} {0.2^2}];
        %--var_b--
        GPHMM_paras{7} = [GPHMM_paras{7} init_GPHMM_paras{7}(best_indx(i))];
        % % %         GPHMM_paras{7} = [GPHMM_paras{7} {[(0.03^2);(0.03^2)]}];
        %--indicator vector--
        GPHMM_paras{8} = [GPHMM_paras{8} init_GPHMM_paras{8}(best_indx(i))];
        %--parameters for observation function and EM algorithm--
        varl = init_GPHMM_paras{6}{best_indx(i)};
        varb = init_GPHMM_paras{7}{best_indx(i)};
        if (sqrt(varl)/0.2>=1.5) && (sqrt(varb(1))/0.03>=1.5)
            %             GPHMM_paras{9} = [GPHMM_paras{9} {[{[1 2 2 2]} 0 0]}];
            GPHMM_paras{9} = [GPHMM_paras{9} {[{[1 1.7 1.7 1.7]} 0 0]}];
        else
            %sometimes when BAF signals are not good, e.g., unsymmetric(GAP
            %data:RMA_468),using 1.4 or 1.5 may lead to obviously wrong annotation.
            %In this case, 2 should be used.
            GPHMM_paras{9} = [GPHMM_paras{9} {[{[1 1.7 1.7 1.7]} 0 0]}];
        end
    end
end
