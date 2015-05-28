function [LL_all,GPHMM_paras,p_states,num_SNP,aCN,p_s] = GPHMM_screening(stepsize1,init_GPHMM_paras, depend_table,depend_table_sc,thres1,max_iter1,verbose1)

%---------------------run the algorithm------------------------------
%1xN cell vectors
prior_all = init_GPHMM_paras{1};
transmat_all = init_GPHMM_paras{2};
o_all = init_GPHMM_paras{3};
a_all = init_GPHMM_paras{4};
w_all = init_GPHMM_paras{5};
varl_all = init_GPHMM_paras{6};
varb_all = init_GPHMM_paras{7};
indivec_all = init_GPHMM_paras{8};
otherparas_all = init_GPHMM_paras{9};

LL_all = [];
GPHMM_paras = cell(1,9);
t_all = 0; 
if nargout >2
    p_states = [];
    num_SNP = [];
    aCN = [];
    p_s = [];
end

for i=1:length(o_all)
    if verbose1
        tic
    end
    
    %1x1 cell
    init_GPHMM_paras(1) = prior_all(i);
    init_GPHMM_paras(2) = transmat_all(i);
    init_GPHMM_paras(3) = o_all(i);
    init_GPHMM_paras(4) = a_all(i);
    init_GPHMM_paras(6) = varl_all(i);
    init_GPHMM_paras(7) = varb_all(i);
    init_GPHMM_paras(8) = indivec_all(i);
    init_GPHMM_paras(9) = otherparas_all(i);
    
    if size(w_all{i},1)> 1  && w_all{i}(2)>0
        init_GPHMM_paras(5) = w_all(:,i);
        [LL, prior, transmat, o,a,w,varl,varb,varb_homo,nrIterations] = GPHMM_EM_Newton(stepsize1,init_GPHMM_paras,depend_table,thres1,max_iter1,verbose1);
    else
%         init_GPHMM_paras(5) = {max(0,1-w_all{i}(1))};
        init_GPHMM_paras(5) = {w_all{i}(1)};
        [LL, prior, transmat, o,a,w,varl,varb,varb_homo,nrIterations] = GPHMM_EM_Newton_single_clone(stepsize1,init_GPHMM_paras,depend_table_sc,thres1,max_iter1,verbose1);
    end
        
    LL_all = [LL_all LL(end)];
    GPHMM_paras{1} = [GPHMM_paras{1} {prior}];
    GPHMM_paras{2} = [GPHMM_paras{2} {transmat}];
    GPHMM_paras{3} = [GPHMM_paras{3} {o}];
    GPHMM_paras{4} = [GPHMM_paras{4} {a}];
    if length(w)>1
        GPHMM_paras{5} = [GPHMM_paras{5} {w}];
    else
%         GPHMM_paras{5} = [GPHMM_paras{5} {[max(0,1-w);0]}];
        GPHMM_paras{5} = [GPHMM_paras{5} {[w;0]}];
    end
    GPHMM_paras{6} = [GPHMM_paras{6} {varl}];
    GPHMM_paras{7} = [GPHMM_paras{7} {[varb;varb_homo]}];
    GPHMM_paras{8} = [GPHMM_paras{8} init_GPHMM_paras(8)];
    GPHMM_paras{9} = [GPHMM_paras{9} init_GPHMM_paras(9)];
    
    if nargout >2
        if length(w)>1
            [temp1,temp2,temp3,temp4] = GPHMM_process_results_new(w,depend_table);
            p_states = [p_states temp1];
        else
            [temp1,temp2,temp3,temp4] = GPHMM_process_results_new(w,depend_table_sc);
            tv = depend_table(:,2)~=0;
            temp = zeros(sum(tv),1);
            temp(depend_table(tv,2)==3) = temp1;
            p_states = [p_states temp];
        end
        num_SNP = [num_SNP temp2];
        aCN = [aCN temp3];
        p_s = [p_s temp4];
    end

    if verbose1
        t = toc;
        disp('--------------- screening report -----------------')
        disp(['run ' num2str(i) ' done, w_0:' num2str(reshape(w_all{i},1,[])) ', w:' ...
            num2str(reshape(w,1,[])) ', LL:' num2str(LL(end),'%5.1f') ', time:' num2str(t,'%5.1f') ', o:' num2str(o,'%5.3f')]);
        disp('--------------- screening report -----------------')
        t_all = t_all+t;
    end
end