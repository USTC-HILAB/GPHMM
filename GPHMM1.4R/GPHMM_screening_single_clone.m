function [LL_all,GPHMM_paras,p_states,num_SNP,aCN,p_s] = GPHMM_screening_single_clone(stepsize1,init_GPHMM_paras,depend_table,thres1,max_iter1,verbose1)

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
    init_GPHMM_paras(5) = w_all(:,i);
    init_GPHMM_paras(6) = varl_all(i);
    init_GPHMM_paras(7) = varb_all(i);
    init_GPHMM_paras(8) = indivec_all(i);
    init_GPHMM_paras(9) = otherparas_all(i);
    
    [LL, prior, transmat, o,a,w,varl,varb,varb_homo,nrIterations] = GPHMM_EM_Newton_single_clone(stepsize1,init_GPHMM_paras,depend_table, thres1, max_iter1,verbose1);
    
    LL_all = [LL_all LL(end)];
    GPHMM_paras{1} = [GPHMM_paras{1} {prior}];
    GPHMM_paras{2} = [GPHMM_paras{2} {transmat}];
    GPHMM_paras{3} = [GPHMM_paras{3} {o}];
    GPHMM_paras{4} = [GPHMM_paras{4} {a}];
    GPHMM_paras{5} = [GPHMM_paras{5} {w}];
    GPHMM_paras{6} = [GPHMM_paras{6} {varl}];
    GPHMM_paras{7} = [GPHMM_paras{7} {[varb;varb_homo]}];
    GPHMM_paras{8} = [GPHMM_paras{8} init_GPHMM_paras(8)];
    GPHMM_paras{9} = [GPHMM_paras{9} init_GPHMM_paras(9)];
    
    if nargout >2
        [temp1,temp2,temp3,temp4] = GPHMM_process_results_new(w,depend_table);
        p_states = [p_states temp1];
        num_SNP = [num_SNP temp2];
        aCN = [aCN temp3];
        p_s = [p_s temp4];
    end

    if verbose1
        t = toc;
        disp('--------------- screening report -----------------')
        disp(['run finished with w_0:' num2str(reshape(w_all{i},1,[])) ',w:' num2str(reshape(w,1,[])) ' ,LL:' num2str(LL(end)) ' ,time used:' num2str(t) ' ,o:' num2str(o,'%5.3g')]);
        disp('--------------- screening report -----------------')
        t_all = t_all+t;
    end
end