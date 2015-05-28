function [LL, prior, transmat, o,a,w,varl,varb,varb_homo,nrIterations] = ...
    GPHMM_EM_Newton_single_clone(stepsize1,init_GPHMM_paras,depend_table, thresh, max_iter,verbose)

%2010/05/15, this new EM algorithm uses univariate method to update global
%parameters sperately and for w, Newton-Raphson method is adopted
global data_lrr_sep
global data_baf_sep
global data_pb_sep
global data_gc_sep
global clamp_thres

previous_loglik = -inf;
converged = 0;
num_iter = 1;
LL = [];
% S = depend_table(end,1); %number of states in the table

if ~iscell(data_baf_sep)
    error('baf data should be stored in cells!');
end

if ~iscell(data_lrr_sep)
    error('lrr data should be stored in cells!');
end

if ~iscell(data_pb_sep)
    error('pb data should be stored in cells!');
end

if ~iscell(data_gc_sep)
    error('GC data should be stored in cells!');
end

prior = init_GPHMM_paras{1};
transmat = init_GPHMM_paras{2};
o = init_GPHMM_paras{3};
a = init_GPHMM_paras{4};
w = init_GPHMM_paras{5};
varl = init_GPHMM_paras{6};
varb = init_GPHMM_paras{7}(1);
varb_homo = init_GPHMM_paras{7}(2);
normal_prior = cell2mat(init_GPHMM_paras{9}(1));
varl_prior = cell2mat(init_GPHMM_paras{9}(2));
varb_prior = cell2mat(init_GPHMM_paras{9}(3));

while (num_iter <= max_iter) && ~converged
    % perform EM algorithm
    [loglik, exp_num_trans, exp_num_visits1,o_u,a_u,w_u,varl_u,varb_u,varb_homo_u] = ...
        GPHMM_compute_ess_single_clone(stepsize1,prior,transmat,o,a,w,varl,varb,depend_table,varl_prior,varb_prior,normal_prior,varb_homo);

    % update parameters
    if init_GPHMM_paras{8}(1)
% %         prior = normalise(exp_num_visits1);
        prior = norm_trans(exp_num_visits1',0)';
    end
    if init_GPHMM_paras{8}(2) && ~isempty(exp_num_trans)
%         clamp_thres = 1-1e-4;
        transmat = norm_trans(exp_num_trans,clamp_thres);
        %===============================================
    end
    if init_GPHMM_paras{8}(3) %update o here, modified by Ao Li
        o = o_u;
    end
    if init_GPHMM_paras{8}(4) %update a here, modified by Ao Li
        a = a_u;
    end
    if init_GPHMM_paras{8}(5) %update w here, modified by Ao Li
        w = w_u;
    end
    if init_GPHMM_paras{8}(6) %update varl here, modified by Ao Li
        varl = varl_u;
    end
    if init_GPHMM_paras{8}(7) %update varb here, modified by Ao Li
        varb = varb_u;
        varb_homo = varb_homo_u;
    end

    if verbose,disp(['o:' num2str(o) ',a:' num2str(a) ',sigmal:' num2str(sqrt(varl)) ',w:' num2str(reshape(w,1,[])) ',sigmab:' num2str(sqrt(varb))]);end
    if verbose, fprintf(1, 'iteration %d, loglik = %f\n', num_iter, loglik); end
    
    num_iter =  num_iter + 1;
    converged = em_converged_m(loglik, previous_loglik, verbose,thresh);
    previous_loglik = loglik;
    LL = [LL loglik];
end
nrIterations = num_iter - 1;

%%%%%%%%%%%%%%%%%%%%%%%

function [loglik, exp_num_trans, exp_num_visits1,o_u,a_u,w_u,varl_u,varb_u,varb_homo_u] = ...
    GPHMM_compute_ess_single_clone(stepsize1,prior,transmat,o,a,w,varl,varb,depend_table,varl_prior,varb_prior,normal_prior,varb_homo)
global data_lrr_sep
global data_baf_sep
global data_pb_sep
global data_gc_sep
global gamma_sep
global condi_probs_sep
global condi_probs_fluct_sep
global tumor_range

numex = length(data_lrr_sep); % each row is a sample
S_all = length(transmat); % number of all states 
exp_num_trans = zeros(S_all,S_all);
exp_num_visits1 = zeros(S_all,1);

%-----------------------E step-----------------------------
gamma_sep = [];
condi_probs_sep = [];
condi_probs_fluct_sep = [];
loglik = 0;
N = 0; % the size of the whole data set

save_memory = 0;

for ex=1:numex %
    if stepsize1 >1 %down_screening
        indx_ds = 1:stepsize1:length(data_lrr_sep{ex});
        obs_baf = data_baf_sep{ex}(indx_ds);
        obs_lrr = data_lrr_sep{ex}(indx_ds);
        obs_pb = data_pb_sep{ex}(indx_ds);
        obs_gc = data_gc_sep{ex}(indx_ds);
        N = N+length(indx_ds);
    else %no ds
        obs_baf = data_baf_sep{ex};
        obs_lrr = data_lrr_sep{ex};
        obs_pb = data_pb_sep{ex};
        obs_gc = data_gc_sep{ex};
        N = N+length(obs_lrr);
    end

    %condi_probs: Pi(G=k|S=j,O)
    [obslik,condi_probs,condi_probs_fluct] = GPHMM_get_obslik_single_clone(obs_baf,obs_lrr,obs_pb,obs_gc,o,a,w,varl+varl_prior,varb+varb_prior,depend_table,normal_prior,varb_homo);
    %     [alpha, beta, gamma, current_ll, xi_summed] = fwdback(prior, transmat, obslik);
    [alpha, gamma, current_ll, beta, xi_summed] = Forward_Backward_Algorithm(prior, transmat, obslik);
    %         [path,current_ll] = viterbi_path(prior, transmat, obslik);
    clear alpha beta obslik;
    loglik = loglik +  current_ll;
    exp_num_trans = exp_num_trans + xi_summed;
    exp_num_visits1 = exp_num_visits1 + gamma(:,1);
    
    if save_memory
        %threshold used to truncate low probabities to zeroes in order to save memory!
        thres = 1e-6;
        gamma(gamma<thres) = 0;
        gamma_sep = [gamma_sep {sparse(gamma)}];
        clear gamma;
        condi_probs(condi_probs<thres) = 0;
        condi_probs_sep = [condi_probs_sep {sparse(condi_probs)}];
        clear condi_probs;
        condi_probs_fluct(condi_probs_fluct<thres) = 0;
        condi_probs_fluct_sep = [condi_probs_fluct_sep {sparse(condi_probs_fluct)}];
        clear condi_probs_fluct;
    else
        gamma_sep = [gamma_sep {gamma}];
        clear gamma;
        condi_probs_sep = [condi_probs_sep {condi_probs}];
        clear condi_probs;
        condi_probs_fluct_sep = [condi_probs_fluct_sep {condi_probs_fluct}];
        clear condi_probs_fluct;
    end
end

%-----------------------M step-----------------------------
%update global parameters
% update w
w_tol = 0.01;
max_iter = 10;
w_u = GPHMM_update_w_single_clone(stepsize1,o,a,w,varl,varb,depend_table,w_tol,max_iter);
if w_u>(1-tumor_range(1)),w_u =(1-tumor_range(1)); end
if w_u<(1-tumor_range(2)),w_u = (1-tumor_range(2)); end

%update o
% o_u = -0.189;
o_u = GPHMM_update_o_single_clone(stepsize1,a,w_u,depend_table);

%update a
a_u = GPHMM_update_a_single_clone(stepsize1,o_u,w_u,depend_table);

%update varl
varl_u = GPHMM_update_varl_single_clone(stepsize1,o_u,a_u,w_u,depend_table);

%update varb
[varb_u,varb_homo_u] = GPHMM_update_varb_single_clone(stepsize1,w_u,depend_table);


%--------------------------------------------------------------------------
function w_u = GPHMM_update_w_single_clone(stepsize1,o,a,w,varl,varb,depend_table,w_tol,max_iter)
global data_lrr_sep
global data_baf_sep
global data_gc_sep
global gamma_sep
global condi_probs_sep
global condi_probs_fluct_sep

numex = length(data_lrr_sep); % each row is a sample
ns = 2; %copy number of stromal cells
Nc = depend_table(:,3)'; %row vector of copy numbers of different entries
mus = 0.5;% baf mean of stromal cells
Muc = depend_table(:,4)'; %row vector of baf means of different entries
Num_US = 20; % the number of unique states regulated by global parameters
tv = (depend_table(:,1)<=Num_US);

iter = 0;
while 1
    Y = w*ns+(1-w)*Nc(tv);
    Z = w*ns*mus+(1-w)*Nc(tv).*Muc(tv);
    Y_d = ns-Nc(tv);
    Z_d = ns*mus-Nc(tv).*Muc(tv);
   
    %first order differential
    ELL_L_D_1 = 0;
    ELL_B_D_1 = 0;
    %second order differential
    ELL_L_D_2 = 0;
    ELL_B_D_2 = 0;

    for ex=1:numex %
        if stepsize1 >1 %down_screening
            indx_ds = 1:stepsize1:length(data_lrr_sep{ex});
            obs_baf = data_baf_sep{ex}(indx_ds);
            obs_lrr = data_lrr_sep{ex}(indx_ds);
            obs_gc = data_gc_sep{ex}(indx_ds);
        else %no ds
            obs_baf = data_baf_sep{ex};
            obs_lrr = data_lrr_sep{ex};
            obs_gc = data_gc_sep{ex};
        end
        post_probs_not_fluct = gamma_sep{ex}(1:length(Y),:).*(1-condi_probs_fluct_sep{ex}(1:length(Y),:));
        post_probs = gamma_sep{ex}(1:length(Y),:).*condi_probs_sep{ex}(1:length(Y),:);
        part1 = o+a*obs_gc-obs_lrr;        
        
        for i=1:length(Y) % 05/28/2010, now only calcualte heter entries in the dependent table
            %LRR
            %first order differential
            ELL_L_D_1 = ELL_L_D_1 + sum(post_probs_not_fluct(i,:).*((2*log10(Y(i)/2)+part1)*Y_d(i)/Y(i)));
            %second order differential
            ELL_L_D_2 = ELL_L_D_2 + sum(post_probs_not_fluct(i,:).*((2/log(10)-(2*log10(Y(i)/2)+part1))*(Y_d(i)/Y(i))^2));
            %BAF
            %first order fifferential
            ELL_B_D_1 = ELL_B_D_1 + sum(post_probs(i,:).*((Z_d(i)*Y(i)-Z(i)*Y_d(i))*(Z(i)-obs_baf.*Y(i))/Y(i)^3));
            %second order differential
            ELL_B_D_2 = ELL_B_D_2 + sum(post_probs(i,:).*((Z_d(i)*Y(i)-Z(i)*Y_d(i))*(Z_d(i)*Y(i)+Y_d(i)*(2*Y(i)*obs_baf-3*Z(i)))/Y(i)^4));
        end
    end

    ELL_L_D_1 = -2/(log(10)*varl)*ELL_L_D_1;
    ELL_B_D_1 = -1/(varb)*ELL_B_D_1;
    ELL_ALL_D_1 = ELL_L_D_1+ELL_B_D_1;

    ELL_L_D_2 = -2/(log(10)*varl)*ELL_L_D_2;
    ELL_B_D_2 = -1/(varb)*ELL_B_D_2;
    ELL_ALL_D_2 = ELL_L_D_2+ELL_B_D_2;

    w_u = w - ELL_ALL_D_1/ELL_ALL_D_2;
    iter = iter+1;
    
    if abs(w_u-w)<w_tol || iter>max_iter
%         disp(['The Newton-Raphson method is done in ' num2str(iter) ' iterations'])
        break;
    else
        w = w_u;
    end
end

%--------------------------------------------------------------------------
function o_u = GPHMM_update_o_single_clone(stepsize1,a,w,depend_table)
global data_lrr_sep
global data_gc_sep
global gamma_sep
global condi_probs_fluct_sep

numex = length(data_lrr_sep); % each row is a sample
ns = 2; %copy number of stromal cells
Nc = depend_table(:,3)'; %row vector of copy numbers of different entries
Num_US = 20; % the number of unique states regulated by global parameters
tv = (depend_table(:,1)<=Num_US);
Y = w*ns+(1-w)*Nc(tv);

o_num = 0;
o_den = 0;
for ex=1:numex %
    if stepsize1 >1 %down_screening
        indx_ds = 1:stepsize1:length(data_lrr_sep{ex});
        obs_lrr = data_lrr_sep{ex}(indx_ds);
        obs_gc = data_gc_sep{ex}(indx_ds);
    else %no ds
        obs_lrr = data_lrr_sep{ex};
        obs_gc = data_gc_sep{ex};
    end
    post_probs_not_fluct = gamma_sep{ex}(1:length(Y),:).*(1-condi_probs_fluct_sep{ex}(1:length(Y),:));
    for i=1:length(Y)
        o_num = o_num + sum(post_probs_not_fluct(i,:).*(obs_lrr-(log10(Y(i)/2)*2+a*obs_gc)));
    end
    o_den = o_den + sum(sum(post_probs_not_fluct));
end
o_u = o_num/o_den;

%--------------------------------------------------------------------------
function a_u = GPHMM_update_a_single_clone(stepsize1,o,w,depend_table)
global data_lrr_sep
global data_gc_sep
global gamma_sep
global condi_probs_fluct_sep

numex = length(data_lrr_sep); % each row is a sample
ns = 2; %copy number of stromal cells
Nc = depend_table(:,3)'; %row vector of copy numbers of different entries
Num_US = 20; % the number of unique states regulated by global parameters
tv = (depend_table(:,1)<=Num_US);
Y = w*ns+(1-w)*Nc(tv);

a_num = 0; %numerator of o estimator
a_den = 0; %denominator of a estimator
for ex=1:numex %
    if stepsize1 >1 %down_screening
        indx_ds = 1:stepsize1:length(data_lrr_sep{ex});
        obs_lrr = data_lrr_sep{ex}(indx_ds);
        obs_gc = data_gc_sep{ex}(indx_ds);
    else %no ds
        obs_lrr = data_lrr_sep{ex};
        obs_gc = data_gc_sep{ex};
    end
    post_probs_not_fluct = gamma_sep{ex}(1:length(Y),:).*(1-condi_probs_fluct_sep{ex}(1:length(Y),:));
    for i=1:length(Y)
        a_num = a_num + sum(post_probs_not_fluct(i,:).*obs_gc.*(obs_lrr-(log10(Y(i)/2)*2+o)));
        a_den = a_den + sum(post_probs_not_fluct(i,:).*obs_gc.*obs_gc);
    end
end
a_u = a_num/a_den;

%--------------------------------------------------------------------------
function varl_u = GPHMM_update_varl_single_clone(stepsize1,o,a,w,depend_table)
global data_lrr_sep
global data_gc_sep
global gamma_sep
global condi_probs_fluct_sep

numex = length(data_lrr_sep); % each row is a sample
ns = 2; %copy number of stromal cells
Nc = depend_table(:,3)'; %row vector of copy numbers of different entries
Num_US = 20; % the number of unique states regulated by global parameters
tv = (depend_table(:,1)<=Num_US);
Y = w*ns+(1-w)*Nc(tv);

varl_num = 0; %numerator of o estimator
varl_den = 0;
for ex=1:numex %
    if stepsize1 >1 %down_screening
        indx_ds = 1:stepsize1:length(data_lrr_sep{ex});
        obs_lrr = data_lrr_sep{ex}(indx_ds);
        obs_gc = data_gc_sep{ex}(indx_ds);
    else %no ds
        obs_lrr = data_lrr_sep{ex};
        obs_gc = data_gc_sep{ex};
    end
    post_probs_not_fluct = gamma_sep{ex}(1:length(Y),:).*(1-condi_probs_fluct_sep{ex}(1:length(Y),:));
    for i=1:length(Y)
        varl_num = varl_num +  sum(post_probs_not_fluct(i,:).*(obs_lrr-(log10(Y(i)/2)*2+a*obs_gc+o)).^2);
    end
    varl_den = varl_den + sum(sum(post_probs_not_fluct));
end
varl_u = varl_num/varl_den;

%--------------------------------------------------------------------------
function [varb_u,varb_homo_u] = GPHMM_update_varb_single_clone(stepsize1,w,depend_table)
global data_baf_sep
global gamma_sep
global condi_probs_sep
global condi_probs_fluct_sep

numex = length(data_baf_sep); % each row is a sample
ns = 2; %copy number of stromal cells
Nc = depend_table(:,3)'; %row vector of copy numbers of different entries
mus = 0.5;% baf mean of stromal cells
Muc = depend_table(:,4)'; %row vector of baf means of different entries
Num_US = 20; % the number of unique states regulated by global parameters
tv = (depend_table(:,1)<=Num_US);
Y = w*ns+(1-w)*Nc(tv);
Z = w*ns*mus+(1-w)*Nc(tv).*Muc(tv);

varb_num = 0;
varb_den = 0;
varb_homo_num = 0;
varb_homo_den = 0;
for ex=1:numex %
    if stepsize1 >1 %down_screening
        indx_ds = 1:stepsize1:length(data_baf_sep{ex});
        obs_baf = data_baf_sep{ex}(indx_ds);
    else %no ds
        obs_baf = data_baf_sep{ex};
    end
    post_probs =  gamma_sep{ex}(1:length(Y),:).*condi_probs_sep{ex}(1:length(Y),:); %post prob. for heter bands
    post_probs_fluct =  gamma_sep{ex}(1:length(Y),:).*condi_probs_fluct_sep{ex}(1:length(Y),:);
    varb_den = varb_den + sum(sum(post_probs));
     for i=1:length(Y)
        %now only use heter bands to estimate varb.Not sure how it works
        varb_num = varb_num + sum(post_probs(i,:).*(obs_baf-Z(i)/Y(i)).^2);
     end    
     post_probs_homo = 1-sum(post_probs+post_probs_fluct);
     tv = obs_baf<1; %only BAF signals < 1 are considered for estimation of variance
     varb_homo_den = varb_homo_den + sum(post_probs_homo(tv));
     varb_homo_num = varb_homo_num + sum(post_probs_homo(tv).*(obs_baf(tv)-1).^2);
end
varb_u = varb_num/varb_den;
varb_homo_u = varb_homo_num/varb_homo_den;
