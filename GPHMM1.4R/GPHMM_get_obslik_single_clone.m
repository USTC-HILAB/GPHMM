function [obslik,condi_probs,condi_probs_fluct] = GPHMM_get_obslik_single_clone(data_baf,data_lrr,Pbs,data_gc,o,a,w,varl,varb,depend_table,normal_prior,varb_homo)
%obslik:
%condi_probs
% magic_num = 0.02;

N = length(data_lrr); %number of data points
w_all = w;%in single clone model, there is only one p

ns = 2; %copy number of stromal cells
Nc = depend_table(:,3)'; %vector of copy numbers of different entries
mus = 0.5;% baf mean of stromal cells
Muc = depend_table(:,4)'; %vector of baf means of different entries
sigmal = sqrt(varl);
sigmab = sqrt(varb);
sigmab_homo = sqrt(varb_homo);
    
%----------------------------------------------------------------------%
%-----------------calculate all of the probablities--------------------%
% fluct_prob = 0.01;
somatic_prob = 0.01; %probability of somatic LOH
temp = 2*Pbs.*(1-Pbs);
Priors = [(1-temp);temp]; % 2*State by Nsample
Num_US = 20; % the number of unique states regulated by global parameters
tv_S = depend_table(:,2)~=0;

Y = w_all(depend_table(tv_S,2)')*ns+(1-w_all(depend_table(tv_S,2)')).*Nc(tv_S);
Z = w_all(depend_table(tv_S,2)')*ns*mus+(1-w_all(depend_table(tv_S,2)')).*Nc(tv_S).*Muc(tv_S);
US_indx = depend_table(tv_S,1);

obslik = zeros(sum(tv_S),N);
condi_probs = zeros(sum(tv_S&(depend_table(:,1)<=Num_US)),N);
condi_probs_fluct = zeros(sum(tv_S&(depend_table(:,1)<=Num_US)),N);

obslik_BAF_Homo = GPHMM_eval_pdf_BAF(data_baf,1,sigmab_homo);%Muc(i) = 1 in all homo bands
for i=1:length(Y)
    if US_indx(i)<=Num_US %normal state
        %LRR
        obslik_LRR = GPHMM_eval_pdf_LRR(data_lrr,log10(Y(i)/2)*2+a*data_gc+o,sigmal);
        %BAF
        temp = find([1 3 8 15] == US_indx(i));
        if ~isempty(temp)
            obslik_BAF_Het = normal_prior(temp)*GPHMM_eval_pdf_BAF(data_baf,Z(i)/Y(i),sigmab);         
        else
            obslik_BAF_Het = GPHMM_eval_pdf_BAF(data_baf,Z(i)/Y(i),sigmab);
        end
        obslik_BAF = somatic_prob*obslik_BAF_Homo+(1-somatic_prob)*(Priors(1,:).*obslik_BAF_Homo+Priors(2,:).*obslik_BAF_Het);
        
        if US_indx(i)==1
            fluct_prob = 0.01;
        else
            fluct_prob = 0.001;
        end
        obslik(i,:) = (1-fluct_prob)*obslik_LRR.*obslik_BAF+fluct_prob*0.1;
        condi_probs(i,:) = (1-fluct_prob)*obslik_LRR.*((1-somatic_prob)*Priors(2,:).*obslik_BAF_Het)./obslik(i,:);
        condi_probs_fluct(i,:) = (fluct_prob*0.1)./obslik(i,:);
    else %the last state 
%   
    end
end

%%%%%%%%%%%%%%%%%%%%%%%

function results = GPHMM_eval_pdf_LRR(data,mu_lrr,Sigma_lrr)

if size(data,1)>size(data,2) %Nx1->1xN
    data = data';
end

if size(mu_lrr,1)>size(mu_lrr,2) %Nx1->1xN
    mu_lrr = mu_lrr';
end

if size(Sigma_lrr,1)>size(Sigma_lrr,2) %Nx1->1xN
    Sigma_lrr = Sigma_lrr';
end

results = normpdf(data,mu_lrr,Sigma_lrr);

%%%%%%%%%%%%%%%%%%%%%%%

function results = GPHMM_eval_pdf_BAF(data,mu_baf,Sigma_baf)

if size(data,1)>size(data,2) %Nx1->1xN
    data = data';
end

if size(mu_baf,1)>size(mu_baf,2) %Nx1->1xN
    mu_baf = mu_baf';
end

if size(Sigma_baf,1)>size(Sigma_baf,2) %Nx1->1xN
    Sigma_baf = Sigma_baf';
end

results = normpdf(data,mu_baf,Sigma_baf);

