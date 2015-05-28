function [ score ] = GPHMM_postprocess( state,het_seq,w,o,a,sigmaBhet,sigmaBhom,sigmaL)
%GPHMM_POSTPROCESS Summary of this function goes here
%   Detailed explanation goes here
global data_lrr_sep
global data_baf_sep
global data_gc_sep

cn = [0.01 1 2 2 3 3 4 4 4 5 5 5 6 6 6 6 7 7 7 7 ];
Bcn= [0 1 1 2 2 3 3 2 4 4 3 5 5 4 3 6 6 5 4 7 ];
score=[];
for i=1:length(state)
    y=2*(1-w)+cn(state{i})*w;
    Lexp=2*log10(y./2)+o+data_gc_sep{i}*a;
    tempL=normpdf(data_lrr_sep{i},Lexp,sigmaL)./normpdf(0,0,sigmaL);
    tempB=normpdf(data_baf_sep{i},1,sigmaBhom)./normpdf(0,0,sigmaBhom);
    tempB(tempB~=1)=1;
    z=1*(1-w)+Bcn(state{i})*w;
    Bexp=z./y;
    tempB(het_seq{i})=normpdf(data_baf_sep{i}(het_seq{i}),Bexp(het_seq{i}),sigmaBhet)./normpdf(0,0,sigmaBhet);
    score=[score {tempB.*tempL}];
    clear y z Lexp Bexp tempB tempL
end

