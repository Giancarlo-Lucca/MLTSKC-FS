clear all
clc
tic
clear

load('flags.mat');

sh = [0.01,0.1,1,10,100];
sk = [1,2,3,4,5,6,7,8,9,10];
%sk = [3];


target(find(target==-1))=0;
data = (mapminmax(data',0,1))';
oldOptmParameter=struct('alpha_searchrange',[sh],'beta_searchrange',[sh],'gamma_searchrange',[sh],...
    'maxIter',100,'minimumLossMargin',0.01,'outputtempresult',0,'drawConvergence',0,'bQuiet',0);
TSKoptions=struct('k_searchrange',[sk],'h_searchrange',[sh]);

[ BestParameter, BestResult ] = ML_TSKFS_adaptive_validate( data, target, oldOptmParameter,TSKoptions);

toc 



