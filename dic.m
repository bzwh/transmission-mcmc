clear;

%%
dflag = 0;
if dflag==0
    outpath = './outputs/fmd_sheep/';
    %outpath = './outputs/fmd_sheep/new_sep/';
    %outpath = './outputs/fmd_sheep/new_com/';
elseif dflag==1
    outpath = './outputs/fmd_pigs/both_exp-bestprior/';
    outpath = './outputs/fmd_pigs/both_exp-combined_lat/';
    outpath = './outputs/fmd_pigs/longer-com/';
elseif dflag==2
    outpath = './outputs/asf_pigs/combined/';
    outpath = './outputs/asf_pigs/sep3/';
elseif dflag==3
    outpath = './outputs/vacc_pigs/';
elseif dflag==4
    outpath = './outputs/eble_pigs/old_priors/';
end

pfiles = dir([outpath 'par_' num2str(0) '*']);
lfiles = dir([outpath 'lhd_*']);
ifiles = dir([outpath 'tinf_*']);
efiles = dir([outpath 'latp_' num2str(0) '*']);
pars = {};
lhds = {};
tinf = {};
lats = {};

% Select chains!
II=[2,3]';
%
z = zeros(10000,1);
for i=1:size(II,1)  % Concatenate
  pars{i} = load([outpath pfiles(II(i)).name]);
  lhds{i} = load([outpath lfiles(II(i)).name]);
  tinf{i} = load([outpath ifiles(II(i)).name]);
  lats{i} = load([outpath efiles(II(i)).name]);
end

ParSamp = {};
for i=1:size(II,1)
  ParSamp{i} = [tinf{i} lats{i} pars{i}];
end

%% ==========================================================================
% COMPUTE THE DIC
% Compute the deviance for each sample
PS=[];  % all parameters all chains (tInf latP pars)
logliks = [];
for i=1:size(II,1),
  logliks=[logliks; lhds{i}(:,2)];
  PS=[PS; ParSamp{i}];
end

% Dump mean of parameters to calc lhood in cpp
dlmwrite('./input/dic.txt',mean(PS,1),'delimiter','\n');
%%

L_mean = load('./outputs/llik.txt');
mean_L = mean(logliks);
pd  =  2*(L_mean-mean_L)
DIC = -2*L_mean + 2*pd

%%
D_mean = load('./outputs/llik.txt');

pdalt = 2*var(logliks)
DICalt=-2*D_mean+pdalt

% Compute the effective number of parameters
%pD=Dbar-Dhat;
%==========================================================================


%%

lik_mean - ( 2*(lik_mean-mean_lik) )

%%
d=[0 2 2 3 7 0
0 2 1 2 5 1
0 2 1 2 6 0
0 2 1 2 4 0
0 2 1 2 6 0
0 4 2 3 9 0
0 4 3 4 11 0
0 4 2 3 7 0
0 4 2 3 8 0
0 4 3 4 8 1
1 2 1 2 6 0
1 2 1 2 5 1
1 2 1 2 5 1
1 2 2 3 6 0
1 2 3 4 7 0
1 4 2 3 9 0
1 4 2 3 8 0
1 4 2 3 8 0
1 4 3 4 9 0
1 4 2 3 6 0
2 2 2 3 9 0
2 1 2 3 9 0
3 2 1 2 6 0
3 1 0 1 4 0
4 2 1 2 6 0
4 1 0 1 8 0
5 2 1 2 7 0
5 1 0 1 5 0
6 2 1 2 7 0
6 1 1 2 9 0];

tNeg = d(:,3);
tPos = d(:,4);
tD = d(:,5);
cens = d(:,6);
iType = d(:,2);
rList = d(:,1);












