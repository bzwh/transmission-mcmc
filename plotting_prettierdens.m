%% 
clear
close
%%
tosave = 1;
if tosave
  vf = 'off';
else
  vf = 'on';
end
%% LOAD THE DATA
dflag = 3;
pflag = dflag;
aflag = 0;
tflag = 0;
II=[1,2,3,4]';
if dflag==0
    outpath = './outputs/fmd_sheep/fixed-sep/';
    pr = load('./input/fmd-sheep/priors_fmdsheep.txt');
    fmdpigs = 0;
    fmd = 1;
elseif dflag==1
    outpath = './outputs/fmd_pigs/both_exp-combined_lat/';
    outpath = './outputs/fmd_pigs/longer-sep/';
    pr = load('./input/fmd-pigs/priors_fmdv_pigs.txt');
    fmdpigs = 1;
    fmd = 1;
elseif dflag==2
    outpath = './outputs/asf_pigs/sep3/';
    %outpath = './outputs/asf_pigs/combined/';
    pr = load('./input/asf-ppc/priors_asf.txt');
    fmdpigs = 1;
    fmd=0;
elseif dflag==3
    outpath = './outputs/vacc_pigs/';
    pr = load('./input/fmd-pigs/priors_fmdv_pigs.txt');
    fmdpigs = 1;
    fmd = 1;
    pflag = 1;
elseif dflag==4
    outpath = './outputs/eble_pigs/old_priors/';
    pr = load('./input/eble-pigs/priors_fmdv_pigs.txt');
    fmdpigs = 1;
    fmd = 1;
    pflag = 1;
elseif dflag>=10
    if aflag
      outpath = ['./outputs/synth-asf/'];
      true_pars = load(['./input/synth-asf/pars_' num2str(dflag) '.txt']);
      pr = load('./input/synth-asf/priors.txt');
      fmdpigs = 1;
      fmd = 1;
      pflag = 2;
    else
      outpath = ['./outputs/synth-fmd/'];
      true_pars = load(['./input/synth-fmd/pars_' num2str(dflag) '.txt']);
      pr = load('./input/synth-fmd/priors.txt');
      fmdpigs = 0;
      fmd = 1;
      pflag = 0;
    end
    tflag = 1;
end

if dflag>=10
    pfiles = dir([outpath 'par_' num2str(dflag) '*']);
else
    pfiles = dir([outpath 'par_' num2str(0) '*']);
end

plah = [];
pars = {};
for i=1:size(II,1)
  pars{i} = load([outpath pfiles(II(i)).name]);
  plah = [plah;pars{i}];
  np = size(pars{i},2);
end

% Parse data
roomID={'room A', 'room B', 'room C', 'room D', 'room E', 'room F'};
routeID={'inoculated', 'within-pen contact', 'between-pen contact'};
switch (size(plah,2))
  case 8
    mflag = 2;
    bflag = 1;
  case 7
    mflag = 2;
    bflag = 0;
  case 6
    mflag = 1;
    bflag = 1;
  case 5
    mflag = 1;
    bflag = 0;
  case 3
    mflag = -1;
    bflag = 0;
  otherwise
    mflag = 0;
    bflag = 0;
end

% Extract priors
p_Ek = pr(1,:);
p_Em = pr(2,:);
p_Ik = pr(3,:);
p_Im = pr(4,:);
p_bt = pr(5,:);

switch (mflag)
  case 2
    lg = {'kE_c','\mu E_c','kE_i','\mu E_i','kI','\mu I','\beta W','\beta_B'};
  case 1
    if (bflag)
      lg = {'kE','\mu E','kI','\mu I','\beta_W','\beta_B'};
    else
      lg = {'kE','\mu E','kI','\mu I','\beta'};
    end
  otherwise
    lg = {'kI','\mu I','\beta'};
end

lw = 1;
cc = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560];
%% Parse parameters
% Extract parameters
if mflag==-1
  kI = plah(:,1); 
  muI= plah(:,2); 
  bW = plah(:,3); 
else
  kEc = plah(:,1);
  muEc= plah(:,2);
  if mflag==2
    kEi = plah(:,3);
    muEi= plah(:,4);
    kI = plah(:,5); 
    muI= plah(:,6); 
    bW = plah(:,7); 
    if (bflag)
      bB = plah(:,8); 
    end
  elseif mflag==1
    kI = plah(:,3); 
    muI= plah(:,4); 
    bW = plah(:,5); 
    if (bflag)
      bB = plah(:,6); 
    end
  else
    kEi = plah(:,1);
    muEi= plah(:,2);
    kI = plah(:,3); 
    muI= plah(:,4); 
    bW = plah(:,5); 
    if (bflag)
      bB = plah(:,6); 
    end
  end
end
%% Plot densities
m=2;
n=3;
fsize = 8;
lw=1;
ccol = [0,0.4470,0.7410];
icol = [0.8500,0.3250,0.0980];
wcol = [0.494,0.184,0.556];
bcol = [0.466,0.674,0.188];
f=figure('outerposition',[10,10,500,400]);

% Latent period shape r-contact b-inoc
subplot(m,n,1); box on
hold on
switch pflag
  case 0
    k = 0:0.1:40;
  case 1
    k = 0:0.1:10;
  case 2
    k = 0:0.1:40;
end
plot(k,gampdf(k,p_Ek(1),p_Ek(2)),'c-','linewidth',lw)
plot(k,pdf(fitdist(kEc,'kernel','support','positive'),k),'Color',ccol,'linewidth',lw)
if (mflag==2)
  plot(k,pdf(fitdist(kEi,'kernel','support','positive'),k),'Color',icol,'linewidth',lw)
end
set(gca,'fontname','arial','fontsize',fsize,'linewidth',1)
xlabel('Latent period shape')
ylabel('Density');

axis manual
if tflag
  vline(true_pars(1),'b');
  vline(true_pars(3),'r');
end

% Mean latent period r-contact b-inoc;
if m==2
  subplot(m,n,4); 
else
  subplot(m,n,2);
end
box on
hold all
switch pflag
  case 0
    i=0:0.05:4;
  case 1
    i=0:0.05:2;
  case 2
    i=0:0.05:10;
end
plot(i,gampdf(i,p_Em(1),p_Em(2)),'c-','linewidth',lw)
plot(i,pdf(fitdist(muEc,'kernel','support','positive'),i),'Color',ccol,'linewidth',lw)
if mflag==2
  plot(i,pdf(fitdist(muEi,'kernel','support','positive'),i),'Color',icol,'linewidth',lw)
end
if tflag
  vline(true_pars(2),'b');
  vline(true_pars(4),'r');
end
%plot(muEc,1e-4*rand([1,size(kEc,1)]),'ro','MarkerSize',2);
set(gca,'fontname','arial','fontsize',fsize,'linewidth',1)
xlabel('Latent period mean')
ylabel('Density');

% Infectious period shape
if m==2
  subplot(m,n,2); box on
else
  subplot(m,n,3); box on
end
hold on
switch pflag
  case 0
    k = 0:0.1:15;
  case 1
    k = 0:0.05:30;
  case 2
    k = 0:0.5:50;
end
plot(k,gampdf(k,p_Ik(1),p_Ik(2)),'c-','linewidth',lw)
plot(k,pdf(fitdist(kI,'kernel','support','positive'),k),'Color',ccol,'linewidth',lw)
set(gca,'fontname','arial','fontsize',fsize,'linewidth',1)
xlabel('Infectious period shape')
ylim manual
if tflag
  vline(true_pars(5),'b');
end


% Mean infectious period
if m==2
  subplot(m,n,5); box on
else
  subplot(m,n,4); box on
end
hold on
switch pflag
  case 0
    i = 0:0.05:25;
  case 1
    i = 0:0.01:10;
  case 2
    i = 0:0.05:10;
end
plot(i,gampdf(i,p_Im(1),p_Im(2)),'c-','linewidth',lw)
plot(i,pdf(fitdist(muI,'kernel','support','positive'),i),'Color',ccol,'linewidth',lw)
%plot(muI,1e-4*rand([1,size(muI,1)]),'ro','MarkerSize',2);
set(gca,'fontname','arial','fontsize',fsize,'linewidth',1)
xlabel('Infectious period mean')
if tflag
  vline(true_pars(6),'b'); %#ok<*UNRCH>
end


% Transmission
if m==2
  subplot(m,n,3); box on
else
  subplot(m,n,5); box on
end
set(gca,'fontname','arial','fontsize',fsize,'linewidth',1)
hold on
if (fmdpigs==0)
  b = 0:0.001:1;
  acol = ccol;
else
  if (fmd==0)
    b = 0:0.1:10;
    acol=ccol;
  else
    b = 0:0.05:5;
    acol = ccol;
  end
end
plot(b,gampdf(b,p_bt(1),p_bt(2)),'c-','linewidth',lw)
if (bflag)
  plot(b,pdf(fitdist(bW,'kernel','support','positive'),b),'Color',wcol,'linewidth',lw)
  plot(b,pdf(fitdist(bB,'kernel','support','positive'),b),'Color',bcol,'linewidth',lw)
else
  plot(b,pdf(fitdist(bW,'kernel','support','positive'),b),'Color',acol,'linewidth',lw)
end
xlabel('Transmission parameter');
axis manual
if tflag
  vline(true_pars(7),'b');
  if bflag
    vline(true_pars(8),'r');
  end
end


% R0
subplot(m,n,6); box on
set(gca,'fontname','arial','fontsize',fsize,'linewidth',1);
hold on
if fmdpigs==0
  r = 0:0.05:5;
else
  if (fmd==0)
    r = 0:0.1:100;
  else
    r = 0:0.5:20;
  end
end

if (bflag)
  plot(r,pdf(fitdist(muI.*bW,'kernel','support','positive'),r),'Color',wcol,'linewidth',lw);
  plot(r,pdf(fitdist(muI.*bB,'kernel','support','positive'),r),'Color',bcol,'linewidth',lw);
else
  plot(r,pdf(fitdist(muI.*bW,'kernel','support','positive'),r),'Color',acol,'linewidth',lw);
end
axis manual
if tflag
  vline(true_pars(8),'b');
end

xlabel('R_0')

%%

rez=300; 
figpos=getpixelposition(f); %dont need to change anything here
resolution=get(0,'ScreenPixelsPerInch'); %dont need to change anything here

set(f,'paperunits','inches','papersize',figpos(3:4)/resolution,'paperposition',[0 0 figpos(3:4)/resolution]); %dont need to change anything here

print(f,[outpath 'transmission-pars'],'-depsc',['-r',num2str(rez)],'-opengl') %save file 
close()







%% LOAD THE DATA
dflag = 0;
pflag = dflag;
aflag = 0;
tflag = 0;
if dflag==0
    outpath = './outputs/fmd_sheep/fixed-com/';
    pr = load('./input/fmd-sheep/priors_fmdsheep.txt');
    fmdpigs = 0;
    fmd = 1;
elseif dflag==1
    outpath = './outputs/fmd_pigs/both_exp-combined_lat/';
    outpath = './outputs/fmd_pigs/both_exp-bestprior/';
    pr = load('./input/fmd-pigs/priors_fmdv_pigs.txt');
    fmdpigs = 1;
    fmd = 1;
elseif dflag==2
    outpath = './outputs/asf_pigs/sep3/';
    %outpath = './outputs/asf_pigs/combined/';
    pr = load('./input/asf-ppc/priors_asf.txt');
    fmdpigs = 1;
    fmd=0;
elseif dflag==3
    outpath = './outputs/vacc_pigs/';
    pr = load('./input/fmd-pigs/priors_fmdv_pigs.txt');
    fmdpigs = 1;
    fmd = 1;
    pflag = 1;
elseif dflag==4
    outpath = './outputs/eble_pigs/old_priors/';
    pr = load('./input/eble-pigs/priors_fmdv_pigs.txt');
    fmdpigs = 1;
    fmd = 1;
    pflag = 1;
elseif dflag>=10
    if aflag
      outpath = ['./outputs/synth-asf/'];
      true_pars = load(['./input/synth-asf/pars_' num2str(dflag) '.txt']);
      pr = load('./input/synth-asf/priors.txt');
      fmdpigs = 1;
      fmd = 1;
      pflag = 2;
    else
      outpath = ['./outputs/synth/' num2str(dflag) '/'];
      true_pars = load(['./input/synth/pars_' num2str(dflag) '.txt']);
      pr = load('./input/synth/priors.txt');
      fmdpigs = 0;
      fmd = 1;
      pflag = 0;
    end
    tflag = 1;
end

if dflag>=10
    pfiles = dir([outpath 'par_' num2str(dflag) '*']);
else
    pfiles = dir([outpath 'par_' num2str(0) '*']);
end
II=[1,2]';
plah = [];
pars = {};
for i=1:size(II,1)
  pars{i} = load([outpath pfiles(i).name]);
  plah = [plah;pars{II(i)}];
  np = size(pars{i},2);
end

% Parse data
roomID={'room A', 'room B', 'room C', 'room D', 'room E', 'room F'};
routeID={'inoculated', 'within-pen contact', 'between-pen contact'};
switch (size(plah,2))
  case 8
    mflag = 2;
    bflag = 1;
  case 7
    mflag = 2;
    bflag = 0;
  case 6
    mflag = 1;
    bflag = 1;
  case 5
    mflag = 1;
    bflag = 0;
  case 3
    mflag = -1;
    bflag = 0;
  otherwise
    mflag = 0;
    bflag = 0;
end

% Extract priors
p_Ek = pr(1,:);
p_Em = pr(2,:);
p_Ik = pr(3,:);
p_Im = pr(4,:);
p_bt = pr(5,:);

switch (mflag)
  case 2
    lg = {'kE_c','\mu E_c','kE_i','\mu E_i','kI','\mu I','\beta W','\beta_B'};
  case 1
    if (bflag)
      lg = {'kE','\mu E','kI','\mu I','\beta_W','\beta_B'};
    else
      lg = {'kE','\mu E','kI','\mu I','\beta'};
    end
  otherwise
    lg = {'kI','\mu I','\beta'};
end

lw = 1;
cc = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560];
%% Parse parameters
% Extract parameters
if mflag==-1
  kI = plah(:,1); 
  muI= plah(:,2); 
  bW = plah(:,3); 
else
  kEc = plah(:,1);
  muEc= plah(:,2);
  if mflag==2
    kEi = plah(:,3);
    muEi= plah(:,4);
    kI = plah(:,5); 
    muI= plah(:,6); 
    bW = plah(:,7); 
    if (bflag)
      bB = plah(:,8); 
    end
  elseif mflag==1
    kI = plah(:,3); 
    muI= plah(:,4); 
    bW = plah(:,5); 
    if (bflag)
      bB = plah(:,6); 
    end
  else
    kEi = plah(:,1);
    muEi= plah(:,2);
    kI = plah(:,3); 
    muI= plah(:,4); 
    bW = plah(:,5); 
    if (bflag)
      bB = plah(:,6); 
    end
  end
end
%% Plot densities
m=2;
n=3;
fsize = 8;
lw=1;
lstyle = '--';
ccol = [0.5,0.5,0.5];%[0,0.4470,0.7410];
icol = [0.5,0.5,0.5];%[0.8500,0.3250,0.0980];
wcol = [0.5,0.5,0.5];%[0.494,0.184,0.556];
bcol = [0.5,0.5,0.5];%[0.466,0.674,0.188];

% Latent period shape r-contact b-inoc
h1=subplot(m,n,1); box on
hold on

switch pflag
  case 0
    k = 0:0.1:10;
  case 1
    k = 0:0.1:10;
  case 2
    k = 0:0.1:50;
end

plot(k,pdf(fitdist(kEc,'kernel','support','positive'),k),'Color',ccol,'linewidth',lw,'linestyle',lstyle)
if (mflag==2)
  plot(k,pdf(fitdist(kEi,'kernel','support','positive'),k),'Color',icol,'linewidth',lw,'linestyle',lstyle)
end
if tflag
  vline(true_pars(1),'b');
  vline(true_pars(3),'r');
end
set(gca,'fontname','arial','fontsize',fsize,'linewidth',1)
xlabel('Latent period shape')

% Mean latent period r-contact b-inoc;
if m==2
  h2=subplot(m,n,4); 
else
  h2=subplot(m,n,2);
end
box on
hold all
switch pflag
  case 0
    i=0:0.05:10;
  case 1
    i=0:0.05:2;
  case 2
    i=0:0.05:10;
end

plot(i,pdf(fitdist(muEc,'kernel','support','positive'),i),'Color',ccol,'linewidth',lw,'linestyle',lstyle)
if mflag==2
  plot(i,pdf(fitdist(muEi,'kernel','support','positive'),i),'Color',icol,'linewidth',lw,'linestyle',lstyle)
end
if tflag
  vline(true_pars(2),'b');
  vline(true_pars(4),'r');
end

%plot(muEc,1e-4*rand([1,size(kEc,1)]),'ro','MarkerSize',2);
set(gca,'fontname','arial','fontsize',fsize,'linewidth',1)
xlabel('Latent period mean')

% Infectious period shape
if m==2
  subplot(m,n,2); box on
else
  subplot(m,n,3); box on
end
hold on
switch pflag
  case 0
    k = 0:0.1:15;
  case 1
    k = 0:0.05:30;
  case 2
    k = 0:0.5:100;
end


plot(k,pdf(fitdist(kI,'kernel','support','positive'),k),'Color',ccol,'linewidth',lw,'linestyle',lstyle)
set(gca,'fontname','arial','fontsize',fsize,'linewidth',1)
xlabel('Infectious period shape')
if tflag
  vline(true_pars(5),'b');
end

% Mean infectious period
if m==2
  subplot(m,n,5); box on
else
  subplot(m,n,4); box on
end
hold on
switch pflag
  case 0
    i = 0:0.05:25;
  case 1
    i = 0:0.01:10;
  case 2
    i = 0:0.05:20;
end

plot(i,pdf(fitdist(muI,'kernel','support','positive'),i),'Color',ccol,'linewidth',lw,'linestyle',lstyle)
%plot(muI,1e-4*rand([1,size(muI,1)]),'ro','MarkerSize',2);
set(gca,'fontname','arial','fontsize',fsize,'linewidth',1)
xlabel('Infectious period mean')
if tflag
  vline(true_pars(6),'b'); %#ok<*UNRCH>
end


% Transmission
if m==2
  subplot(m,n,3); box on
else
  subplot(m,n,5); box on
end
hold on
if (fmdpigs==0)
  b = 0:0.001:1;
  acol = ccol;
else
  if (fmd==0)
    b = 0:0.1:10;
    acol=ccol;
  else
    b = 0:0.05:2;
    acol = ccol;
  end
end

if (bflag)
  plot(b,pdf(fitdist(bW,'kernel','support','positive'),b),'Color',wcol,'linewidth',lw,'linestyle',lstyle)
  plot(b,pdf(fitdist(bB,'kernel','support','positive'),b),'Color',bcol,'linewidth',lw,'linestyle',lstyle)
else
  plot(b,pdf(fitdist(bW,'kernel','support','positive'),b),'Color',acol,'linewidth',lw,'linestyle',lstyle)
end
xlabel('Transmission parameter');
set(gca,'fontname','arial','fontsize',fsize,'linewidth',1)
if tflag
  vline(true_pars(7),'b');
  if bflag
    vline(true_pars(8),'r');
  end
end


% R0
subplot(m,n,6); box on
hold on
if fmdpigs==0
  r = 0:0.05:5;
else
  if (fmd==0)
    r = 0:0.1:100;
  else
    r = 0:0.5:5;
  end
end

if (bflag)
  plot(r,pdf(fitdist(muI.*bW,'kernel','support','positive'),r),'Color',wcol,'linewidth',lw,'linestyle',lstyle);
  plot(r,pdf(fitdist(muI.*bB,'kernel','support','positive'),r),'Color',bcol,'linewidth',lw,'linestyle',lstyle);
else
  plot(r,pdf(fitdist(muI.*bW,'kernel','support','positive'),r),'Color',acol,'linewidth',lw,'linestyle',lstyle);
end

  
if tflag
  vline(true_pars(8),'b');
end


% p1=get(h1,'position');
% p2=get(h2,'position');
% height = p1(2)+p1(4)-p2(2);
% h3 = axes('position',[p2(1) p2(2) p2(3) height],'visible','off');
% h_label=ylabel('Density','visible','on')
set(gca,'fontname','arial','fontsize',fsize,'linewidth',1);
xlabel('R_0')






