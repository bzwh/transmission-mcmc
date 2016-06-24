%% 
clear
%% LOAD THE DATA
% Load the challenge experiment data
df = 0;
%outpath = './output-pigs/';%exp2-6custard/';
%outpath = 'd:/ben/projects/transmission_cpp/output-sheep/';
outpath = 'd:/ben/projects/transmission_cpp/outputs/fmd_sheep1/';
outpath = 'd:/ben/projects/transmission_cpp/outputs/asf_pigs/';

fmdpigs = 1;
fmd=0;

pfiles = dir([outpath 'par_' num2str(df) '*']);
bfiles = dir([outpath 'brn_' num2str(df) '*']);
brnflag=1;
plah = [];
blah = [];
pars = {};
parb = {};
for i=1:size(pfiles,1)
  pars{i} = load([outpath pfiles(i).name]);
  np = size(pars{i},2);
  if brnflag
    parb{i} = load([outpath bfiles(i).name]);
    parb{i}=parb{i}(:,end-np+1:end);
    blah = [blah;parb{i}];
  end
  %plah = [plah;pars{i}];
end
% Load the priors used
if (fmd==1)
  if(fmdpigs==1)
    pr = load('d:/ben/projects/transmission_cpp/input/priors_fmdv_pigs.txt');
  else
    pr = load('d:/ben/projects/transmission_cpp/input/fmdv-sheep/priors_fmdsheep.txt');
  end
else
  pr = load('d:/ben/projects/transmission_cpp/input/asfv-ppc/priors_asf.txt');
end

% Parse data
roomID={'room A', 'room B', 'room C', 'room D', 'room E', 'room F'};
routeID={'inoculated', 'within-pen contact', 'between-pen contact'};
switch (size(pars{1},2))
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

%%  plot burn-in traces too?
lw = 2;
II=[1:4]';
plah = [];

h=figure('outerposition',[-900 0 450 1600],'PaperPositionMode','auto');
for i=1:size(pars{1},2)
  subplot(size(pars{1},2),1,i)
  hold all
  %for cno=1:size(pars,2)
  for j=1:size(II);
    cno=II(j);
    aplot=plot([parb{cno}(:,i); pars{cno}(:,i)],'linewidth',lw);
    aplot.Color(4)=0.5/cno;
  end
  xlabel(lg(i))
  set(gca,'fontname','arial','fontsize',8,'linewidth',1)
end
legend('0','1','2','3','Location','northwest');
print(h,'-dpng','-r0',[outpath 'pars-burn.png'])

%% Parse parameters
for i=1:size(II,1)
  plah = [plah;pars{II(i)}];
end
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
%% Parameter traces
lw = 2;
h=figure('outerposition',[-450 0 450 1600],'PaperPositionMode','auto');

for i=1:size(pars{1},2)
  subplot(size(pars{1},2),1,i)
  hold all
  for cno=1:size(pars,2)
    aplot=plot(pars{cno}(:,i),'linewidth',lw);
    aplot.Color(4)=0.5/cno;
  end
  xlabel(lg(i))
  set(gca,'fontname','arial','fontsize',8,'linewidth',1)
end
legend('0','1','2','3','Location','northwest');
print(h,'-dpng','-r0',[outpath 'pars-traces.png'])


%% Plot densities
fsize = 8;
tflag = 0;
lw=2;
ccol = [0,0.4470,0.7410];
icol = [0.8500,0.3250,0.0980];

figure('units','normalized','outerposition',[0 0 1 1])
if (mflag>0)
% Latent period shape r-contact b-inoc
subplot(2,3,1)
hold on
if (fmdpigs)
  k = 0:0.5:50;
else
  k = 0:0.1:10;
end
plot(k,pdf(fitdist(kEc,'kernel','support','positive'),k),'Color',ccol,'linewidth',lw)
if (mflag==2)
  plot(k,pdf(fitdist(kEi,'kernel','support','positive'),k),'Color',icol,'linewidth',lw)
end
plot(k,gampdf(k,p_Ek(1),p_Ek(2)),'c-','linewidth',lw)
set(gca,'fontname','arial','fontsize',fsize,'linewidth',1)
xlabel('Latent period shape')
if mflag==2
  legend('Contact','Inoculated','Prior');
end
if tflag
  vline(true_pars(1),'Color',ccol);
  vline(true_pars(3),'Color',icol);
end

% Mean latent period r-contact b-inoc;
subplot(2,3,2)
hold all
m=0:0.1:10;
plot(m,pdf(fitdist(muEc,'kernel','support','positive'),m),'Color',ccol,'linewidth',lw)
if mflag==2
  plot(m,pdf(fitdist(muEi,'kernel','support','positive'),m),'Color',icol,'linewidth',lw)
end
plot(m,gampdf(m,p_Em(1),p_Em(2)),'c-','linewidth',lw)
%plot(muEc,1e-4*rand([1,size(kEc,1)]),'ro','MarkerSize',2);
set(gca,'fontname','arial','fontsize',fsize,'linewidth',1)
xlabel('Latent period mean')
if tflag
  vline(true_pars(2),'Color',ccol);
  vline(true_pars(4),'Color',icol);
end
end

% Transmission
subplot(2,3,3)
hold on
if (fmdpigs==0)
  b = 0:0.005:1;
else
  b = 0:0.1:10;
end
plot(b,pdf(fitdist(bW,'kernel','support','positive'),b),'Color',icol,'linewidth',lw)
if (bflag)
  plot(b,pdf(fitdist(bB,'kernel','support','positive'),b),'Color',ccol,'linewidth',lw)
end
plot(b,gampdf(b,p_bt(1),p_bt(2)),'c-','linewidth',lw)
%plot(b,unifpdf(b,p_bt(1),p_bt(2)),'c-','linewidth',lw)
if (fmd==0)
  legend('Within-pen','Between-pen','Prior');
end
xlabel('Transmission parameter');
set(gca,'fontname','arial','fontsize',fsize,'linewidth',1)
if tflag
  vline(true_pars(7),'Color',ccol);
  vline(true_pars(8),'Color',icol);
end

% Infectious period shape
subplot(2,3,4)
hold on
if (fmdpigs)
  k = 0:0.5:50;
else
  k = 0:0.1:10;
end
plot(k,pdf(fitdist(kI,'kernel','support','positive'),k),'b','linewidth',lw)
plot(k,gampdf(k,p_Ik(1),p_Ik(2)),'c-','linewidth',lw)
set(gca,'fontname','arial','fontsize',fsize,'linewidth',1)
xlabel('Infectious period shape')
if tflag
  vline(true_pars(5),'b');
end

% Mean infectious period
subplot(2,3,5)
hold on
m = 0:0.5:25;
plot(m,pdf(fitdist(muI,'kernel','support','positive'),m),'b','linewidth',lw)
plot(m,gampdf(m,p_Im(1),p_Im(2)),'c-','linewidth',lw)
%plot(muI,1e-4*rand([1,size(muI,1)]),'ro','MarkerSize',2);
set(gca,'fontname','arial','fontsize',fsize,'linewidth',1)
xlabel('Infectious period mean')
if tflag
  vline(true_pars(6),'b'); %#ok<*UNRCH>
end

% R0
subplot(2,3,6)
hold on
if fmdpigs==0
  r = 0:0.005:5;
else
  r = 0:0.1:50;
end
plot(r,pdf(fitdist(muI.*bW,'kernel','support','positive'),r),'Color',icol,'linewidth',lw);
if (bflag)
  plot(r,pdf(fitdist(muI.*bB,'kernel','support','positive'),r),'Color',ccol,'linewidth',lw);
end
set(gca,'fontname','arial','fontsize',fsize,'linewidth',1);
xlabel('R_0')
%suptitle(num2str(df));
saveas(gcf,[outpath 'pars-dens.png'])

%% Load infection time data-------------------------------------------------------------------------
tflag=0;
ifiles = dir([outpath 'tinf_*']);
tInf = [];
tinf = {};
nb = 50;
for i=1:size(pfiles,1)
  tinf{i} = load([outpath ifiles(i).name]);
  tInf = [tInf;tinf{i}];
end
if (tflag)
  tinf_true = load(['../transmission_gen/data/' num2str(df) '-tinf.txt']);
end
n = size(tinf{1},1);
na = size(tinf{1},2);
%% Infection times - Traces
figure('units','normalized','outerposition',[0 0 1 1])
nsep = 5;
for i=1:na
  subplot(nsep,ceil(na/5),i)
  hold on
  for cno=1:size(tinf,2)
    aplot=plot(tinf{cno}(:,i));
    aplot.Color(4)=1/cno;
  end
end
suptitle('Infection Time');
saveas(gcf,[outpath 'tinf-traces.png'])
%% Infection times of contact transmission - densities
figure('units','normalized','outerposition',[0 0 1 1])
nsp = 5;
for i=1:na
  subplot(nsp,ceil(na/nsp),i)
  hold on
  [h,c] = hist(tInf(:,i),nb);
  plot(c,smooth(h/size(tInf,1)),'linewidth',2)
  for cno=1:size(tinf,2)
    [h,c] = hist(tinf{cno}(:,i),nb);
    l = plot(c,smooth(h/n));
    l.Color(4)=0.5;
  end
  if (tflag)
    vline(tinf_true(i),'-');
  end
end
saveas(gcf,[outpath 'tinf-hist.png'])

%% Load latent period data
lfiles = dir([outpath 'latp_' num2str(df) '*']);
lat = {};
Lat = [];
for i=1:size(lfiles,1)
  lat{i} = load([outpath lfiles(i).name]);
  Lat = [Lat;lat{i}];
end
%% Latent period traces
figure('units','normalized','outerposition',[0 0 1 1])
n = size(lat{1},1);
for i=1:size(lat{1},2)
  subplot(ceil(size(lat{1},2)/nsp),nsp,i)
  hold on
  for cno=1:size(lat,2)
    aplot = plot(lat{cno}(:,i));
    aplot.Color(4)=1/cno;
  end
end
suptitle('Latent period');
saveas(gcf,[outpath 'latp-traces.png'])
%% Latent period durations
figure('units','normalized','outerposition',[0 0 1 1])
nb = 50;
nsp = 5;
for i=1:size(lat{1},2)
  %subplot(4,ceil(size(lat{1},2)/4),i)
  subplot(ceil(size(lat{1},2)/nsp),nsp,i)
  hold on
  [h,c] = hist(Lat(:,i),nb);
  plot(c,smooth(h/size(Lat,1)),'linewidth',2);
  for cno=1:size(lat,2)
    [h,c] = hist(lat{cno}(:,i),nb);
    l = plot(c,smooth(h/n));
    l.Color(4)=0.5;
  end
end
saveas(gcf,[outpath 'latp-hist.png'])

%% Load infectious period data - not actually inferred = (t_E - t_I - E)
ifiles = dir([outpath 'infp_' num2str(df) '*']);
infp = {};
Infp = [];
for i=1:size(pfiles,1)
  infp{i} = load([outpath ifiles(i).name]);
  Infp = [Infp;infp{i}];
end
n = size(infp{1},1);
nsp = 5;
%% Infectious period traces
figure('units','normalized','outerposition',[0 0 1 1])
for i=1:size(infp{1},2)
  subplot(ceil(size(infp{1},2)/nsp),nsp,i)
  hold on
  for cno=1:size(infp,2)
    plot(infp{cno}(:,i))
  end
end
suptitle('Infectious period');
saveas(gcf,[outpath 'infp-traces.png'])
%% Infectious period durations
figure('units','normalized','outerposition',[0 0 1 1])
nb = 25;
for i=1:size(infp{1},2)
  subplot(ceil(size(infp{1},2)/nsp),nsp,i)
  hold on
  [h,c] = hist(Infp(:,i),nb);
  plot(c,smooth(h/size(Infp,1)),'linewidth',2);
  for cno=1:size(infp,2)
    [h,c] = hist(infp{cno}(:,i),nb);
    plot(c,smooth(h/n));
  end
end
%suptitle('Infectious period');
saveas(gcf,[outpath 'infp-hist.png'])

%% correlation plots
np = size(plah,2);
nb = 50;

switch (np)
  case 3
    lbls = {'kI','\mu I','\beta'};
  case 5
    lbls = {'kE','\mu E','kI','\mu I','\beta'};
  case 7
    lbls = {'kEc','\mu Ec','kEi','\mu Ei','kI','\mu I','\beta'};
  case 8
    lbls = {'kEc','\mu Ec','kEi','\mu Ei','kI','\mu I','\beta_W','\beta_B'};
end
np = size(plah,2);
figure('units','normalized','outerposition',[0 0 1 1])
for i=1:size(plah,2)
  for j=1:size(plah,2)
    if i==j
      subplot(np,np,(j-1)*np+i)
      [h,c]=hist(plah(:,i),nb);
      plot(c,h/size(plah,1));
      xlabel(lbls{i})
    elseif i>j
      subplot(np,np,(j-1)*np+i)
      hh=hist3(plah(:,[j,i]),[nb,nb]);
      h = hh';
      h(size(hh,1) + 1, size(hh,2) + 1) = 0;
      h(h==0) = NaN;
      yb = linspace(min(plah(:,i)),max(plah(:,i)),size(hh,1)+1);
      xb = linspace(min(plah(:,j)),max(plah(:,j)),size(hh,1)+1);
      p=pcolor(xb,yb,h);
      set(p,'Edgecolor','none');
    end
  end
end
saveas(gcf,[outpath 'pars-corr.png'])

%% Autocorrelation plots for model parameters
figure('units','normalized','outerposition',[0 0 1 1])
nsp = 2;
tn = 1;
for i=1:size(plah,2)
  subplot(nsp,ceil(size(plah,2)/nsp),i)
  autocorr(plah(1:tn:end,i),100)
end
saveas(gcf,[outpath 'pars-autocor.png'])






