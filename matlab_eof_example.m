% PURPOSE:
% This script looks for explained variance in simulated melt using a EOF rotation scheme.

%% Prepare simulated data
clear all
close all
clc

homedir = '/Users/karllapo/Dropbox/';
workdir = 'C:\Users\Karl\Dropbox\';
atmosdir = '/home/disk/p/lapok/ATMOS/';

if isunix
	curdir = atmosdir;
	BECKdir = 'MetData/SenatorBeck/';
	LOSdir = 'MetData/PNW_SNOTEL/';
	YOSdir = 'MetData/Yosemite/';
	CSLdir = 'MetData/CSL/';
	RMEdir = 'MetData/Reynolds/WY2009/';
	MODISdir = 'RadiationData/MODIS.SW/';
	CERESdir = 'RadiationData/CERES.SYN/';
	MODdir = 'SnowHydrology/SW_Intercomparison/';
	PRINTdir = 'SnowHydrology/SW_Intercomparison/Graphics';
	if ismac
		curdir = homedir;
	end
elseif ispc
	curdir = workdir;
	BECKdir = 'MetData\SenatorBeck\';
	LOSdir = 'MetData\PNW_SNOTEL\';
	YOSdir = 'MetData\Yosemite\';
	CSLdir = 'MetData\CSL\';
	RMEdir = 'MetData\Reynolds\WY2009\';
	MODISdir = 'RadiationData\MODIS.SW\';
	CERESdir = 'RadiationData\CERES.SYN\';
	MODdir = 'SnowHydrology\SW_Intercomparison\';
	PRINTdir = 'SnowHydrology\SW_Intercomparison\Graphics';
end

% Get MET data
cd([curdir,MODdir])
load 3hrAggMET.mat

% Get UEB model output
cd UEB_Results
s = dir;

for n = 1:length(s)

	if strcmp(s(n).name(1),'.')
		continue
	end

	nind = strfind(s(n).name,'.');
	site = s(n).name(1:nind(1)-1);
	SWtype = s(n).name(nind(1)+1:nind(2)-1);

	UEBMOD.(site).(SWtype) = load(s(n).name);

end

s = fieldnames(UEBMOD);

%%%%%%%%%%%%%%%%%%
%% SNOW METRICS %%
%%%%%%%%%%%%%%%%%%
for n = 1:length(s)
	% UEB
	[UEBMOD.(s{n}).ObsSW.n_MAX,UEBMOD.(s{n}).ObsSW.MAXSWE,UEBMOD.(s{n}).ObsSW.t_MAX,...
		UEBMOD.(s{n}).ObsSW.SDD,UEBMOD.(s{n}).ObsSW.t_SDD,UEBMOD.(s{n}).ObsSW.AvgMelt,...
		UEBMOD.(s{n}).ObsSW.MELT,UEBMOD.(s{n}).ObsSW.ACC] = SnowMetrics(MET.(s{n}).t,UEBMOD.(s{n}).ObsSW.statev(:,2));
end


%% EOF analysis

%for n = 1:length(s)
% Limit focus right now to just one site
% Grab the hourly data
cd([curdir,BECKdir])
load SWA.mat
% Check out other sites and times
cd([curdir,RMEdir])
load RMEWY2009.mat



% Variables to correlate with melt:
% Other ideas: solar zenith angle, time since last snow fall

% remove mean and standard deviation from each factor considered - variance of 1  
CF = interp1(SWA.t(:,7),SWA.CF_Ratio,MET.SWA.t(:,7)); 			% Interpolate ratio derived cloud fraction to the model timestep
CF = RemoveVariance(CF); 
T = RemoveVariance(MET.SWA.T); 
RH = RemoveVariance(MET.SWA.RH); 
WIND = RemoveVariance(MET.SWA.WIND); 
SWdwn = RemoveVariance(MET.SWA.SWdwn); 
LW_Emp = RemoveVariance(MET.SWA.LW_Emp);
% use days with melt.
md = find(UEBMOD.SWA.ObsSW.MELT ~= 0);
% MELT = UEBMOD.SWA.ObsSW.MELT;
MELT = -RemoveVariance(UEBMOD.SWA.ObsSW.MELT(md));               % Divide by mean so that zero values are unaffected

t = MET.SWA.t(:,7);

%%%%%%%%%
%% MCA %%
%%%%%%%%%
% Not so useful? My covariance matrix is 7x1 - only one siginficant
% eigenvalue.
xm = [CF(md), T(md), RH(md), WIND(md), SWdwn(md), LW_Emp(md)];
[rmsq,cx,sx,cy,sy,xy,ux_cov,uy_cov,s_cov] = MaxCovAnalysis(xm,MELT);
% Heterogeneous mapping
xt = MELT*ux_cov(:,1)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EOF/PC and Rotated EOFs %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make into one big data matrix
N = length(MELT);
M = 7;
xm = [MELT,xm];

[C,ralph,E,EM,L,D,Z,ZM] = EOF_PC_ANALYSIS(xm');
ind = M+1-[1:M];

% Plot Eigenvalues
figure
errorbar(ind,D,L)
set(get(1,'CurrentAxes'),'FontSize',14)
xlabel('Index','Fontsize',14)
ylabel('Eigenvalue','Fontsize',14)

% Plot EOFs
figure
plot(E(:,M-2),'g-.')
hold on
plot(E(:,M-1),'b--')
plot(E(:,M),'r')
set(gca,'XTick',1:8,'XTickLabel',{'Melt','CF','T','RH','WIND','SW','LW','Precip'})
legend('3','2','1','Location','EastOutside')
axis tight
grid on

%% Plot PCs

% Look at the melt season (May -> June)
tind = find(datenum(2009,5,15,0,0,0) == t):find(datenum(2009,6,7,0,0,0) == t);
% Convert from melt only to all time steps
PCs = NaN(7,length(t));
PCs(:,md) = Z(:,:);

figure
subaxis(4,1,1,'P',0,'MT',.06,'ML',.12,'sv',.06,'mr',0)
hold on
plot(t(tind),PCs(M-1,tind),'b--')
plot(t(tind),PCs(M,tind),'r')
legend('2','1','Location','East')
axis tight
grid on
datetick('x','keeplimits')
ylabel('PC Amplitude')
subaxis(2)                                                  %%% EOF 2
plot(t(tind),CF(tind))
hold all
plot(t(tind),RH(tind))
plot(t(tind),LW_Emp(tind))
legend('CF','RH','LW','Location','East')
axis tight
grid on
datetick('x','keeplimits')
ylabel('Standardized')
subaxis(3)                                                 %%% EOF 1
plot(t(tind),T(tind))
hold all
plot(t(tind),SWdwn(tind))
legend('T','SW','Location','East')
axis tight
grid on
datetick('x','keeplimits')
ylabel('Standardized')
subaxis(4)                                                 %%% MELT
plot(t(tind),-UEBMOD.SWA.ObsSW.MELT(tind))
axis tight
grid on
datetick('x','keeplimits')
ylabel('Simulated Melt (mm)')

% OK, let's try rotating these eigenvectors
% Determines number of EOFs to rotate
ii=2;                                   % Only rotate the significant ones from above
lambda=E(:,M-ii+1:M);
[ER, V] = varimax(lambda,1.,1.0e-6);    % Check and see if this is an oblique or orthogonal rotation
%**************************************************

figure
plot(ER(:,end),'y','LineWidth',2)
hold on
plot(ER(:,end-1),'g-.','LineWidth',2)
set(gca,'XTick',1:8,'XTickLabel',{'Melt','CF','T','RH','WIND','SW','LW','Precip'})



