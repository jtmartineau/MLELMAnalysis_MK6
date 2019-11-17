% this function tests mu_nPointSpectrum written with mu_Cools

clc
clear

SpectralStructure=load('F:\Google Drive\PSFInformationAnalysis\ReticleFunctionAnalysis\Matlab Files\BornWolfPSF\FSB_YG_OR_DR\SpectralStructure_230BFSB_515_580_680.mat');
SpectralCell=struct2cell(SpectralStructure);

analysisfile='F:\Google Drive\PSFInformationAnalysis\ReticleFunctionAnalysis\Matlab Files\BornWolfPSF\MLELMAnalysisMK1\analysisfile--2019-03-25-14-27-03';%uigetdir();
datafile='datafile--2019-03-25-14-27-03.xls';%uigetfile;
dataxls=fullfile(analysisfile,datafile);

K_X=xlsread(dataxls,'BasicParameters','C2');
K_Y=xlsread(dataxls,'BasicParameters','D2');
PIF=xlsread(dataxls,'BasicParameters','E2');
NA=xlsread(dataxls,'BasicParameters','G2');
N_photon=1000;
[num,str,raw]=xlsread(dataxls,'LookupTable','G2:G10');
lambda_vec=num;
nw=xlsread(dataxls,'BasicParameters','H2');
np=xlsread(dataxls,'BasicParameters','I2');
M=xlsread(dataxls,'BasicParameters','J2');
x=xlsread(dataxls,'x_table');

y=xlsread(dataxls,'y_table');

%rho=xlsread(dataxls,'rho_calib');
%theta=xlsread(dataxls,'theta_calib');
%THETA=xlsread(dataxls,'THETA');
Cools_struct=Cools_rhotheta_S2();
rhotheta=Cools_struct.rhotheta;
Bstar=Cools_struct.Bstar;
max_stroke=xlsread(dataxls,'PupilPlane','H2');
THETA_Cools=Create_THETA_Cools(analysisfile,datafile);
%THETA_Cools=zeros(size(THETA_Cools));

x0=0.0*M;
y0=0.0*M;
z0=0.0*(M^2);
%%
SpectralIndex=2;
Wavelength=SpectralCell{SpectralIndex}.Wavelength;
Spectrum=SpectralCell{SpectralIndex}.Spectrum;
Filter=SpectralCell{SpectralIndex}.Filter;
n=20;
delta_lambda=0.002;

tic
mu=mu_nPointSpectrum(K_X,K_Y,PIF,N_photon,NA,n,delta_lambda,Wavelength,Spectrum,Filter,x,x0,M,y,y0,rhotheta,Bstar,nw,z0,np,max_stroke,THETA_Cools);
toc
imagesc(mu),axis image, colormap(gray),colorbar;