% test GLPSF

clc
clear

GLMPP=GLMPPClass;
GLMPP=GLMPP.init('GLMPP');
%GLMPP=GLMPP.update('/media/strgA/July2019LiveCellImaging_exported_tiffs/CAL_20190911040227/GLMPP_XPP_20190911040227.csv');
PupilPlane=PupilPlaneClass;
PupilPlane=PupilPlane.init('PupilPlane');
%PupilPlane=PupilPlane.update('/media/strgA/July2019LiveCellImaging_exported_tiffs/PupilPlane_Optimalfor_SEMPNTB_AXF0647_AXF0660_AXF0680_AXF0700_imgaussblur1.csv');
SpectralStructure=load('/media/strgA/July2019LiveCellImaging_exported_tiffs/SpectralStructure_SEMPNTB_AXF0647_AXF0660_AXF0680_AXF0700.mat');

fnames=fieldnames(SpectralStructure);
Wavelength=zeros(numel(fnames),numel(SpectralStructure.(fnames{1}).Wavelength));
Spectrum=zeros(size(Wavelength));
Filter=zeros(size(Spectrum));
for ii=1:numel(fnames)
    Wavelength(ii,:)=SpectralStructure.(fnames{ii}).Wavelength;
    Spectrum(ii,:)=SpectralStructure.(fnames{ii}).Spectrum;
    Filter(ii,:)=SpectralStructure.(fnames{ii}).Filter;
end

K=GLMPP.K;
PixelSize=GLMPP.PixelSize;
M=GLMPP.M;

oiin=0.1*ones(K);
variin=0.0*ones(K);
giin=1*ones(K);

x0in=0; 
y0in=0;
zpin=0;
lambdain=0.760;%SCOM(filtersin.*spectrain,wavelengthsin);
Nphoton=3000;

GLMPP.DELTAt=0;

GibsonLaniPSF=GibsonLaniPSFClass;
GibsonLaniPSF=GibsonLaniPSF.init(GLMPP,PupilPlane);
%tic
%GibsonLaniPSF=GibsonLaniPSF.intens(x0in, y0in, zpin, lambdain);
%toc
%figure(1)
%imagesc((Nphoton*GibsonLaniPSF.PSF)), axis image, colormap gray;
%tic
%GibsonLaniPSF=GibsonLaniPSF.bead(x0in,y0in,zpin,lambdain);
%toc
%%
specidx=4;
tic
GibsonLaniPSF=GibsonLaniPSF.spectralintens(x0in, y0in, zpin,...
    Wavelength(specidx,:), Spectrum(specidx,:), Filter(specidx,:));
toc
figure(2)
imagesc((1000*GibsonLaniPSF.spectralPSF)), axis image, colormap gray;
%%
residual=abs(GibsonLaniPSF.spectralPSF-GibsonLaniPSF.PSF);
figure(3)
subplot(1,3,1)
imagesc((Nphoton*GibsonLaniPSF.PSF)), axis image, colormap gray; colorbar;
subplot(1,3,2)
imagesc(Nphoton*residual), axis image, colormap gray; colorbar;
subplot(1,3,3)
imagesc((Nphoton*GibsonLaniPSF.spectralPSF)), axis image, colormap gray; colorbar;
%%
GibsonLaniPSF=GibsonLaniPSF.spectralintensstandard(x0in, y0in, zpin, wavelengthsin, spectrain, filtersin); 
figure(4)
imagesc(Nphoton*GibsonLaniPSF.PSF), axis image, colormap gray; colorbar
%%
tic
GibsonLaniPSF=GibsonLaniPSF.ddx0(x0in, y0in, zpin, lambdain);
toc
tic
GibsonLaniPSF=GibsonLaniPSF.ddy0(x0in, y0in, zpin, lambdain);
toc
tic
GibsonLaniPSF=GibsonLaniPSF.ddzp(x0in, y0in, zpin, lambdain);
toc
tic
GibsonLaniPSF=GibsonLaniPSF.ddlambda(x0in, y0in, zpin, lambdain);
toc
tic
GibsonLaniPSF=GibsonLaniPSF.d2dx02(x0in, y0in, zpin, lambdain);
toc
tic
GibsonLaniPSF=GibsonLaniPSF.d2dy02(x0in, y0in, zpin, lambdain);
toc
tic
GibsonLaniPSF=GibsonLaniPSF.d2dzp2(x0in, y0in, zpin, lambdain);
toc
tic
GibsonLaniPSF=GibsonLaniPSF.d2dlambda2(x0in, y0in, zpin, lambdain);
toc
%figure('Name','intensity'); imagesc(GibsonLaniPSF.PSF), axis image, colormap gray, title("intensity");
%figure('Name','ddx0'); imagesc(GibsonLaniPSF.ddx0PSF), axis image, colormap gray, title("ddx0");
%figure('Name','ddy0'); imagesc(GibsonLaniPSF.ddy0PSF), axis image, colormap gray, title("ddy0");
%figure('Name','ddzp'); imagesc(GibsonLaniPSF.ddzpPSF), axis image, colormap gray, title("ddzp");
%figure('Name','ddlambda'); imagesc(GibsonLaniPSF.ddlambdaPSF), axis image, colormap gray, title("ddlambda");
%figure('Name','d2dx02'); imagesc(GibsonLaniPSF.d2dx02PSF), axis image, colormap gray, title("d2dx02");
%figure('Name','d2dy02'); imagesc(GibsonLaniPSF.d2dy02PSF), axis image, colormap gray, title("d2dy02");
%figure('Name','d2dzp2'); imagesc(GibsonLaniPSF.d2dzp2PSF), axis image, colormap gray, title("d2dzp2");
%figure('Name','d2dlambda2'); imagesc(GibsonLaniPSF.d2dlambda2PSF), axis image, colormap gray, title("d2dlambda2");
tic
GibsonLaniPSF=GibsonLaniPSF.SimData(Nphoton, variin, giin, oiin, x0in, y0in, zpin, lambdain);
toc
figure('Name','intensity'); imagesc(GibsonLaniPSF.simDATA), axis image, colormap gray, title("intensity"),colorbar;
tic
GibsonLaniPSF=GibsonLaniPSF.loglikelihood(GibsonLaniPSF.simDATA, GibsonLaniPSF.simBACKGROUND, x0in, y0in, zpin, lambdain);
toc
tic
GibsonLaniPSF=GibsonLaniPSF.llsCMOS(GibsonLaniPSF.simDATA,...
                GibsonLaniPSF.simBACKGROUND, ...
                variin, giin, oiin, x0in, y0in, zpin, lambdain);
[xcom,ycom]=COM((GibsonLaniPSF.simDATA-GibsonLaniPSF.oi)./GibsonLaniPSF.gi-GibsonLaniPSF.simBACKGROUND,...
    GLMPP.PixelSize,GLMPP.M);
xyzlin=[xcom,ycom,0.650,0.630];
iternumin=150;
%tic
%GibsonLaniPSF = GibsonLaniPSF.NRfitxyzlllsCMOS(GibsonLaniPSF.simDATA,...
%    GibsonLaniPSF.simBACKGROUND, variin, giin, oiin, xyzlin, iternumin);
%toc
tic
GibsonLaniPSF = GibsonLaniPSF.fminsearchxyzllsCMOS(GibsonLaniPSF.simDATA,...
                GibsonLaniPSF.simBACKGROUND, variin, giin, oiin, xyzlin, iternumin);
toc           
display(GibsonLaniPSF.xyzl);
GibsonLaniPSF=GibsonLaniPSF.intens(GibsonLaniPSF.xyzl(1),...
     GibsonLaniPSF.xyzl(2), GibsonLaniPSF.xyzl(3), GibsonLaniPSF.xyzl(4));
figure('Name','intensity'); imagesc(GibsonLaniPSF.Nphoton.*GibsonLaniPSF.PSF.*...
    GibsonLaniPSF.gi+GibsonLaniPSF.simBACKGROUND+GibsonLaniPSF.oi), axis image, colormap gray, title("intensity"), colorbar;