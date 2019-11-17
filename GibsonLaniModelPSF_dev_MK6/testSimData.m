%test SimData
clc
clear

PIXSTATFileNamein = '/media/strgA/July2019LiveCellImaging_exported_tiffs/PIXSTAT';
GLMPPFileNamein = '/media/strgA/July2019LiveCellImaging_exported_tiffs/CAL_20190911040227/GLMPP_XPP_20190911040227.csv';
PupilPlaneFileNamein = '/media/strgA/July2019LiveCellImaging_exported_tiffs/CAL_20190911040227/PupilPlane_XPP_20190911040227.csv';
bgmeanin = 25;
bgsigmain = 5;
Nphotonvecin = [100,200,300,500,600,800,1000,2000,3000,5000,6000,8000,10000,20000,30000,50000,60000,80000,100000];
SpectralStructureFileNamein = '/media/strgA/July2019LiveCellImaging_exported_tiffs/SpectralStructure_SEMPNTB_AXF0647_AXF0660_AXF0680_AXF0700.mat';
OptimalPupilPlaneFileNamein = '/media/strgA/July2019LiveCellImaging_exported_tiffs/PupilPlane_Optimalfor_SEMPNTB_AXF0647_AXF0660_AXF0680_AXF0700_imgaussblur1_1.csv';
SS = load(SpectralStructureFileNamein);
fnames=fieldnames(SS);
fluortagsin=fnames';
filterwidthvecinin = ones(1,10);
mincutoffin = 0.650;
maxcutoffin = 0.850;
deltax0 = single(0.05);
deltay0 = single(0.05);
deltazp = single(0.05);
x0range = single([-0.05,0.05]);
y0range = single([-0.05,0.05]);
zprange = single([-0.05,0.05]);
xoffsetin = 0;
yoffsetin = 0;
thetaoffsetin = 0;
liveactuatorsin = [1,1,1,1];
actuatorheightsin = [0,2.4430,2.4355,2.1515];

SimData=SimDataClass;
SimData=SimData.init(PIXSTATFileNamein,GLMPPFileNamein,...
                PupilPlaneFileNamein,bgmeanin,bgsigmain,Nphotonvecin,fluortagsin,...
                SpectralStructureFileNamein,...
                filterwidthvecinin,...
                mincutoffin,...
                maxcutoffin,...
                deltax0,...
                deltay0,...
                deltazp,...
                x0range,...
                y0range,...
                zprange,...
                xoffsetin,...
                yoffsetin,...
                thetaoffsetin,...
                liveactuatorsin,...
                actuatorheightsin);

%%
x0in = 0;
y0in = 0;
zpin = 0;
Nphotonin = 1000;
lbin = [0,0,0,0]; %plate thickness in microns
ubin = [0,5,5,5];


GLMPP=GLMPPClass;
GLMPP=GLMPP.init('GLMPP_temp.csv');
%%GLMPP=GLMPP.update(GLMPPFileNamein);

PupilPlane=PupilPlaneClass;
PupilPlane=PupilPlane.init('PupilPlane_temp.csv');
%PupilPlane=PupilPlane.update(OptimalPupilPlaneFileNamein);

GibsonLaniPSF=GibsonLaniPSFClass;
GibsonLaniPSF=GibsonLaniPSF.init(...
    GLMPP,PupilPlane);

SimData.SpectralStructure=...
                load(SimData.SpectralStructureFileName);
          
            
%out = SimData.PIXCON_OBJ(PupilPlane,...
%                GibsonLaniPSF, x0in, y0in, zpin, Nphotonin,...
%                [0,.29,.29,0]);
%
fun=@(actuatorheightsin)SimData.PIXCON_OBJ(PupilPlane,...
                GibsonLaniPSF, x0in, y0in, zpin, Nphotonin,...
                actuatorheightsin);

options = optimoptions(@particleswarm,...
    'PlotFcn', 'pswplotbestf', 'HybridFcn', @fmincon,...
    'UseParallel', true, 'SwarmSize', 100);

SimData.optimalactuatorheights=...
    particleswarm(fun,numel(lbin),lbin,ubin,options); 
        
PupilPlane.THETAplate=XPP(PupilPlane.rho,...
                PupilPlane.theta,...
                SimData.xoffset,...
                SimData.yoffset,...
                SimData.thetaoffset,...
                SimData.liveactuators,...
                SimData.optimalactuatorheights);

PupilPlane.write(OptimalPupilPlaneFileNamein);
%%            
%SimData=SimData.PIXCON_MIN(OptimalPupilPlaneFileNamein,...
%                x0in, y0in, zpin, Nphotonin, lbin, ubin, 'ga');   
 
%%    

PIXSTATFileNamein = SimData.PIXSTATFileName;
GLMPPFileNamein = SimData.GLMPPFileName;
PupilPlaneFileNamein = OptimalPupilPlaneFileNamein;
bgmeanin = SimData.bgmean;
bgsigmain = SimData.bgsigma;
Nphotonvecin = SimData.Nphotonvec;
SpectralStructureFileNamein = SimData.SpectralStructureFileName;
SS = load(SpectralStructureFileNamein);
fnames=fieldnames(SS);
fluortagsin=fnames';
LTFilterSetFileNamein = '/media/strgA/July2019LiveCellImaging_exported_tiffs/LookupTable_FitlerSet.mat';
LTPhasePlateFileNamein = '/media/strgA/July2019LiveCellImaging_exported_tiffs/LookupTable_PhasePlate.mat';
deltax0 = single(0.05);
deltay0 = single(0.05);
deltazp = single(0.05);
x0range = single([-0.05,0.05]);
y0range = single([-0.05,0.05]);
zprange = single([-0.05,0.05]);
xoffsetin = SimData.xoffset;
yoffsetin = SimData.yoffset;
thetaoffsetin = SimData.thetaoffset;
liveactuatorsin = SimData.liveactuators;
actuatorheightsin = SimData.actuatorheights;

SimData=SimDataClass;
SimData=SimData.init(PIXSTATFileNamein,GLMPPFileNamein,...
                PupilPlaneFileNamein,bgmeanin,bgsigmain,Nphotonvecin,fluortagsin,...
                SpectralStructureFileNamein,...
                filterwidthvecinin,...
                mincutoffin,...
                maxcutoffin,...
                deltax0,...
                deltay0,...
                deltazp,...
                x0range,...
                y0range,...
                zprange,...
                xoffsetin,...
                yoffsetin,...
                thetaoffsetin,...
                liveactuatorsin,...;
                actuatorheightsin);
%
SimData=SimData.createfiltersLookupTable(LTFilterSetFileNamein);
%
LT=load(SimData.LTFilterSetFileName);

for ii=1:size(LT.LookupTable,3)
    imagesc(LT.LookupTable(:,:,ii)), axis image, colormap gray
    pause(0.03)
end
%
SimData=SimData.createphaseplateLookupTable(LTPhasePlateFileNamein);
%
LT=load(SimData.LTPhasePlateFileName);

for ii=1:size(LT.LookupTable,3)
    imagesc(poissrnd(1000*LT.LookupTable(:,:,ii))), axis image, colormap gray
    pause(0.1)
end
%test SimData

PIXSTATFileNamein = '/media/strgA/July2019LiveCellImaging_exported_tiffs/PIXSTAT';
GLMPPFileNamein = '/media/strgA/July2019LiveCellImaging_exported_tiffs/CAL_20190911040227/GLMPP_XPP_20190911040227.csv';
PupilPlaneFileNamein = OptimalPupilPlaneFileNamein;
bgmeanin = 25;
bgsigmain = 5;
Nphotonvecin = [100,200,300,500,600,800,1000,2000,3000,5000,6000,8000,10000,20000,30000,50000,60000,80000,100000];
numbertrialsin=100000;
SpectralStructureFileNamein = '/media/strgA/July2019LiveCellImaging_exported_tiffs/SpectralStructure_SEMPNTB_AXF0647_AXF0660_AXF0680_AXF0700.mat';
SS = load(SpectralStructureFileNamein);
fnames=fieldnames(SS);
fluortagsin=fnames';
LTFilterSetFileNamein =  '/media/strgA/July2019LiveCellImaging_exported_tiffs/LookupTable_FitlerSet.mat';
LTPhasePlateFileNamein = '/media/strgA/July2019LiveCellImaging_exported_tiffs/LookupTable_PhasePlate.mat';
 
SimData=SimDataClass;
SimData=SimData.init(PIXSTATFileNamein,GLMPPFileNamein,...
                PupilPlaneFileNamein,bgmeanin,bgsigmain,Nphotonvecin,fluortagsin,...
                LTFilterSetFileNamein,...
                LTPhasePlateFileNamein);


LT=load(SimData.LTFilterSetFileName);

for ii=1:size(LT.LookupTable,3)
    imagesc(LT.LookupTable(:,:,ii)), axis image, colormap gray
    pause(0.1)
end

LT=load(SimData.LTPhasePlateFileName);

for ii=1:size(LT.LookupTable,3)
    imagesc(LT.LookupTable(:,:,ii)), axis image, colormap gray
    pause(0.1)
end

SimData=SimData.spectralMonteCarlo(numbertrialsin,'FilterSet');
SimData=SimData.maketruthtablearray;

SimData=SimData.spectralMonteCarlo(numbertrialsin,'PhasePlate');
SimData=SimData.maketruthtablearray; 

%save the data before plotting
mkdir(strcat(datestr(now,'YYYYMMDDmmSS'),'MonteCarlo'))
save(strcat(datestr(now,'YYYYMMDDmmSS'),'MonteCarlo','/results','.m'))

figure('Name','MonteCarloResults')
for ii=1:numel(SimData.SpecCategories)
    semilogx(...
            (SimData.Nphotonvec),...
            100*(reshape(SimData.filtersettruthtablearray(ii,ii,:),[1,numel(SimData.Nphotonvec)]))/numbertrialsin,'--');
    %legend(strcat('Filter Set','_',SimData.fluortags{ii}))
    hold on
end
for ii=1:numel(SimData.SpecCategories)
    semilogx(...
            (SimData.Nphotonvec),...
            100*(reshape(SimData.phaseplatetruthtablearray(ii,ii,:),[1,numel(SimData.Nphotonvec)]))/numbertrialsin,'-');
            %legend(strcat('Phase Plate',' ',SimData.fluortags{ii}))
end
legendarray = cell(1,numel(SimData.fluortags));
for ii = 1:numel(legendarray)
    legendarray{ii} = SimData.fluortags{ii};
end
legendarray=cat(2,legendarray,legendarray);
legend(legendarray)
hold off;
ylabel('% correct per spectrum');

