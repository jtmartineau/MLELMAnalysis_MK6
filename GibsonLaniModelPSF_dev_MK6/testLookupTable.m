%this script tests LookupTableClass
clc
clear

SpectralStructureFileName=...
    '/media/strgA/July2019LiveCellImaging_exported_tiffs/SpectralStructure_BRUVUT_SJF0549_AXF0647.mat';
FluorTagsin={'SJF0549','AXF0647'};

GLMPP=GLMPPClass;
GLMPP=GLMPP.init('GLMPP.csv');
GLMPP=GLMPP.update(...
    '/media/strgA/July2019LiveCellImaging_exported_tiffs/CAL_20191117210929/GLMPP_XPP_20191117210929.csv');
GLMPP=GLMPP.write('GLMPP.csv');

PupilPlane=PupilPlaneClass;
PupilPlane=PupilPlane.init('PupilPlane.csv');
PupilPlane=PupilPlane.update(...
    '/media/strgA/July2019LiveCellImaging_exported_tiffs/CAL_20191117210929/PupilPlane_XPP_20191117210929.csv');
PupilPlane=PupilPlane.write('PupilPlane.csv');

%%

filterwidthvecin=ones(1,6);
mincutoffin=0.400;
maxcutoffin=0.850;
Nphoton=1000;

GibsonLaniPSF=GibsonLaniPSFClass;
GibsonLaniPSF=GibsonLaniPSF.init(GLMPP,PupilPlane);

LookupTable=LookupTableClass;
LookupTable=LookupTable.init(SpectralStructureFileName,...
    FluorTagsin,GLMPP);

%%%%%% this block creates the lookup table first for the phase plate signal
%%%%%% and then a filterset signal second. these lookup tables will
%%%%%% be used to compare the performace of a phase plate and a filter set
LookupTable=LookupTable.spectralcreate(GibsonLaniPSF);   
LookupTableFileNamePhasePlate=LookupTable.write('/media/strgA/July2019LiveCellImaging_exported_tiffs/LookupTable_PhasePlate.mat');

LookupTable=LookupTable.filtersetstandardcreate(GLMPP,...
            GibsonLaniPSF,filterwidthvecin, mincutoffin, maxcutoffin);
LookupTableFileNameFilterSet=LookupTable.write('/media/strgA/July2019LiveCellImaging_exported_tiffs/LookupTable_FitlerSet.mat');
%%%%%%
%%%%%%
%%%%%%
%%
LT=load(LookupTableFileNameFilterSet);

for ii=1:size(LT.LookupTable,3)
    imagesc(LT.LookupTable(:,:,ii)), axis image, colormap winter
    pause(0.01)
end
%%
LT=load('/media/strgA/July2019LiveCellImaging_exported_tiffs/CAL_20191117210929/LookupTable_XPP_20191117210929');

for ii=1:size(LT.LookupTable,3)
    imagesc(LT.LookupTable(:,:,ii)), axis image, colormap gray
    pause(0.05)
end
%%

GibsonLaniPSF=GibsonLaniPSFClass;
GibsonLaniPSF=GibsonLaniPSF.init(GLMPP,PupilPlane);

LookupTable=LookupTableClass;
LookupTable=LookupTable.init(SpectralStructureFileName,...
    FluorTagsin,GLMPP);
LookupTable=LookupTable.create(GibsonLaniPSF);
%LookupTable=LookupTable.createvarblur(GibsonLaniPSF);
LookupTableFileName=LookupTable.write('/strgA/July2019LiveCellImaging_exported_tiffs/CAL_20190911040227/LookupTable_XPP_20190911040227.mat');

LT=load(LookupTableFileName);

for ii=1:size(LT.LookupTable,3)
    imagesc(LT.LookupTable(:,:,ii)), axis image, colormap winter
    pause(0.010)
end

%%
  
for ii=25:size(LookupTable,3)
    imagesc(LookupTable(:,:,ii)), axis image, colormap winter
    pause(0.010)
end