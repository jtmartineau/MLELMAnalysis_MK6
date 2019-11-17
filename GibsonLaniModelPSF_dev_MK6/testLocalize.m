% this script tests LocalizeClass
clc
clear 

GLMPP=GLMPPClass;
GLMPP=GLMPP.init('GLMPP');
GLMPP=GLMPP.update('/strgA/July2019LiveCellImaging_exported_tiffs/GLMPP_XPP_20190827002755');

PupilPlane=PupilPlaneClass;
PupilPlane=PupilPlane.init('PupilPlane');
PupilPlane=PupilPlane.update('/strgA/July2019LiveCellImaging_exported_tiffs/PupilPlane_XPP_20190827002755');

LookupTableStructure=load('/strgA/July2019LiveCellImaging_exported_tiffs/LookupTable_XPP_20190827002755.mat');

numbertrials=250;
error=cell(1,numbertrials);

iternumin=50;
Nphotonin=1000;
oiin=0.1*ones(GLMPP.K);
variin=0.1*ones(GLMPP.K);
giin=2*ones(GLMPP.K);
parfor ii=1:numbertrials
xyzlin=0.3*rand(4,1);
xyzlin(4)=0.2*rand(1)+0.5;

GibsonLaniPSF=GibsonLaniPSFClass;
GibsonLaniPSF=GibsonLaniPSF.init(GLMPP,PupilPlane);
GibsonLaniPSF=GibsonLaniPSF.SimData(Nphotonin, variin, giin, oiin,...
    xyzlin(1), xyzlin(2), xyzlin(3), xyzlin(4));

Localize=LocalizeClass;
Localize=Localize.init(GibsonLaniPSF.simDATA,...
    GibsonLaniPSF.simBACKGROUND,variin,oiin,giin,...
    LookupTableStructure,iternumin);

Localize=Localize.guess(GLMPP);
Localize=Localize.data(GibsonLaniPSF);

error{ii}=abs(Localize.ftvl1-xyzlin);
end
