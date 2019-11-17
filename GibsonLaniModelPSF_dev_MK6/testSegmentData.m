%this script tests SegmentDataClass
clc
clear

GLMPPFileName='/strgA/July2019LiveCellImaging_exported_tiffs/GLMPP_XPP_20190827002755';
LookupTableFileName='/strgA/July2019LiveCellImaging_exported_tiffs/LookupTable_XPP_20190827002755.csv';
PupilPlaneFileName='/strgA/July2019LiveCellImaging_exported_tiffs/PupilPlane_XPP_20190827002755.csv';
SpectralStructureFileName='/strgA/July2019LiveCellImaging_exported_tiffs/SpectralStructure_BRUVUT_SJF0549_AXF0647.mat';
oiFileName='/strgA/PIXSTAT/offset_map.tiff';
variFileName='/strgA/PIXSTAT/variance_map.tiff';
giFileName='/strgA/PIXSTAT/gain_map.tiff';
DATAdir='/strgA/July2019LiveCellImaging_exported_tiffs/20190726_XPSF_LiveCell_SetA_tiff';

maindirectory=pwd;
cd(DATAdir);
filelisting=dir('*.tif');
cd(maindirectory);

framenumber=10;
framein=cast(imread(fullfile(DATAdir,filelisting(framenumber).name)),'single');
oiin=cast(imread(oiFileName),'single');
variin=cast(imread(variFileName),'single');
giin=cast(imread(giFileName),'single');
imagesc(framein), axis image, colormap(gray); colorbar

GLMPP=GLMPPClass;
GLMPP=GLMPP.init('GLMPP');
GLMPP=GLMPP.update(GLMPPFileName);

SegmentData=SegmentDataClass;
SegmentData=SegmentData.init(GLMPP,framein,oiin,variin,giin);
tic
SegmentData=SegmentData.selectROIs(GLMPP);
toc
SegmentData=SegmentData.drawROIs;
SegmentData=SegmentData.plotROIs;

