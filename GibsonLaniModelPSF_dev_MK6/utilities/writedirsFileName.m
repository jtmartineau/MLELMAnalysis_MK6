% this script writes the needed filenames to a csv file named FileNamedir
clc
clear

timestamp='20190827002755';

EXPdir='/strgA/July2019LiveCellImaging_exported_tiffs';
GLMPPFileName='/strgA/July2019LiveCellImaging_exported_tiffs/GLMPP_XPP_20190827002755';
LookupTableFileName='/strgA/July2019LiveCellImaging_exported_tiffs/LookupTable_XPP_20190827002755.mat';
PupilPlaneFileName='/strgA/July2019LiveCellImaging_exported_tiffs/PupilPlane_XPP_20190827002755';
DATASpectralStructureFileName='/strgA/July2019LiveCellImaging_exported_tiffs/SpectralStructure_BRUVUT_SJF0549_AXF0647.mat';
PSFDATASpectralStructureFileName='/strgA/July2019LiveCellImaging_exported_tiffs/SpectralStructure_BRUVUT_FSB0515_FSB0580_FSB0680.mat';
oiFileName='/strgA/PIXSTAT/offset_map.tiff';
variFileName='/strgA/PIXSTAT/variance_map.tiff';
giFileName='/strgA/PIXSTAT/gain_map.tiff';
DATAdir='/strgA/July2019LiveCellImaging_exported_tiffs/20190730_XPSF_LiveCell_SetAH_tiff';
PSFDATAdir='/strgA/July2019LiveCellImaging_exported_tiffs/20190726_XPSF_Calibration_tiff';

FileIndex={'GLMPP',GLMPPFileName;...
    'LookupTable',LookupTableFileName;...
    'PupilPlane',PupilPlaneFileName;...
    'DATASpectralStructure',DATASpectralStructureFileName;...
    'PSFDATASpectralStructure',PSFDATASpectralStructureFileName;...
    'oi',oiFileName;...
    'vari',variFileName;...
    'gi',giFileName;...
    'DATA',DATAdir;...
    'PSFDATA',PSFDATAdir};

outdir=strcat(EXPdir,'/','OUT','_',timestamp);
mkdir(outdir);
FileIndexFileName=strcat(outdir,'/FileIndex','.csv');
writecell(FileIndex,FileIndexFileName,'Delimiter','comma');