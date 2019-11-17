%this script tests AnalyzeClass
clc 
clear 

EXPdirin='/strgA/July2019LiveCellImaging_exported_tiffs';
DATAdirin='/20190730_XPSF_LiveCell_SetAH_tiff';
PSFDATAdirin='/20190726_XPSF_Calibration_tiff';
DATASpectralStructureFileNamein='/SpectralStructure_BRUVUT_SJF0549_AXF0647';
PSFDATASpectralStructureFileNamein='/SpectralStructure_BRUVUT_FSB0515_FSB0580_FSB0680';
PhasePlateTypein='XPP';
timestampin='20190827002755';
startingframein=1;
numberframesin=0;
iternumin=100;
plotyesnoin='yes';
excitlambdain=[0.640,0.561];
numberposin=41;
framesatposin=10;
Analyze=AnalyzeClass;
Analyze=Analyze.initanalysis(EXPdirin,timestampin);
tic
Analyze=Analyze.frames(startingframein,numberframesin,iternumin,...
    plotyesnoin);
toc
Analyze=Analyze.prepareparticlescsv;
Analyze=Analyze.preparedataxml;

