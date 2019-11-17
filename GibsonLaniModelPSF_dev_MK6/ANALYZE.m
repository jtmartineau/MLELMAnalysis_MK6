clc
clear

EXPdirin='/media/strgA/July2019LiveCellImaging_exported_tiffs';
caltimestampin='20190911040227';
DATASpectralStructureFileNamein='/SpectralStructure_BRUVUT_SJF0549_AXF0647.mat';
DATAdirin='/20190730_XPSF_LiveCell_SetAA_tiff';

Analyze=AnalyzeClass;
Analyze=Analyze.initanalysis(EXPdirin,caltimestampin,...
                DATASpectralStructureFileNamein,DATAdirin);
 
startingframein=1;
numberframesin=9999;
iternumin=100;
plotyesnoin='no';

tic            
Analyze=Analyze.frames(startingframein,numberframesin,...
                iternumin,plotyesnoin);
toc
Analyze=Analyze.prepareparticlescsv;
Analyze=Analyze.preparedataxml;

cd(Analyze.outdir);