%this script tests AnalyzeClass
clc 
clear 

EXPdirin='/media/strgA/July2019LiveCellImaging_exported_tiffs';
PSFDATAdirin='/20190726_XPSF_Calibration_tiff';
DATASpectralStructureFileNamein='/SpectralStructure_BRUVUT_SJF0549_AXF0647';
PSFDATASpectralStructureFileNamein='/SpectralStructure_BRUVUT_FSB0515_FSB0580_FSB0680';
PhasePlateTypein='XPP';

Analyze=AnalyzeClass;
Analyze=Analyze.initcalibrate(EXPdirin,...
            PSFDATAdirin,DATASpectralStructureFileNamein,...
            PSFDATASpectralStructureFileNamein,PhasePlateTypein);
        
beadtagsin={'FSB0680','FSB0580'};
FluorTagsin={'AXF0647','SJF0549'};
excitlambdain=[0.640,0.561];
numberposin=41;
framesatposin=10;

Analyze=Analyze.calibrate(beadtagsin,FluorTagsin,excitlambdain,...
                            numberposin,framesatposin);
                     
cd(Analyze.CALdir);                        
