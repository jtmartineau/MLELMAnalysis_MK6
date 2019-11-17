%this script calibrates a PSF.
clc 
clear

EXPdir='/strgA/July2019LiveCellImaging_exported_tiffs';
DATAdir='/strgA/July2019LiveCellImaging_exported_tiffs/20190726_XPSF_Calibration_tiff';
PIXSTATdir='/strgA/PIXSTAT';
PhasePlateType='XPP';
timestamp=string(datetime('now','TimeZone','UTC','Format','yMMddHHmmss'));
SpectralStructure=load(strcat(EXPdir,'/SpectralStructure_BRUVUT_FSB0515_FSB0580_FSB0680.mat'));
beadtagsin={'FSB0680','FSB0580'};
excitlambdain=[0.640,0.561];
numberposin=41;
framesatposin=10;

%initialize the basic microscope parameters class as GLMPP
%(GibsonLaniModelPSFParameters).
GLMPPFileName=strcat('GLMPP','_',PhasePlateType,'_',timestamp,'.csv');
GLMPP=GLMPPClass;
GLMPP=GLMPP.init(GLMPPFileName);

%initialize the class that handles the phase structure of the pupil plane
%as Pupil Plane.
PupilPlaneCalibrationFileName=strcat('PupilPlane','_',PhasePlateType,'_',timestamp,'.csv');
PupilPlane=PupilPlaneClass;
PupilPlane=PupilPlane.init(PupilPlaneCalibrationFileName);

%initialize the class that handles simulating the PSF as GibsonLaniPSF.
GibsonLaniPSF=GibsonLaniPSFClass;
GibsonLaniPSF=GibsonLaniPSF.init(GLMPP,PupilPlane);

%initialize the class that handles the PSF calibration as CalibratePSF
CalibratePSF=CalibratePSFClass;
CalibratePSF=CalibratePSF.init(GLMPP,PupilPlane,...
                SpectralStructure,beadtagsin,excitlambdain,...
                numberposin,framesatposin);
            
%import PSF calibration data.
CalibratePSF=CalibratePSF.importdata(DATAdir,PIXSTATdir);

%select ROIs.
CalibratePSF=CalibratePSF.selectROI(GLMPP);

%fit the model to the calibration data.
CalibratePSF=CalibratePSF.fminsearchllsCOMS(GibsonLaniPSF,...
    GLMPP,PupilPlane);

%transfer the pupil plane values from the fit, from CalibratePSF to
%PupilPlane.
PupilPlane.THETAplate(:)=reshape(CalibratePSF.THETAplate(1,1,1:end),...
    1,numel(CalibratePSF.THETAplate(1,1,1:end)));
PupilPlane.THETAZernike(:)=reshape(CalibratePSF.THETAZernike(1,1,1:end),...
    1,numel(CalibratePSF.THETAZernike(1,1,1:end)));

%write the propery values of GLMPP and PupilPlane to a calibration file.
GLMPP=GLMPP.write(GLMPPFileName);
PupilPlane=PupilPlane.write(PupilPlaneCalibrationFileName);

CALdir=strcat(EXPdir,'/CAL','_',timestamp);
copyfile(strcat(GLMPPFileName,'.csv'),strcat(CALdir,GLMPPFileName,'.csv'));
copyfile(strcat(PupilPlaneCalibrationFileName,'.csv'),strcat(CALdir,PupilPlaneCalibrationFileName,'.csv'));

%%
%display calibrated PSF.
CalibratePSF=CalibratePSF.genmodelarray(GibsonLaniPSF,GLMPP);
figure('Name','data and fit');
for ii=1:(size(CalibratePSF.MODELarray,3))
    subplot(1,2,1)
    imagesc(CalibratePSF.SIGNALarray(:,:,ii)), axis image, colormap winter;
    colorbar;
    subplot(1,2,2)
    imagesc(CalibratePSF.MODELarray(:,:,ii)), axis image, colormap winter;
    colorbar;
    pause(0.1);
end