% this script tests CalibratePSFClass
clc 
clear

SpectralStructure=load('SpectralStructure_BRUFSB_515_580_680.mat');
beadtagsin={'FSB0680','FSB0580'};
excitlambdain=[0.640,0.561];
numberposin=41;
framesatposin=10;

GLMPP=GLMPPClass;
GLMPP=GLMPP.init('GLMPP');

PupilPlane=PupilPlaneClass;
PupilPlane=PupilPlane.init('PupilPlane');

GibsonLaniPSF=GibsonLaniPSFClass;
GibsonLaniPSF=GibsonLaniPSF.init(GLMPP,PupilPlane);

CalibratePSF=CalibratePSFClass;

CalibratePSF=CalibratePSF.init(GLMPP,PupilPlane,...
                SpectralStructure,beadtagsin,excitlambdain,...
                numberposin,framesatposin);  

%hold on
%for ii=1:size(CalibratePSF.wavelengths,1)
%    plot(CalibratePSF.wavelengths(ii,:),CalibratePSF.spectra(ii,:));
%end
%for ii=1:size(CalibratePSF.wavelengths,1)
%    plot(CalibratePSF.wavelengths(ii,:),CalibratePSF.filters(ii,:),'k');
%end
%hold off

DATAdir='/strgA/July2019LiveCellImaging_exported_tiffs/20190726_XPSF_Calibration_tiff';
PIXSTATdir='/strgA/PIXSTAT';
CalibratePSF=CalibratePSF.importdata(DATAdir,PIXSTATdir);
%figure('Name','g_i')
%imagesc(CalibratePSF.gi), axis image, colormap gray;

%figure('Name','var_i')
%imagesc(CalibratePSF.vari), axis image, colormap gray;

%figure('Name','o_i')
%imagesc(CalibratePSF.oi), axis image, colormap gray;

%figure('Name','imported DATA');
%for ii=1:size(CalibratePSF.DATAarray,3)
%    imagesc(CalibratePSF.DATAarray(:,:,ii)), axis image, colormap winter;
%    pause(0.3);
%end

CalibratePSF=CalibratePSF.selectROI(GLMPP);
%figure('Name','selected ROIs');
%for ii=1:size(CalibratePSF.ROIarray,3)
%    imagesc(CalibratePSF.ROIarray(:,:,ii)), axis image, colormap winter;
%    colorbar;
%    pause(0.1);
%end

%figure('Name','gi');
%for ii=1:size(CalibratePSF.ROIarray,3)
%    imagesc(CalibratePSF.giarray(:,:,ii)), axis image, colormap winter;
%    colorbar;
%    pause(0.1);
%end

%figure('Name','vari');
%for ii=1:size(CalibratePSF.ROIarray,3)
%    imagesc(CalibratePSF.variarray(:,:,ii)), axis image, colormap winter;
%    colorbar;
%    pause(0.1);
%end

%figure('Name','oi');
%for ii=1:size(CalibratePSF.ROIarray,3)
%    imagesc(CalibratePSF.oiarray(:,:,ii)), axis image, colormap winter;
%    colorbar;
%    pause(0.1);
%end

%figure('Name','background');
%for ii=1:size(CalibratePSF.ROIarray,3)
%    imagesc(CalibratePSF.BACKGROUNDarray(:,:,ii)), axis image, colormap winter;
%    colorbar;
%    pause(0.1);
%end

Zweightsin=CalibratePSF.Zweights(:,1);
xoffsetin=0;
yoffsetin=0;
thetaoffsetin=0;
CalibratePSF=CalibratePSF.setpupil(PupilPlane,Zweightsin, xoffsetin, yoffsetin,...
                thetaoffsetin);
            
%x0in=0;
%y0in=0;
%zpin=0;
%lambdain=0.575;
%DELTAtin=0.5;
%CalibratePSF=CalibratePSF.genmodel(GibsonLaniPSF,...
%                x0in, y0in, zpin, lambdain, DELTAtin);
%figure('Name','model')
%imagesc(CalibratePSF.MODEL), axis image, colormap gray;

CalibratePSF=CalibratePSF.genmodelarray(GibsonLaniPSF,GLMPP);
%figure('Name','model array');
%for ii=1:size(CalibratePSF.MODELarray,3)
%    subplot(1,3,1)
%    imagesc(CalibratePSF.SIGNALarray(:,:,ii)), axis image, colormap winter;
%    colorbar;
%    subplot(1,3,2)
%    imagesc(CalibratePSF.MODELarray(:,:,ii)), axis image, colormap winter;
%    colorbar;
%        subplot(1,3,1)
%    imagesc(CalibratePSF.SIGNALarray(:,:,ii)), axis image, colormap winter;
%    colorbar;
%    subplot(1,3,3)
%    imagesc(CalibratePSF.BACKGROUNDarray(:,:,ii)), axis image, colormap winter;
%    colorbar;
%    pause(0.1);
%end

%CalibratePSF=CalibratePSF.gallsCOMS(GibsonLaniPSF,...
%    GLMPP,PupilPlane);
CalibratePSF=CalibratePSF.fminsearchllsCOMS(GibsonLaniPSF,...
    GLMPP,PupilPlane);
%CalibratePSF=CalibratePSF.fminuncllsCOMS(GibsonLaniPSF,...
%    GLMPP,PupilPlane);
%CalibratePSF=CalibratePSF.fminconllsCOMS(GibsonLaniPSF,...
%    GLMPP,PupilPlane);

CalibratePSF=CalibratePSF.genmodelarray(GibsonLaniPSF,GLMPP);
figure('Name','data and fit');
for ii=1:size(CalibratePSF.MODELarray,3)
    subplot(1,2,1)
    imagesc(CalibratePSF.SIGNALarray(:,:,ii)), axis image, colormap winter;
    colorbar;
    subplot(1,2,2)
    imagesc(CalibratePSF.MODELarray(:,:,ii)), axis image, colormap winter;
    colorbar;
    pause(0.1);
end