classdef CalibratePSFClass
    properties
        DELTAtrange
        DELTAtstpsz
        DELTAtvec
        excitlambda
        beadtags
        numberpos
        framesatpos
        names
        numsamples
        DATAarray
        ROIarray
        BACKGROUNDarray
        SIGNALarray
        xxxarray
        giarray
        variarray
        oiarray
        lbr
        liveactuators
        actuatorheights
        rho
        theta
        Bstar
        THETAplate 
        THETAZernike 
        xoffset
        yoffset
        thetaoffset
        numberZweights
        Zweights
        Zmodes
        MODEL
        MODELarray
        Nphotonarray
        maindirectory
        filelisting
        filetemplate
        filechecklist
        gi
        vari
        oi
        x0
        y0
        zp
        lambda
        x0vec
        y0vec
        zpvec
        lambdavec
        xROICOMvec
        yROICOMvec
        SpecMax
        spectra
        filters
        wavelengths
        termA
        XPhasePlateThickness
    end
    methods
        function out = init(CalibratePSF,GLMPP,PupilPlane,...
                SpectralStructure,beadtagsin,excitlambdain,...
                numberposin,framesatposin)          
            CalibratePSF.numberpos=single(numberposin);
            CalibratePSF.DELTAtrange=single([-1,1]);
            CalibratePSF.DELTAtstpsz=single(0.05);
            CalibratePSF.DELTAtvec=linspace(...
                CalibratePSF.DELTAtrange(1),CalibratePSF.DELTAtrange(2),...
                CalibratePSF.numberpos);
            CalibratePSF.excitlambda=single(excitlambdain);
            CalibratePSF.beadtags=beadtagsin; 
            CalibratePSF.framesatpos=single(framesatposin);           
            CalibratePSF.names=fieldnames(SpectralStructure);
            CalibratePSF.numsamples=single(20);
            CalibratePSF.DATAarray=single(0.0);
            CalibratePSF.ROIarray=single(0.0);
            CalibratePSF.BACKGROUNDarray=single(0.0);
            CalibratePSF.SIGNALarray=single(0.0);
            CalibratePSF.xxxarray=single(0.0);
            CalibratePSF.MODEL=single(0.0);
            CalibratePSF.MODELarray=single(0.0);
            CalibratePSF.Nphotonarray=single(0.0);
            CalibratePSF.gi=single(0.0);
            CalibratePSF.vari=single(0.0);
            CalibratePSF.oi=single(0.0);
            CalibratePSF.lbr=single(0.0);
            CalibratePSF.liveactuators=single([0,1,1,0]);
            CalibratePSF.XPhasePlateThickness=0.96;
            CalibratePSF.actuatorheights=CalibratePSF.XPhasePlateThickness*CalibratePSF.liveactuators;
            
            CalibratePSF.rho=ones(GLMPP.K,GLMPP.K,numel(PupilPlane.rho));
            CalibratePSF.theta=ones(GLMPP.K,GLMPP.K,numel(PupilPlane.rho));
            CalibratePSF.Bstar=ones(GLMPP.K,GLMPP.K,numel(PupilPlane.rho));
            CalibratePSF.THETAplate=ones(GLMPP.K,GLMPP.K,numel(PupilPlane.rho));
            CalibratePSF.THETAZernike=ones(GLMPP.K,GLMPP.K,numel(PupilPlane.rho));
            for ii=1:numel(PupilPlane.rho)
                CalibratePSF.rho(:,:,ii)=CalibratePSF.rho(:,:,ii)*PupilPlane.rho(ii);
                CalibratePSF.theta(:,:,ii)=CalibratePSF.theta(:,:,ii)*PupilPlane.theta(ii);
                CalibratePSF.Bstar(:,:,ii)=CalibratePSF.Bstar(:,:,ii)*PupilPlane.Bstar(ii);
                CalibratePSF.THETAplate(:,:,ii)=CalibratePSF.THETAplate(:,:,ii)*PupilPlane.THETAplate(ii);
                CalibratePSF.THETAZernike(:,:,ii)=CalibratePSF.THETAZernike(:,:,ii)*PupilPlane.THETAZernike(ii);
            end
            
            CalibratePSF.xoffset=single(0.0);
            CalibratePSF.yoffset=single(0.0);
            CalibratePSF.thetaoffset=single(0.0);          
            CalibratePSF.Zweights=single(repmat(...
                [0;...%horizontal tilt
                0;...%vertical tilt  
                0;...%oblique astigmatism
                0;...%defocus 
                0;...%vertical astigmatism
                0;...%vertical trefoil
                0;...%vertical coma
                0;...%horizontal coma
                0;...%oblique trefoil
                0;...%oblique quadrafoil
                0;...%oblique secondary astigmatism
                0;...%primary sperical
                0;...%vertical secondary astigmatism
                0],....%vertical quadrafoil
                1,numel(PupilPlane.rho))...
                );
            CalibratePSF.numberZweights=size(CalibratePSF.Zweights,1);
            CalibratePSF.Zmodes=single(...
                zeros(...
                CalibratePSF.numberZweights,numel(PupilPlane.rho))...
                );
            for ii=1:CalibratePSF.numberZweights
                CalibratePSF.Zmodes(ii,:)=ZernikePolynomial(ii,...
                    PupilPlane.rho,PupilPlane.theta);
            end
            CalibratePSF.maindirectory=pwd;
            CalibratePSF.filelisting='nul';
            CalibratePSF.filetemplate='img000_000_000000_0000000000_0000000000_0.tif';
            CalibratePSF.filechecklist=char.empty;
            %img000_000_zpositionnumber_framenumber_probenumber_0
            %zpositionnumber <- filetemplate(12:17)
            %franenumber <-filetemplate(19:28)
            %probenumber <-filetemplate(30:39)

            CalibratePSF.x0=single(0.0);
            CalibratePSF.y0=single(0.0);
            CalibratePSF.zp=single(GLMPP.beadsz/2);
            CalibratePSF.lambda=single(0.0);
            CalibratePSF.x0vec=single(zeros(1,CalibratePSF.numberpos));
            CalibratePSF.y0vec=single(zeros(1,CalibratePSF.numberpos));
            CalibratePSF.zpvec=single(zeros(1,CalibratePSF.numberpos));
            CalibratePSF.lambdavec=single(zeros(1,numel(...
                CalibratePSF.beadtags)));
            CalibratePSF.xROICOMvec=single(0.0);
            CalibratePSF.yROICOMvec=single(0.0);
            CalibratePSF.spectra=zeros(numel(CalibratePSF.beadtags),...
                numel(SpectralStructure.(CalibratePSF.names{1}).Spectrum));
            
            CalibratePSF.filters=zeros(numel(CalibratePSF.beadtags),...
                numel(SpectralStructure.(CalibratePSF.names{1}).Filter));          
        
            CalibratePSF.wavelengths=zeros(numel(CalibratePSF.beadtags),...
                numel(SpectralStructure.(CalibratePSF.names{1}).Wavelength));    
            
            for ii=1:numel(CalibratePSF.beadtags)
                for jj=1:numel(CalibratePSF.names)
                    if strcmp(CalibratePSF.beadtags{ii},CalibratePSF.names{jj})
                        CalibratePSF.spectra(ii,:)=single(SpectralStructure.(CalibratePSF.beadtags{ii}).Spectrum);
                        CalibratePSF.filters(ii,:)=single(SpectralStructure.(CalibratePSF.beadtags{ii}).Filter);
                        CalibratePSF.wavelengths(ii,:)=single(SpectralStructure.(CalibratePSF.beadtags{ii}).Wavelength);
                    end
                end
            end
                        
            CalibratePSF.SpecMax=zeros(1,numel(CalibratePSF.beadtags));
            for ii=1:numel(CalibratePSF.beadtags)
                CalibratePSF.SpecMax(ii)=...
                    CalibratePSF.wavelengths(ii,(max(CalibratePSF.spectra(ii,:).*...
                    CalibratePSF.filters(ii,:))==(CalibratePSF.spectra(ii,:).*...
                    CalibratePSF.filters(ii,:))));
            end
            CalibratePSF.termA=single(0.0);
            out = CalibratePSF;
        end
        function out = importdata(CalibratePSF,DATAdir,PIXSTATdir)           
            CalibratePSF.gi=cast(imread(fullfile(PIXSTATdir,'gain_map.tiff')),'single');
            CalibratePSF.vari=cast(imread(fullfile(PIXSTATdir,'variance_map.tiff')),'single');
            CalibratePSF.oi=cast(imread(fullfile(PIXSTATdir,'offset_map.tiff')),'single');
            
            CalibratePSF.DATAarray=single(...
                zeros(size(CalibratePSF.gi,1),size(CalibratePSF.gi,2),...
                numel(CalibratePSF.beadtags)*CalibratePSF.numberpos));
            
            cd(DATAdir);
            CalibratePSF.filelisting=dir('*.tif');
            cd(CalibratePSF.maindirectory);
            
            for ii=1:CalibratePSF.numberpos
                for jj=1:numel(CalibratePSF.beadtags)

                    %filetemplate=img000_000_zpositionnumber_framenumber_probenumber_0

                    %zpositionnumber <- filetemplate(12:17)
                    zpositionnumber=strcat('00000',int2str(ii-1));
                    zpositionnumber=zpositionnumber(...
                        numel(zpositionnumber)-5:end);

                    %framenumber <-filetemplate(19:28)
                    % choose the last data frame from each position
                    framenumber=strcat('000000000',int2str(CalibratePSF.framesatpos-1));
                    framenumber=framenumber(...
                        numel(framenumber)-9:end);

                    %probenumber <-filetemplate(30:39)
                    probenumber=strcat('000000000',int2str(jj-1));
                    probenumber=probenumber(...
                        numel(probenumber)-9:end);
                    
                    filetemp=strcat('img000_000_',zpositionnumber,'_',...
                        framenumber,'_',probenumber,'_0');
                    
                    CalibratePSF.filechecklist=...
                        cat(1,CalibratePSF.filechecklist,filetemp);
                end
            end
            
            for ii=1:numel(CalibratePSF.filelisting)
                for jj=1:size(CalibratePSF.filechecklist,1)
                    if strcmp(CalibratePSF.filelisting(ii).name(12:39),...
                            CalibratePSF.filechecklist(jj,12:39))
                        probe=cast(str2double(...
                            CalibratePSF.filechecklist(jj,30:39)),'uint64');
                        position=cast(str2double(...
                            CalibratePSF.filechecklist(jj,12:17)),'uint64');
                        CalibratePSF.DATAarray(:,:,...
                            (position+1)+probe*...
                            cast(CalibratePSF.numberpos,'uint64'))=...
                            cast(imread(fullfile(DATAdir,...
                            CalibratePSF.filelisting(ii).name)),'single');
                    end    
                end
            end
            
            out = CalibratePSF;
        end
        function out = selectROI(CalibratePSF,GLMPP)
            h=figure('Name','select ROI');
            imagesc(CalibratePSF.DATAarray(:,:,1)),...
                axis image, colormap winter;
            [xx,yy]=getpts(gca);
            close(h);
            
            CalibratePSF.xROICOMvec=zeros(1,numel(xx));
            CalibratePSF.yROICOMvec=zeros(1,numel(yy));
            
            %initialize CalibratePSF.ROIarray
            CalibratePSF.ROIarray=single(zeros(...
                GLMPP.K*numel(xx),GLMPP.K,size(CalibratePSF.DATAarray,3)));
            CalibratePSF.BACKGROUNDarray=single(zeros(...
                GLMPP.K*numel(xx),GLMPP.K,size(CalibratePSF.DATAarray,3)));
            CalibratePSF.giarray=single(zeros(...
                GLMPP.K*numel(xx),GLMPP.K,size(CalibratePSF.DATAarray,3)));
            CalibratePSF.variarray=single(zeros(...
                GLMPP.K*numel(xx),GLMPP.K,size(CalibratePSF.DATAarray,3)));
            CalibratePSF.oiarray=single(zeros(...
                GLMPP.K*numel(xx),GLMPP.K,size(CalibratePSF.DATAarray,3)));
            CalibratePSF.MODELarray=single(zeros(...
                GLMPP.K*numel(xx),GLMPP.K,size(CalibratePSF.DATAarray,3)));
            CalibratePSF.Nphotonarray=single(zeros(numel(xx),...
                size(CalibratePSF.DATAarray,3)));
            CalibratePSF.xxxarray=single(zeros(...
                GLMPP.K*numel(xx),GLMPP.K,size(CalibratePSF.DATAarray,3)));
            CalibratePSF.termA=single(zeros(...
                GLMPP.K*numel(xx),GLMPP.K,size(CalibratePSF.DATAarray,3)));
            %multiple ROIs are fit simultaneously instead of one by one.
            %ROIs are stacked along the first dimension. the last three
            %two dimensional arrays in the the third direction are the
            %pixel statistic arrays associated with the ROI.
            
            for ii=1:numel(xx)            
                    %assign data ROI
                    CalibratePSF.ROIarray(GLMPP.K*(ii-1)+1:GLMPP.K*ii,...
                        1:GLMPP.K,:)=...
                        CalibratePSF.DATAarray(...
                        round(linspace(yy(ii)-(GLMPP.K/2),yy(ii)+(GLMPP.K/2),GLMPP.K)),...
                        round(linspace(xx(ii)-(GLMPP.K/2),xx(ii)+(GLMPP.K/2),GLMPP.K)),...
                        :);      
                for jj=1:size(CalibratePSF.DATAarray,3)    
                   %calculate background for each PSF for each frame
                   fullROI=...
                        CalibratePSF.DATAarray(...
                        round(linspace(yy(ii)-(GLMPP.K/2+1),yy(ii)+(GLMPP.K/2+1),GLMPP.K+2)),...
                        round(linspace(xx(ii)-(GLMPP.K/2+1),xx(ii)+(GLMPP.K/2+1),GLMPP.K+2)),...
                        jj);  
                    
                   fullgi=...
                        CalibratePSF.gi(...
                        round(linspace(yy(ii)-(GLMPP.K/2+1),yy(ii)+(GLMPP.K/2+1),GLMPP.K+2)),...
                        round(linspace(xx(ii)-(GLMPP.K/2+1),xx(ii)+(GLMPP.K/2+1),GLMPP.K+2)));  
                    
                   fulloi=...
                        CalibratePSF.oi(...
                        round(linspace(yy(ii)-(GLMPP.K/2+1),yy(ii)+(GLMPP.K/2+1),GLMPP.K+2)),...
                        round(linspace(xx(ii)-(GLMPP.K/2+1),xx(ii)+(GLMPP.K/2+1),GLMPP.K+2)));  
                    
                    CalibratePSF=CalibratePSF.ROIbackground(fullROI,...
                        fullgi,fulloi); 
                    
                    %assign corresponding background ROI.
                    %not in terms of photon count
                    CalibratePSF.BACKGROUNDarray(GLMPP.K*(ii-1)+1:GLMPP.K*ii,...
                        1:GLMPP.K,jj)=CalibratePSF.lbr;
                    
                end                   
                %assign corresponding gi ROI
                CalibratePSF.giarray(GLMPP.K*(ii-1)+1:GLMPP.K*ii,...
                    1:GLMPP.K,1:size(CalibratePSF.DATAarray,3))=...
                    repmat(CalibratePSF.gi(...
                    round(linspace(yy(ii)-(GLMPP.K/2),yy(ii)+(GLMPP.K/2),GLMPP.K)),...
                    round(linspace(xx(ii)-(GLMPP.K/2),xx(ii)+(GLMPP.K/2),GLMPP.K))),...
                    1,1,size(CalibratePSF.DATAarray,3));

                %assign corresponding vari ROI
                CalibratePSF.variarray(GLMPP.K*(ii-1)+1:GLMPP.K*ii,...
                    1:GLMPP.K,1:size(CalibratePSF.DATAarray,3))=...
                    repmat(CalibratePSF.vari(...
                    round(linspace(yy(ii)-(GLMPP.K/2),yy(ii)+(GLMPP.K/2),GLMPP.K)),...
                    round(linspace(xx(ii)-(GLMPP.K/2),xx(ii)+(GLMPP.K/2),GLMPP.K))),...
                    1,1,size(CalibratePSF.DATAarray,3));

                %assign corresponding oi ROI
                CalibratePSF.oiarray(GLMPP.K*(ii-1)+1:GLMPP.K*ii,...
                    1:GLMPP.K,1:size(CalibratePSF.DATAarray,3))=...
                    repmat(CalibratePSF.oi(...
                    round(linspace(yy(ii)-(GLMPP.K/2),yy(ii)+(GLMPP.K/2),GLMPP.K)),...
                    round(linspace(xx(ii)-(GLMPP.K/2),xx(ii)+(GLMPP.K/2),GLMPP.K))),...
                    1,1,size(CalibratePSF.DATAarray,3));        
                
            end
            
            %calculate the pixel statistic adjusted signal and xxx value.
            CalibratePSF.SIGNALarray=(CalibratePSF.ROIarray-CalibratePSF.oiarray)./CalibratePSF.giarray;
            
            CalibratePSF.xxxarray=CalibratePSF.SIGNALarray+CalibratePSF.variarray./((CalibratePSF.giarray).^2);
            for ii=1:numel(xx)
                for jj=1:size(CalibratePSF.DATAarray,3)
                    CalibratePSF.Nphotonarray(ii,jj)=sum(sum(...
                        CalibratePSF.SIGNALarray((ii-1)*GLMPP.K+1:ii*GLMPP.K,1:GLMPP.K,jj)-...
                        CalibratePSF.BACKGROUNDarray((ii-1)*GLMPP.K+1:ii*GLMPP.K,1:GLMPP.K,jj)...
                        ));                       
                end
            end
            
            %calculate the center of mass for each ROI
            for ii=1:numel(CalibratePSF.xROICOMvec)
                [CalibratePSF.xROICOMvec(ii),...
                    CalibratePSF.yROICOMvec(ii)]=...
                    COM(CalibratePSF.ROIarray(...
                    (ii-1)*GLMPP.K+1:ii*GLMPP.K,1:GLMPP.K,...
                    round(CalibratePSF.numberpos/2)),...
                    GLMPP.PixelSize,GLMPP.M);
            end
            out = CalibratePSF;
        end
        function out = ROIbackground(CalibratePSF,fullROI,...
                fullgi,fulloi)
            
            %background is in photon counts
            fullROI=(fullROI-fulloi)./fullgi;
            
            %calculate average of perimieter pixels of meanfullROI
            Kfull=size(fullROI,2);
            background=((1/(4*size(fullROI,2)-2+4))*(sum(fullROI(1,1:Kfull))+...
                sum(fullROI(Kfull,1:Kfull))+sum(fullROI(2:(Kfull-1),1))+...
                sum(fullROI(2:(Kfull-1),Kfull))))*ones(Kfull-2);

            CalibratePSF.lbr = abs(background);
            
            out = CalibratePSF;
        end
        function out = setpupil(CalibratePSF,PupilPlane,...
                Zweightsin, xoffsetin, yoffsetin,...
                thetaoffsetin)
            %Zweightsin is a column vector
            CalibratePSF.Zweights(:,:)=single(...
                repmat(Zweightsin,1,numel(PupilPlane.rho))...
                );
            CalibratePSF.xoffset=single(xoffsetin);
            CalibratePSF.yoffset=single(yoffsetin);
            CalibratePSF.thetaoffset=single(thetaoffsetin);
            
            PupilPlane.THETAZernike=sum(...
                CalibratePSF.Zweights.*CalibratePSF.Zmodes,1);
            
            PupilPlane.THETAplate=XPP(PupilPlane.rho,...
                PupilPlane.theta,...
                CalibratePSF.xoffset,...
                CalibratePSF.yoffset,...
                CalibratePSF.thetaoffset,...
                CalibratePSF.liveactuators,...
                CalibratePSF.actuatorheights);
            
            for ii=1:numel(PupilPlane.rho)
                CalibratePSF.THETAplate(:,:,ii)=PupilPlane.THETAplate(ii);
                CalibratePSF.THETAZernike(:,:,ii)=PupilPlane.THETAZernike(ii);
            end
            
            out = CalibratePSF;
        end
        function out = genmodel(CalibratePSF,GibsonLaniPSF,...
                x0in, y0in, zpin, lambdain, DELTAtin)
                        
            GibsonLaniPSF.DELTAt(:)=single(DELTAtin);
            GibsonLaniPSF.THETAZernike(:)=CalibratePSF.THETAZernike(:);
            GibsonLaniPSF.THETAplate(:)=CalibratePSF.THETAplate(:);
            GibsonLaniPSF=GibsonLaniPSF.bead(x0in, y0in, zpin, lambdain); 
            CalibratePSF.MODEL=GibsonLaniPSF.PSF;
            
            out = CalibratePSF;
        end
        function out = genmodelarray(CalibratePSF,GibsonLaniPSF,GLMPP)
            for ii=1:numel(CalibratePSF.SpecMax)
                for jj=1:numel(CalibratePSF.xROICOMvec)
                    for kk=1:CalibratePSF.numberpos
                        CalibratePSF=CalibratePSF.genmodel(...
                            GibsonLaniPSF,CalibratePSF.x0vec(jj),...
                            CalibratePSF.y0vec(jj),CalibratePSF.zp,...
                            CalibratePSF.SpecMax(ii),...
                            CalibratePSF.DELTAtvec(kk));
                        
                        CalibratePSF.MODELarray((jj-1)*GLMPP.K+1:...
                        jj*GLMPP.K,1:GLMPP.K,kk+...
                        (ii-1)*CalibratePSF.numberpos)=...
                        CalibratePSF.Nphotonarray(jj,kk+...
                        (ii-1)*CalibratePSF.numberpos)*CalibratePSF.MODEL+...
                        CalibratePSF.BACKGROUNDarray((jj-1)*GLMPP.K+1:jj*GLMPP.K,...
                        1:GLMPP.K,kk+(ii-1)*CalibratePSF.numberpos);
                  
                    end
                end
            end
            
            out = CalibratePSF;
        end
        function out = ObjFunllsCMOS(CalibratePSF,GibsonLaniPSF,GLMPP,...
                PupilPlane,ftval)
            ftval=single(reshape(ftval,numel(ftval),1));
                        
            CalibratePSF=CalibratePSF.setpupil(PupilPlane,...
                ftval(1:CalibratePSF.numberZweights),...
                CalibratePSF.xoffset,...
                CalibratePSF.yoffset,...
                CalibratePSF.thetaoffset);
            
            CalibratePSF.x0vec=ftval(CalibratePSF.numberZweights+1:...
                CalibratePSF.numberZweights+numel(CalibratePSF.xROICOMvec));
            CalibratePSF.y0vec=ftval(CalibratePSF.numberZweights+numel(CalibratePSF.xROICOMvec)+1:...
                CalibratePSF.numberZweights+numel(CalibratePSF.xROICOMvec)+numel(CalibratePSF.yROICOMvec));...   
            CalibratePSF.zp=ftval(CalibratePSF.numberZweights+...
                numel(CalibratePSF.xROICOMvec)+...
                numel(CalibratePSF.yROICOMvec)+1);
           
            CalibratePSF=CalibratePSF.genmodelarray(GibsonLaniPSF,GLMPP);
            
            CalibratePSF.termA=CalibratePSF.MODELarray+...
                (CalibratePSF.variarray./((...
                 CalibratePSF.giarray).^2));
            
            out = double(sum(sum(sum(...
                 abs(CalibratePSF.termA)-CalibratePSF.xxxarray.*...
                 log(abs(CalibratePSF.termA))))));
            
        end   
        function out = gallsCOMS(CalibratePSF,GibsonLaniPSF,GLMPP,...
                PupilPlane)
            ftval0=double([CalibratePSF.Zweights(:,1);...
                CalibratePSF.xROICOMvec';...
                CalibratePSF.yROICOMvec';...
                CalibratePSF.zp]);
            
            rng('default')
            
            fun=@(ftval)CalibratePSF.ObjFunllsCMOS(GibsonLaniPSF,...
                GLMPP,PupilPlane,ftval);
            
            options=optimoptions(@ga,'PlotFcn', {@gaplotbestf, @gaplotrange}, 'UseParallel',true);
            
            fminconOptions = optimoptions(@fmincon,'Display','iter','Algorithm','sqp',...
                'UseParallel',true,'MaxIterations',500,'MaxFunctionEvaluations',1000,...
                'OptimalityTolerance',1e-20,'StepTolerance',1e-20);
            
            options = optimoptions(options,'HybridFcn',{@fmincon, fminconOptions});
            
            ftval1=ga(fun,numel(ftval0),[],[],[],[],[-0.3*ones(1,numel(ftval0)-1),-0.5, 0.5, 0.5],...
                                                    [ 0.3*ones(1,numel(ftval0)-1), 0.5, 0.5, 0.5],[],options);
            reshape(ftval1,numel(ftval1),1);
            poolobj = gcp('nocreate');
            delete(poolobj);
            ftval1=single(ftval1);
            
            CalibratePSF=CalibratePSF.setpupil(PupilPlane,...
                ftval1(1:CalibratePSF.numberZweights),...
                CalibratePSF.xoffset,...
                CalibratePSF.yoffset,...
                CalibratePSF.thetaoffset);
            
            CalibratePSF.x0vec=ftval1(CalibratePSF.numberZweights+1:...
                CalibratePSF.numberZweights+numel(CalibratePSF.xROICOMvec));
            CalibratePSF.y0vec=ftval1(CalibratePSF.numberZweights+numel(CalibratePSF.xROICOMvec)+1:...
                CalibratePSF.numberZweights+numel(CalibratePSF.xROICOMvec)+numel(CalibratePSF.yROICOMvec));...   
            CalibratePSF.zp=ftval1(CalibratePSF.numberZweights+...
                numel(CalibratePSF.xROICOMvec)+...
                numel(CalibratePSF.yROICOMvec)+1);
            
            out = CalibratePSF;
        end
        function out = fminsearchllsCOMS(CalibratePSF,GibsonLaniPSF,GLMPP,...
                PupilPlane)
            ftval0=double([CalibratePSF.Zweights(:,1);...
                CalibratePSF.xROICOMvec';...
                CalibratePSF.yROICOMvec';...
                CalibratePSF.zp]);
            
            fun=@(ftval)CalibratePSF.ObjFunllsCMOS(GibsonLaniPSF,...
                GLMPP,PupilPlane,ftval);
            
            options=optimset('PlotFcns',@optimplotfval,...
                'MaxFunEvals',10000*numel(ftval0),'MaxIter',10000,...
                'TolFun',1e-4,'TolX',1e-4);
            
            ftval1=fminsearch(fun,ftval0,options);
            reshape(ftval1,numel(ftval1),1);
            poolobj = gcp('nocreate');
            delete(poolobj);
            ftval1=single(ftval1);
            
            CalibratePSF=CalibratePSF.setpupil(PupilPlane,...
                ftval1(1:CalibratePSF.numberZweights),...
                CalibratePSF.xoffset,...
                CalibratePSF.yoffset,...
                CalibratePSF.thetaoffset);
            
            CalibratePSF.x0vec=ftval1(CalibratePSF.numberZweights+1:...
                CalibratePSF.numberZweights+numel(CalibratePSF.xROICOMvec));
            CalibratePSF.y0vec=ftval1(CalibratePSF.numberZweights+numel(CalibratePSF.xROICOMvec)+1:...
                CalibratePSF.numberZweights+numel(CalibratePSF.xROICOMvec)+numel(CalibratePSF.yROICOMvec));...   
            CalibratePSF.zp=ftval1(CalibratePSF.numberZweights+...
                numel(CalibratePSF.xROICOMvec)+...
                numel(CalibratePSF.yROICOMvec)+1);
            
            out = CalibratePSF;
        end
        function out = fminuncllsCOMS(CalibratePSF,GibsonLaniPSF,GLMPP,...
                PupilPlane)
            ftval0=double([CalibratePSF.Zweights(:,1);...
                CalibratePSF.xROICOMvec';...
                CalibratePSF.yROICOMvec';...
                CalibratePSF.zp]);
                        
            fun=@(ftval)CalibratePSF.ObjFunllsCMOS(GibsonLaniPSF,...
                GLMPP,PupilPlane,ftval);
            
            options=optimoptions(@fminunc,'Algorithm','quasi-newton',...
                'MaxFunctionEvaluations',10000,'MaxIterations',500,...
                'OptimalityTolerance',1e-20,'StepTolerance',1e-20,...
                'UseParallel',true);
            
            ftval1=fminunc(fun,ftval0,options);
            reshape(ftval1,numel(ftval1),1);
            poolobj = gcp('nocreate');
            delete(poolobj);
            ftval1=single(ftval1);
            
            CalibratePSF=CalibratePSF.setpupil(PupilPlane,...
                ftval1(1:CalibratePSF.numberZweights),...
                CalibratePSF.xoffset,...
                CalibratePSF.yoffset,...
                CalibratePSF.thetaoffset);
            
            CalibratePSF.x0vec=ftval1(CalibratePSF.numberZweights+1:...
                CalibratePSF.numberZweights+numel(CalibratePSF.xROICOMvec));
            CalibratePSF.y0vec=ftval1(CalibratePSF.numberZweights+numel(CalibratePSF.xROICOMvec)+1:...
                CalibratePSF.numberZweights+numel(CalibratePSF.xROICOMvec)+numel(CalibratePSF.yROICOMvec));...   
            CalibratePSF.zp=ftval1(CalibratePSF.numberZweights+...
                numel(CalibratePSF.xROICOMvec)+...
                numel(CalibratePSF.yROICOMvec)+1);
            
            out = CalibratePSF;
        end
        function out = fminconllsCOMS(CalibratePSF,GibsonLaniPSF,GLMPP,...
                PupilPlane)
            ftval0=double([CalibratePSF.Zweights(:,1);...
                CalibratePSF.xROICOMvec';...
                CalibratePSF.yROICOMvec';...
                CalibratePSF.zp]);       
            
            fun=@(ftval)CalibratePSF.ObjFunllsCMOS(GibsonLaniPSF,...
                GLMPP,PupilPlane,ftval);
            
            options=optimoptions(@fmincon,'Algorithm','interior-point',...
                'MaxFunctionEvaluations',10,'MaxIterations',500,...
                'OptimalityTolerance',1e-20,'StepTolerance',1e-20,...
                'UseParallel',true,'Display','iter');
            
            ftval1=fmincon(fun,ftval0,[],[],[],[],...
                [-0.3*ones(1,CalibratePSF.numberZweights),-0.5*ones(1,numel(CalibratePSF.xROICOMvec)),-0.5*ones(1,numel(CalibratePSF.xROICOMvec)), -0.5],...
                [ 0.3*ones(1,CalibratePSF.numberZweights), 0.5*ones(1,numel(CalibratePSF.xROICOMvec)), 0.5*ones(1,numel(CalibratePSF.xROICOMvec)),  0.5],[],options);
            reshape(ftval1,numel(ftval1),1);
            poolobj = gcp('nocreate');
            delete(poolobj);
            ftval1=single(ftval1);
            
            CalibratePSF=CalibratePSF.setpupil(PupilPlane,...
                ftval1(1:CalibratePSF.numberZweights),...
                CalibratePSF.xoffset,...
                CalibratePSF.yoffset,...
                CalibratePSF.thetaoffset);
            
            CalibratePSF.x0vec=ftval1(CalibratePSF.numberZweights+1:...
                CalibratePSF.numberZweights+numel(CalibratePSF.xROICOMvec));
            CalibratePSF.y0vec=ftval1(CalibratePSF.numberZweights+numel(CalibratePSF.xROICOMvec)+1:...
                CalibratePSF.numberZweights+numel(CalibratePSF.xROICOMvec)+numel(CalibratePSF.yROICOMvec));...   
            CalibratePSF.zp=ftval1(CalibratePSF.numberZweights+...
                numel(CalibratePSF.xROICOMvec)+...
                numel(CalibratePSF.yROICOMvec)+1);
                        
            out = CalibratePSF;
        end
    end
end

function ret=ZernikePolynomial(index,rho,theta)
%I'm using the OSA/ANSI standdard for the index here
%index=((n*(n+2)+m))/2 normalized such that the integral of the dquareo of
%the mode over rho and theta is equalt to pi in a pupil plane whos mazimum
%value in rho is one.
%
%unline teh theta functions for phase plates or deformable mirrors these
%THETA functions are not such that their max value is one. Therefore the
%max stroke should be taken to be one.
switch index
    case 0
        phase_profile=ones(size(rho));
    case 1
        phase_profile=2*(rho.*sin(theta));
    case 2
        phase_profile=2*(rho.*cos(theta));
    case 3
        phase_profile=sqrt(6)*(rho.^2.*sin(2*theta));
    case 4
        phase_profile=sqrt(3)*(2*rho.^2-1);
    case 5
        phase_profile=sqrt(6)*(rho.^2.*cos(2*theta));
    case 6
        phase_profile=sqrt(8)*((rho.^3).*sin(3*theta));
    case 7
        phase_profile=sqrt(8)*((3*rho.^3-2*rho).*sin(theta));
    case 8
        phase_profile=sqrt(8)*((3*rho.^3-2*rho).*cos(theta));
    case 9
        phase_profile=sqrt(8)*((rho.^3).*cos(3*theta));
    case 10
        phase_profile=sqrt(10)*((rho.^4).*sin(4*theta));
    case 11
        phase_profile=sqrt(10)*((4*rho.^4-3*rho.^2).*sin(2*theta));
    case 12
        phase_profile=sqrt(5)*(6*rho.^4-6*rho.^2+1);
    case 13
        phase_profile=sqrt(10)*((4*rho.^4-3*rho.^2).*cos(2*theta));
    case 14
        phase_profile=sqrt(10)*((rho.^4).*cos(4*theta));
    case 15
        phase_profile=sqrt(12)*((rho.^5).*sin(5*theta));
    case 16
        phase_profile=sqrt(12)*((5*rho.^5-4*rho.^3).*sin(3*theta));
    case 17
        phase_profile=sqrt(12)*((10*rho.^5-12*rho.^3+3*rho).*sin(theta));
    case 18
        phase_profile=sqrt(12)*((10*rho.^5-12*rho.^3+3*rho).*cos(theta));
    case 19
        phase_profile=sqrt(12)*((5*rho.^5-4*rho.^3).*cos(3*theta));
    case 20
        phase_profile=sqrt(12)*((rho.^5).*cos(5*theta));
    case 21
        phase_profile=sqrt(14)*((rho.^6).*sin(6*theta));
    case 22
        phase_profile=sqrt(14)*((6*rho.^6-5*rho.^4).*sin(4*theta));
    case 23
        phase_profile=sqrt(14)*((15*rho.^6-20*rho.^4+6*rho.^2).*sin(2*theta));
    case 24
        phase_profile=sqrt(7)*(20*rho.^6-30*rho.^4+12*rho.^2-1);
    case 25
        phase_profile=sqrt(14)*((15*rho.^6-20*rho.^4+6*rho.^2).*cos(2*theta));
    case 26
        phase_profile=sqrt(14)*((6*rho.^6-5*rho.^4).*cos(4*theta));
    case 27
        phase_profile=sqrt(14)*((rho.^6).*cos(6*theta));
    case 28
        phase_profile=4*((rho.^7).*sin(7*theta));
    case 29
        phase_profile=4*((7*rho.^7-6*rho.^5).*sin(5*theta));
    case 30
        phase_profile=4*((21*rho.^7-30*rho.^5+10*rho.^3).*sin(3*theta));
    case 31
        phase_profile=4*((35*rho.^7-60*rho.^5+30*rho.^3-4*rho).*sin(theta));
    case 32
        phase_profile=4*((35*rho.^7-60*rho.^5+30*rho.^3-4*rho).*cos(theta));
    case 33
        phase_profile=4*((21*rho.^7-30*rho.^5+10*rho.^3).*cos(3*theta));
    case 34
        phase_profile=4*((7*rho.^7-6*rho.^5).*cos(5*theta));
    case 35
        phase_profile=4*((rho.^7).*sin(7*theta));
    otherwise
        error('subscript must be a positive integer less than or equatl to 35');
end
ret=phase_profile;
end

function ret=XPP(rho, theta, xoffset, yoffset, thetaoffset,...
    liveactuators, actuatorheights)

    liveactuators=reshape(liveactuators,...
        sqrt(numel(liveactuators)),...
        sqrt(numel(liveactuators)));

    actuatorheights=reshape(actuatorheights,...
        sqrt(numel(actuatorheights)),...
        sqrt(numel(actuatorheights)));

    x_pupil=rho.*cos(theta);
    y_pupil=rho.*sin(theta);

    liveactuators=rot90(liveactuators,2);
    actuatorheights=rot90(actuatorheights,2);

    rotation_matrix=[cos(thetaoffset),-sin(thetaoffset);
                     sin(thetaoffset),cos(thetaoffset)];
    %%%%% rotate pupil coordinates %%%%%
    %                                  % 
    pupil_coords=[reshape(x_pupil,[1,numel(x_pupil)]);...
                  reshape(y_pupil,[1,numel(y_pupil)])];
    pupil_coords_rot=zeros(size(pupil_coords));
    for ii=1:size(pupil_coords,2)
        pupil_coords_rot(:,ii)=rotation_matrix*pupil_coords(:,ii);
    end
    x_pupil_rot=reshape(pupil_coords_rot(1,:),[size(x_pupil,1),size(x_pupil,2)]);
    y_pupil_rot=reshape(pupil_coords_rot(2,:),[size(y_pupil,1),size(y_pupil,2)]);
    %                                  % 
    %%%%% rotate pupil coordinates %%%%%

    K_X_actuators=size(liveactuators,2);
    X_interactuator_distance=2/K_X_actuators;
    X_actuator_positions=linspace(-1+X_interactuator_distance/2,1-X_interactuator_distance/2,K_X_actuators);
    K_Y_actuators=size(liveactuators,1);
    Y_interactuator_distance=2/K_Y_actuators;
    Y_actuator_positions=linspace(-1+Y_interactuator_distance/2,1-Y_interactuator_distance/2,K_Y_actuators);
    [X_points,Y_points]=meshgrid(X_actuator_positions,Y_actuator_positions);

    alpha_x=100;
    alpha_y=100;
    phase_plate=zeros(size(y_pupil,1),size(x_pupil,2));
    for ii=1:size(X_points,1)
        for jj=1:size(Y_points,2)
            if liveactuators(ii,jj)==1           
                if (ii==1) && (jj==size(Y_points,1)) 
                    phase_plate=phase_plate+actuatorheights(ii,jj)*...
                        custom_sigmf(-(x_pupil_rot-xoffset),[alpha_x, X_points(ii,jj)-X_interactuator_distance/2]).*...
                        custom_sigmf(y_pupil_rot-yoffset,[alpha_y, -Y_points(ii,jj)-Y_interactuator_distance/2]);
                elseif (ii==size(X_points,2)) && (jj==1)
                    phase_plate=phase_plate+actuatorheights(ii,jj)*...
                        custom_sigmf((x_pupil_rot-xoffset),[alpha_x, -X_points(ii,jj)-X_interactuator_distance/2]).*...
                        custom_sigmf(-(y_pupil_rot-yoffset),[alpha_y, Y_points(ii,jj)-Y_interactuator_distance/2]);
                else
                    phase_plate=phase_plate+actuatorheights(ii,jj)*...
                        custom_sigmf(x_pupil_rot-xoffset,[alpha_x, -X_points(ii,jj)-X_interactuator_distance/2]).*...
                        custom_sigmf(-(x_pupil_rot-xoffset),[alpha_x, X_points(ii,jj)-X_interactuator_distance/2]).*...
                        custom_sigmf(y_pupil_rot-yoffset,[alpha_y, -Y_points(ii,jj)-Y_interactuator_distance/2]).*...
                        custom_sigmf(-(y_pupil_rot-yoffset),[alpha_y, Y_points(ii,jj)-Y_interactuator_distance/2]);     
                end
            end
        end
    end
    ret=phase_plate;
end

% this is a custom sigmf function 
function ret=custom_sigmf(x,params_vec)
    a=params_vec(1);
    c=params_vec(2);
    ret=1./(1+exp(-a*(x-c)));
end

