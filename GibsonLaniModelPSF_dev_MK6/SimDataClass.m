classdef SimDataClass
    properties
        PIXSTATFileName
        gi
        vari
        oi
        GLMPPFileName
        PupilPlaneFileName
        
        LTFilterSetFileName
        LTPhasePlateFileName
        LTFilterSet
        LTPhasePlate
        
        SpectralStructureFileName
        SpectralStructure
        fluortags
        SpecCategories
        
        filterwidthvecin
        mincutoffin
        maxcutoffin
        deltax0
        deltay0
        deltazp
        x0range
        y0range
        zprange
        x0
        y0
        zp
        
        xoffset
        yoffset
        thetaoffset
        liveactuators
        actuatorheights
        spectralsignalarray
        pixelconfusion
        permutationarray
        P_con_ml
        optimalactuatorheights
        OptimalPupilPlaneFileName
        
        bgmean
        bgsigma
        Nphotonvec
        Nphoton 
        spectralresults
        numbertrials
        type
        filtersettruthtablearray
        phaseplatetruthtablearray
        
    end
    methods
        function out = init(SimData,PIXSTATFileNamein,GLMPPFileNamein,...
                PupilPlaneFileNamein,bgmeanin,bgsigmain,Nphotonvecin,...
                fluortagsin, varargin)
            SimData.PIXSTATFileName = PIXSTATFileNamein;
            
            SimData.GLMPPFileName = GLMPPFileNamein;
            GLMPP=GLMPPClass;
            GLMPP=GLMPP.init('GLMPP_temp.csv');
            GLMPP=GLMPP.update(SimData.GLMPPFileName);
            GLMPP=GLMPP.write('GLMPP_temp.csv');
            
            SimData.PupilPlaneFileName = PupilPlaneFileNamein;
            SimData.bgmean = single(bgmeanin);
            SimData.bgsigma = single(bgsigmain);
            SimData.Nphotonvec = single(Nphotonvecin);
            SimData.numbertrials = 100;
            SimData.spectralresults = [];
            SimData.fluortags = fluortagsin;
            
            SimData.gi=cast(...
                imread(strcat(SimData.PIXSTATFileName,'/','gain_map.tiff')),...
                'single');
            SimData.gi=SimData.gi(1:GLMPP.K,1:GLMPP.K);
            
            SimData.vari=cast(...
                imread(strcat(SimData.PIXSTATFileName,'/','variance_map.tiff')),...
                'single');
            SimData.vari=SimData.vari(1:GLMPP.K,1:GLMPP.K);
            
            SimData.oi=cast(...
                imread(strcat(SimData.PIXSTATFileName,'/','offset_map.tiff')),...
                'single');
            SimData.oi=SimData.oi(1:GLMPP.K,1:GLMPP.K);
            
            if (~isempty(varargin))
                numbervarin = numel(varargin);
                switch numbervarin
                    case 2
                        SimData.LTFilterSetFileName = varargin{1};
                        SimData.LTFilterSet = load(...
                            SimData.LTFilterSetFileName);
                        SimData.LTPhasePlateFileName = varargin{2};
                        SimData.LTPhasePlate = load(...
                            SimData.LTPhasePlateFileName);
                        
                    case 10
                        SimData.SpectralStructureFileName = varargin{1};
                        SimData.filterwidthvecin = single(varargin{2});
                        SimData.mincutoffin = single(varargin{3});
                        SimData.maxcutoffin = single(varargin{4});
                        SimData.deltax0 = single(varargin{5});
                        SimData.deltay0 = single(varargin{6});
                        SimData.deltazp = single(varargin{7});
                        SimData.x0range = single(varargin{8});
                        SimData.y0range = single(varargin{9});
                        SimData.zprange = single(varargin{10});
                        
                    case 15
                        SimData.SpectralStructureFileName = varargin{1};
                        SimData.filterwidthvecin = single(varargin{2});
                        SimData.mincutoffin = single(varargin{3});
                        SimData.maxcutoffin = single(varargin{4});
                        SimData.deltax0 = single(varargin{5});
                        SimData.deltay0 = single(varargin{6});
                        SimData.deltazp = single(varargin{7});
                        SimData.x0range = single(varargin{8});
                        SimData.y0range = single(varargin{9});
                        SimData.zprange = single(varargin{10});
                        SimData.xoffset = single(varargin{11});
                        SimData.yoffset = single(varargin{12});
                        SimData.thetaoffset = single(varargin{13});
                        SimData.liveactuators = single(varargin{14});
                        SimData.actuatorheights = single(varargin{15});
                        SimData.spectralsignalarray= zeros(...
                            GLMPP.K, GLMPP.K, numel(SimData.fluortags));
                        SimData.permutationarray=nchoosek(...
                            1:1:numel(SimData.fluortags),2);
                        SimData.P_con_ml=zeros(GLMPP.K,GLMPP.K,...
                            size(SimData.permutationarray,1));
                        
                end
            end
            
            out = SimData;
        end
        function out = createfiltersLookupTable(SimData,LTFilterSetFileNamein)
            GLMPP=GLMPPClass;
            GLMPP=GLMPP.init('GLMPP_temp.csv');
            GLMPP=GLMPP.update(SimData.GLMPPFileName);
            GLMPP=GLMPP.write('GLMPP_temp.csv');
            
            PupilPlane=PupilPlaneClass;
            PupilPlane=PupilPlane.init('PupilPlane_temp.csv');
            %PupilPlane=PupilPlane.update(SimData.PupilPlaneFileName);
            %PupilPlane=PupilPlane.write('PupilPlane_temp.csv');
            
            GibsonLaniPSF=GibsonLaniPSFClass;
            GibsonLaniPSF=GibsonLaniPSF.init(...
                GLMPP,PupilPlane);
            
            LookupTable=LookupTableClass;
            LookupTable=LookupTable.init(...
                            SimData.SpectralStructureFileName,...
                            SimData.fluortags,GLMPP);
            LookupTable=LookupTable.updategrid(SimData.deltax0,...
                SimData.deltay0,...
                SimData.deltazp,...
                SimData.x0range,...
                SimData.y0range,...
                SimData.zprange);
            
            display(LookupTable)
            
            LookupTable=LookupTable.filtersetstandardcreate(...
                GLMPP,GibsonLaniPSF, SimData.filterwidthvecin,...
                SimData.mincutoffin, SimData.maxcutoffin);   
            
            SimData.LTFilterSetFileName=LookupTable.write(LTFilterSetFileNamein);
            
            out = SimData;
        end
        function out = createphaseplateLookupTable(SimData,LTPhasePlateFileNamein)
            GLMPP=GLMPPClass;
            GLMPP=GLMPP.init('GLMPP_temp.csv');
            GLMPP=GLMPP.update(SimData.GLMPPFileName);
            GLMPP=GLMPP.write('GLMPP_temp.csv');
            
            PupilPlane=PupilPlaneClass;
            PupilPlane=PupilPlane.init('PupilPlane_temp.csv');
            PupilPlane=PupilPlane.update(SimData.PupilPlaneFileName);
            PupilPlane=PupilPlane.write('PupilPlane_temp.csv');
            
            GibsonLaniPSF=GibsonLaniPSFClass;
            GibsonLaniPSF=GibsonLaniPSF.init(...
                GLMPP,PupilPlane);
            
            LookupTable=LookupTableClass;
            LookupTable=LookupTable.init(...
                            SimData.SpectralStructureFileName,...
                            SimData.fluortags,GLMPP);
            LookupTable=LookupTable.updategrid(SimData.deltax0,...
                SimData.deltay0,...
                SimData.deltazp,...
                SimData.x0range,...
                SimData.y0range,...
                SimData.zprange);
              
            LookupTable=LookupTable.spectralcreate(GibsonLaniPSF);   
            SimData.LTPhasePlateFileName=LookupTable.write(LTPhasePlateFileNamein);
            
            out = SimData;
        end
        function out = createCUSTOMphaseplateLookupTable(SimData,LTPhasePlateFileNamein)
            
            GLMPP=GLMPPClass;
            GLMPP=GLMPP.init('GLMPP_temp.csv');
            GLMPP=GLMPP.update(SimData.GLMPPFileName);
            GLMPP=GLMPP.write('GLMPP_temp.csv');
            
            PupilPlane=PupilPlaneClass;
            PupilPlane=PupilPlane.init('PupilPlane_temp.csv');
            PupilPlane.THETAplate=XPP(PupilPlane.rho,...
                PupilPlane.theta,...
                SimData.xoffset,...
                SimData.yoffset,...
                SimData.thetaoffset,...
                SimData.liveactuators,...
                SimData.actuatorheights);
            
            %PupilPlane=PupilPlane.update(SimData.PupilPlaneFileName);
            %PupilPlane=PupilPlane.write('PupilPlane_temp.csv');
            
            GibsonLaniPSF=GibsonLaniPSFClass;
            GibsonLaniPSF=GibsonLaniPSF.init(...
                GLMPP,PupilPlane);
            
            LookupTable=LookupTableClass;
            LookupTable=LookupTable.init(...
                            SimData.SpectralStructureFileName,...
                            SimData.fluortags,GLMPP);
            LookupTable=LookupTable.updategrid(SimData.deltax0,...
                SimData.deltay0,...
                SimData.deltazp,...
                SimData.x0range,...
                SimData.y0range,...
                SimData.zprange);         
            
            LookupTable=LookupTable.spectralcreate(GibsonLaniPSF);
            
            SimData.LTPhasePlateFileName=LookupTable.write(LTPhasePlateFileNamein);
            
            out = SimData;
        end
        function out = PIXCON_OBJ(SimData,PupilPlane,...
                GibsonLaniPSF, x0in, y0in, zpin, Nphotonin,...
                actuatorheightsin)
            
            SimData.x0 = single(x0in);
            SimData.y0 = single(y0in);
            SimData.zp = single(zpin);
            SimData.Nphoton = single(Nphotonin);
            SimData.actuatorheights = single(actuatorheightsin);
                
            PupilPlane.THETAplate=XPP(PupilPlane.rho,...
                PupilPlane.theta,...
                SimData.xoffset,...
                SimData.yoffset,...
                SimData.thetaoffset,...
                SimData.liveactuators,...
                SimData.actuatorheights);
            GibsonLaniPSF=GibsonLaniPSF.updatePupilPlane(PupilPlane);
            for ii=1:numel(SimData.fluortags)
                GibsonLaniPSF = GibsonLaniPSF.spectralintens(...
                    x0in,...
                    y0in,...
                    zpin,...
                    SimData.SpectralStructure.(SimData.fluortags{ii}).Wavelength,...
                    SimData.SpectralStructure.(SimData.fluortags{ii}).Spectrum,...
                    SimData.SpectralStructure.(SimData.fluortags{ii}).Filter);
                SimData.spectralsignalarray(:,:,ii)=imgaussfilt(GibsonLaniPSF.spectralPSF,1);
                SimData.spectralsignalarray(:,:,ii)=SimData.Nphoton*...
                    SimData.spectralsignalarray(:,:,ii)./sum(sum(SimData.spectralsignalarray(:,:,ii)));
            end
        
            for ii=1:size(SimData.P_con_ml,3)
                mm=SimData.permutationarray(ii,1);
                ll=SimData.permutationarray(ii,2);
                SimData.P_con_ml(:,:,ii)=...
                    exp(-(SimData.spectralsignalarray(:,:,mm)+SimData.spectralsignalarray(:,:,ll))).*...
                    besseli(0,2*sqrt(SimData.spectralsignalarray(:,:,mm).*SimData.spectralsignalarray(:,:,ll)));
            end
            
            out = (1/nchoosek(numel(SimData.fluortags),2))*...
                sum(sum(sum(SimData.P_con_ml)));
        end
        function out = PIXCON_MIN(SimData, OptimalPupilPlaneFileNamein,...
                x0in, y0in, zpin, Nphotonin, lbin, ubin, type)
            SimData.OptimalPupilPlaneFileName=OptimalPupilPlaneFileNamein;
            
            GLMPP=GLMPPClass;
            GLMPP=GLMPP.init('GLMPP_temp.csv');
            
            PupilPlane=PupilPlaneClass;
            PupilPlane=PupilPlane.init('PupilPlane_temp.csv');
            
            GibsonLaniPSF=GibsonLaniPSFClass;
            GibsonLaniPSF=GibsonLaniPSF.init(...
                GLMPP,PupilPlane);
            
            SimData.SpectralStructure=...
                load(SimData.SpectralStructureFileName);
            
            fun=@(actuatorheightsin)SimData.PIXCON_OBJ(PupilPlane,...
                GibsonLaniPSF, x0in, y0in, zpin, Nphotonin,...
                actuatorheightsin);
            
            switch type
                case 'fminsearch'
                    options = optimset('PlotFcns','optimplotfval');
                    
                    SimData.optimalactuatorheights=...
                        fminsearch(fun,zeros(size(lbin)),options);
                    
                case 'pso'
            
                    options = optimoptions(@particleswarm,...
                        'PlotFcn', 'pswplotbestf', 'HybridFcn', @fmincon,...
                        'UseParallel', true, 'SwarmSize', 100);

                    SimData.optimalactuatorheights=...
                        particleswarm(fun,numel(lbin),lbin,ubin,options);
                    
                case 'ga'
                    
                    options = optimoptions(@ga,...
                        'PlotFcn', @gaplotbestf, 'HybridFcn', @fmincon,...
                        'UseParallel', true);
                    
                    SimData.optimalactuatorheights=...
                        ga(fun,numel(lbin),[],[],[],[],lbin,ubin,[],options);
            end
            
            PupilPlane.THETAplate=XPP(PupilPlane.rho,...
                PupilPlane.theta,...
                SimData.xoffset,...
                SimData.yoffset,...
                SimData.thetaoffset,...
                SimData.liveactuators,...
                SimData.optimalactuatorheights);
            
            PupilPlane.write(SimData.OptimalPupilPlaneFileName);
            
            out = SimData;           
        end
        function out = spectralMonteCarlo(SimData,numbertrialsin,typein)
            SimData.numbertrials = numbertrialsin;
            SimData.type = typein;
            switch SimData.type
                case 'FilterSet'
                    SimData.LTFilterSet = load(SimData.LTFilterSetFileName);
                    LTtemp = SimData.LTFilterSet.LookupTable;
                    codextemp = SimData.LTFilterSet.codex(4,:);
                    SimData.SpecCategories = unique(codextemp);
                    
                    gitemp = SimData.gi;
                    varitemp = SimData.vari;
                    oitemp = SimData.oi;
                    
                    oitemp=repmat(oitemp,[1,(size(LTtemp(:,:,1),2)/size(LTtemp(:,:,1),1))]);
                    varitemp=repmat(varitemp,[1,(size(LTtemp(:,:,1),2)/size(LTtemp(:,:,1),1))]);
                    gitemp=repmat(gitemp,[1,(size(LTtemp(:,:,1),2)/size(LTtemp(:,:,1),1))]);
                    
                    bgsigmatemp = SimData.bgsigma;
                    bgmeantemp = SimData.bgmean;
                    Nphotonvectemp = SimData.Nphotonvec;
                    numbertrialstemp = SimData.numbertrials;
                    
                    SimData.spectralresults=[];
                    spectralresultstemp = cell(1,numel(SimData.Nphotonvec));
                    
                    parfor ii=1:numel(Nphotonvectemp)
                        spectralresultstemp{ii}=single(zeros(numbertrialstemp,2));
                        for jj=1:numbertrialstemp
                            trueidx=randi(size(LTtemp,3));
                            signal=Nphotonvectemp(ii)*LTtemp(:,:,trueidx);
                            truecolor=codextemp(trueidx);
                            background=normrnd(bgmeantemp,...
                                bgsigmatemp)*ones(size(signal));
                            signal=repmat((poissrnd(signal+background)),...
                                [1,1,numel(codextemp)]);
                            termA=repmat(varitemp./(gitemp.^2),...
                                [1,1,numel(codextemp)]);
                            xxx=signal+termA;
                            model=Nphotonvectemp(ii)*LTtemp+...
                                repmat(background,[1,1,numel(codextemp)]);
                            llscmos=sum(sum(...
                                (model+termA)-xxx.*log(abs(model+termA))...
                                ));
                            estcolor=codextemp(llscmos==min(min(llscmos)));
                            spectralresultstemp{ii}(jj,1)=single(min(truecolor));
                            spectralresultstemp{ii}(jj,2)=single(min(estcolor));
                        end
                    end
                    SimData.spectralresults=spectralresultstemp;
                case 'PhasePlate'
                    SimData.LTPhasePlate = load(SimData.LTPhasePlateFileName);
                    LTtemp = SimData.LTPhasePlate.LookupTable;
                    codextemp = SimData.LTPhasePlate.codex(4,:);
                    SimData.SpecCategories = unique(codextemp);
                    
                    gitemp = SimData.gi;
                    varitemp = SimData.vari;
                    oitemp = SimData.oi;
                                        
                    bgsigmatemp = SimData.bgsigma;
                    bgmeantemp = SimData.bgmean;
                    Nphotonvectemp = SimData.Nphotonvec;
                    numbertrialstemp = SimData.numbertrials;
                    
                    SimData.spectralresults=[];
                    spectralresultstemp = cell(1,numel(SimData.Nphotonvec));
                    
                    parfor ii=1:numel(Nphotonvectemp)
                        spectralresultstemp{ii}=single(zeros(numbertrialstemp,2));
                        for jj=1:numbertrialstemp
                            trueidx=randi(size(LTtemp,3));
                            signal=Nphotonvectemp(ii)*LTtemp(:,:,trueidx);
                            truecolor=codextemp(trueidx);
                            background=abs(normrnd(bgmeantemp,...
                                bgsigmatemp)*ones(size(signal)));
                            signal=repmat((poissrnd(signal+background)),...
                                [1,1,numel(codextemp)]);
                            termA=repmat(varitemp./(gitemp.^2),...
                                [1,1,numel(codextemp)]);
                            xxx=signal+termA;
                            model=Nphotonvectemp(ii)*LTtemp+...
                                repmat(background,[1,1,numel(codextemp)]);
                            llscmos=sum(sum(...
                                (model+termA)-xxx.*log(abs(model+termA))...
                                ));
                            estcolor=codextemp(llscmos==min(min(llscmos)));
                            spectralresultstemp{ii}(jj,1)=single(min(truecolor));
                            spectralresultstemp{ii}(jj,2)=single(min(estcolor));
                        end
                    end
                    SimData.spectralresults=spectralresultstemp;
            end
            
            out = SimData;
        end
        function out = maketruthtablearray(SimData)
            switch SimData.type
                case 'FilterSet'
                    SimData.filtersettruthtablearray = single(zeros(...
                    numel(SimData.SpecCategories),...
                    numel(SimData.SpecCategories),...
                    numel(SimData.Nphotonvec)));
                    for ii = 1:numel(SimData.spectralresults)
                        for jj = 1:size(SimData.spectralresults{ii},1)
                            %truth corresponds to the column index
                            %estimation corersponds to the row index
                            SimData.filtersettruthtablearray(...
                                (SimData.spectralresults{ii}(jj,2)==...
                                SimData.SpecCategories),...
                                (SimData.spectralresults{ii}(jj,1)==...
                                SimData.SpecCategories),ii)=...
                            SimData.filtersettruthtablearray(...
                                (SimData.spectralresults{ii}(jj,2)==...
                                SimData.SpecCategories),...
                                (SimData.spectralresults{ii}(jj,1)==...
                                SimData.SpecCategories),ii)+1;
                        end
                    end
                case 'PhasePlate'
                    SimData.phaseplatetruthtablearray = single(zeros(...
                    numel(SimData.SpecCategories),...
                    numel(SimData.SpecCategories),...
                    numel(SimData.Nphotonvec)));
                    for ii = 1:numel(SimData.spectralresults)
                        for jj = 1:size(SimData.spectralresults{ii},1)
                            %truth corresponds to the column index
                            %estimation corersponds to the row index
                            truthidx=find(SimData.spectralresults{ii}(jj,2)==...
                                SimData.SpecCategories);
                            estidx=find(SimData.spectralresults{ii}(jj,1)==...
                                SimData.SpecCategories);
                            SimData.phaseplatetruthtablearray(estidx,truthidx,ii)=...
                            SimData.phaseplatetruthtablearray(estidx,truthidx,ii)+1;
                        end
                    end
            end
            
            out = SimData;
       
        end         
    end
end
