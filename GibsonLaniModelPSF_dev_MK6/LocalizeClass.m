classdef LocalizeClass
    properties
        DATAarray
        BACKGROUNDarray
        FITBACKGROUND
        variarray
        oiarray
        giarray
        xxxarray
        termAarray
        MODELarray
        Nphoton
        nllarray
        LookupTable
        codex
        ftvl0
        ftvl1
        iternum
    end
    methods
        function out = init(Localize,DATAin,BACKGROUNDin,variin,oiin,giin,...
                LookupTableStructure,iternumin)
            Localize.LookupTable=single(LookupTableStructure.LookupTable);
            Localize.codex=single(LookupTableStructure.codex);
            Localize.DATAarray=single(repmat(DATAin,1,1,size(Localize.codex,2)));
            Localize.BACKGROUNDarray=single(repmat(BACKGROUNDin,1,1,size(Localize.codex,2)));
            Localize.variarray=single(repmat(variin,1,1,size(Localize.codex,2)));
            Localize.oiarray=single(repmat(oiin,1,1,size(Localize.codex,2)));
            Localize.giarray=single(repmat(giin,1,1,size(Localize.codex,2)));
            Localize.Nphoton=single(...
                sum(sum(((DATAin-oiin)./giin)-BACKGROUNDin)));
            Localize.nllarray=zeros(1,size(Localize.codex,2));
            Localize.ftvl0=zeros(4,1);
            Localize.ftvl1=zeros(4,1);
            Localize.iternum=single(iternumin);
            
            out = Localize;
        end
        function out = guess(Localize,GLMPP)            
            Localize.termAarray=Localize.variarray./((Localize.giarray).^2);
            Localize.xxxarray=(Localize.DATAarray-Localize.oiarray)./...
                Localize.giarray+Localize.termAarray;
            
            Localize.MODELarray=...
                (Localize.Nphoton*ones(size(Localize.LookupTable))).*...
                Localize.LookupTable+Localize.BACKGROUNDarray;
            
            Localize.nllarray(:)=...
                reshape(...
                    sum(sum(...
                    (Localize.MODELarray+Localize.termAarray)-...
                    Localize.xxxarray.*log(...
                    Localize.MODELarray+Localize.termAarray)...
                    )),...
                    1,size(Localize.codex,2));
                
            Localize.ftvl0=...
                Localize.codex(:,(Localize.nllarray==min(Localize.nllarray)));
            
            %[Localize.ftvl0(1),Localize.ftvl0(2)]=COM(...
            %    (Localize.DATAarray(:,:,1)...
            %    -Localize.oiarray(:,:,1))./...
            %    Localize.giarray(:,:,1)-...
            %    Localize.BACKGROUNDarray(:,:,1),...
            %   GLMPP.PixelSize,GLMPP.M);
                        
            out = Localize;
            
        end
        function out = data(Localize,GibsonLaniPSF)
            GibsonLaniPSF=GibsonLaniPSF.fminsearchxyzllsCMOS(...
                Localize.DATAarray(:,:,1),...
                Localize.BACKGROUNDarray(:,:,1),...
                Localize.variarray(:,:,1),...
                Localize.giarray(:,:,1),...
                Localize.oiarray(:,:,1),...
                Localize.ftvl0,...
                Localize.iternum);
            
            Localize.ftvl1(1)=GibsonLaniPSF.xyzl(1);
            Localize.ftvl1(2)=GibsonLaniPSF.xyzl(2);
            Localize.ftvl1(3)=GibsonLaniPSF.xyzl(3);
            Localize.ftvl1(4)=GibsonLaniPSF.xyzl(4);
            
            out = Localize;
        end
        function out = datavarblurxyzlb(Localize,GibsonLaniPSF)
            GibsonLaniPSF=GibsonLaniPSF.fminsearchvarblurxyzllblsCMOS(...
                Localize.DATAarray(:,:,1),...
                Localize.BACKGROUNDarray(:,:,1),...
                Localize.variarray(:,:,1),...
                Localize.giarray(:,:,1),...
                Localize.oiarray(:,:,1),...
                Localize.ftvl0,...
                Localize.iternum);
            
            Localize.ftvl1(1)=GibsonLaniPSF.xyzl(1);
            Localize.ftvl1(2)=GibsonLaniPSF.xyzl(2);
            Localize.ftvl1(3)=GibsonLaniPSF.xyzl(3);
            Localize.ftvl1(4)=GibsonLaniPSF.xyzl(4);
            Localize.FITBACKGROUND=GibsonLaniPSF.BACKGROUND;
            Localize.Nphoton=GibsonLaniPSF.Nphoton;
            
            out = Localize;
        end
    end
end
            
            
            

            
        