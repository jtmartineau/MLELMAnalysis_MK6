classdef GLMPPClass
    properties       
        NA
        M
        ns
        ni
        DELTAt
        zd
        zdstar
        np
        na
        beadsz
        K
        PixelSize
    end
    methods
        function out = init(GLMPP,FileName)        
            % intialize the properties of the object
            GLMPP.NA=single(1.2);
            GLMPP.M=single(67); % measured with AFM grid
            GLMPP.ns=single(1.38);
            GLMPP.ni=single(1.46);
            GLMPP.DELTAt=single(0.0);
            GLMPP.zd=single(142e3);
            GLMPP.zdstar=single(142e3);
            GLMPP.np=single(1.46);
            GLMPP.na=single(1.00);
            GLMPP.beadsz=single(0.100);
            GLMPP.K=single(20);
            GLMPP.PixelSize=single(6.5);
            
            % define cell array to write to .csv file
            GLMPPArray={'NA',GLMPP.NA;...
            'M',GLMPP.M;...
            'ns',GLMPP.ns;...
            'ni',GLMPP.ni;...
            'DELTAt',GLMPP.DELTAt;... % in um in image space
            'zd',GLMPP.zd;...
            'zdstar',GLMPP.zdstar;... % in um in image space
            'np',GLMPP.np;...
            'na',GLMPP.na;...
            'beadsz',single(0.200);...
            'K',single(20);...
            'PixelSize',single(6.5)}; % in um in object space: diameter of fluorescent bead
            
            % write cell array values to file
            GLMPPFileName=strcat(FileName);
            writecell(GLMPPArray,GLMPPFileName,'Delimiter','comma');
            
            out=GLMPP;
        end
        function out = update(GLMPP,FileName)          
            % read in .csv file
            GLMPPFileName=strcat(FileName);
            GLMPPArray=readcell(GLMPPFileName);
            
            % update property values of GLMPP
            for ii=1:size(GLMPPArray,1)
                string=GLMPPArray{ii,1};
                switch string
                    case 'NA'
                        GLMPP.NA=GLMPPArray{ii,2};
                    case 'M'
                        GLMPP.M=GLMPPArray{ii,2};
                    case 'ns'
                        GLMPP.ns=GLMPPArray{ii,2};
                    case 'ni'
                        GLMPP.ni=GLMPPArray{ii,2};
                    case 'DELTAt'
                        GLMPP.DELTAt=GLMPPArray{ii,2};
                    case 'zd'
                        GLMPP.zd=GLMPPArray{ii,2};
                    case 'zdstar'
                        GLMPP.zdstar=GLMPPArray{ii,2};
                    case 'np'
                        GLMPP.np=GLMPPArray{ii,2};
                    case 'na'
                        GLMPP.na=GLMPPArray{ii,2};
                    case 'beadsz'
                        GLMPP.beadsz=GLMPPArray{ii,2};
                    case 'K'
                        GLMPP.K=GLMPPArray{ii,2};
                    case 'PixelSize'
                        GLMPP.PixelSize=GLMPPArray{ii,2};
                end
                
                out = GLMPP;
            end
        end
        function out = write(GLMPP,FileName)
            GLMPPArray={'NA',GLMPP.NA;...
            'M',GLMPP.M;...
            'ns',GLMPP.ns;...
            'ni',GLMPP.ni;...
            'DELTAt',GLMPP.DELTAt;... % in um in image space
            'zd',GLMPP.zd;...
            'zdstar',GLMPP.zdstar;... % in um in image space
            'np',GLMPP.np;...
            'na',GLMPP.na;...
            'beadsz',GLMPP.beadsz;...
            'K',GLMPP.K;...
            'PixelSize',GLMPP.PixelSize}; % in um in object space: diameter of fluorescent bead
        
            % write cell array values to file
            GLMPPFileName=strcat(FileName);
            writecell(GLMPPArray,GLMPPFileName,'Delimiter','comma');
            
            out = GLMPP;
        end
    end
end
