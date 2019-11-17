classdef SegmentDataClass
    properties
        frame
        sizeframe
        oi
        vari
        gi
        ROIstruct
        Kfull
        PSF_half_width_measure
        half_width_multiplier
        q
        q_multiplier
        photonthreshold
    end
    methods
        function out = init(SegmentData,GLMPP,framein,oiin,variin,giin)         
            SegmentData.oi=single(oiin);
            SegmentData.vari=single(variin);
            SegmentData.gi=single(giin);
            SegmentData.frame=single(framein);
            SegmentData.sizeframe=size(SegmentData.frame);
            SegmentData.ROIstruct=struct(...
                'ROI',single.empty,...
                'BACKGROUND',single.empty,...
                'ROIpositon',single.empty,...
                'oi',single.empty,...
                'vari',single.empty,...
                'gi',single.empty,...
                'ftvl0',single.empty,...
                'ftvl1',single.empty,...
                'Nphoton',single.empty);
            %SegmentData.ROIstruct=struct(...
            %    'ROI',single.empty,...
            %    'BACKGROUND',single.empty,...
            %    'ROIpositon',single.empty,...
            %    'oi',single.empty,...
            %    'vari',single.empty,...
            %    'gi',single.empty,...
            %    'ftvl0',single.empty,...
            %    'ftvl1',single.empty,...
            %    'Nphoton',single.empty,...
            %    'LLR',single.empty,...
            %    'Chi2',single.empty,...
            %    'varmin',single.empty);
            SegmentData.Kfull=single(GLMPP.K+2);
            SegmentData.PSF_half_width_measure=single(2);
            SegmentData.half_width_multiplier=single(2);
            SegmentData.q=SegmentData.half_width_multiplier*...
                SegmentData.PSF_half_width_measure;
            SegmentData.q_multiplier=single(2);
            SegmentData.photonthreshold=700;
            
            out = SegmentData;
        end
        function out = selectROIs(SegmentData,GLMPP)
            %PSF_half_width_measure=2;
            %half_width_multiplier=6;
            %q=half_width_multiplier*PSF_half_width_measure+1;
            %q_multiplier=2;

            A1=uniformfilter(SegmentData.frame,SegmentData.q_multiplier*SegmentData.q)-...
                uniformfilter(SegmentData.frame,SegmentData.q);
            %display(A1)
            %figure('Name','A1');
            %imagesc(A1),axis image
            A2=max_filter(A1,5*SegmentData.PSF_half_width_measure);
            %figure('Name','A2');
            %imagesc(A2), axis image;
            A3=(A1==A2);
            A3=A3(SegmentData.q_multiplier*SegmentData.q:end-SegmentData.q_multiplier*SegmentData.q,...
                  SegmentData.q_multiplier*SegmentData.q:end-SegmentData.q_multiplier*SegmentData.q);
            %figure('Name','A3');
            %imagesc(A3), axis image; 
            [I,J]=ind2sub(size(A3),find(A3));
            I=I+SegmentData.q_multiplier*SegmentData.q;
            J=J+SegmentData.q_multiplier*SegmentData.q;
            %figure('Name','overlay');
            %imagesc(SegmentData.frame);
            %axis image, colormap(gray),...
            %hold on, scatter(J,I); hold off;
            Icorner=I-GLMPP.K/2;
            Jcorner=J-GLMPP.K/2;
            Iwidth=GLMPP.K*ones(size(Icorner));
            Jwidth=GLMPP.K*ones(size(Icorner));

            ROIposition=single([Jcorner,Icorner,Jwidth,Iwidth]);
            warning('off','all');

            count = 0;
            for ii=1:size(ROIposition,1)                
                if ((ROIposition(ii,1)+ROIposition(ii,3)+1)<(SegmentData.sizeframe(2))) && ...
                        ((ROIposition(ii,2)+ROIposition(ii,4)+1)<(SegmentData.sizeframe(1))) &&...
                        (ROIposition(ii,1)>1) && (ROIposition(ii,2)>1) && (~isempty(ROIposition))                  
                                        
                    fullROI=SegmentData.frame(ROIposition(ii,2):ROIposition(ii,2)+ROIposition(ii,4)+1,...
                        ROIposition(ii,1):ROIposition(ii,1)+ROIposition(ii,3)+1);
                    fulloi=    SegmentData.oi(ROIposition(ii,2):ROIposition(ii,2)+ROIposition(ii,4)+1,...
                        ROIposition(ii,1):ROIposition(ii,1)+ROIposition(ii,3)+1);
                    fullvari=SegmentData.vari(ROIposition(ii,2):ROIposition(ii,2)+ROIposition(ii,4)+1,...
                        ROIposition(ii,1):ROIposition(ii,1)+ROIposition(ii,3)+1);
                    fullgi=    SegmentData.gi(ROIposition(ii,2):ROIposition(ii,2)+ROIposition(ii,4)+1,...
                        ROIposition(ii,1):ROIposition(ii,1)+ROIposition(ii,3)+1);

                    temp=processROI(fullROI,fulloi,fullvari,fullgi);
                     
                    Nphoton=sum(sum((temp.ROI-temp.oi)./temp.gi-temp.BACKGROUND));
                    
                    r1r2varval=r1r2vartest(...
                        temp.ROI,temp.gi,temp.oi,temp.BACKGROUND);

                    if (r1r2varval==1) && Nphoton>SegmentData.photonthreshold && (~isempty(fullROI))
                        
                        count = count + 1;
                        
                        SegmentData.ROIstruct(count).ROI=temp.ROI;
                        SegmentData.ROIstruct(count).BACKGROUND=temp.BACKGROUND;
                        SegmentData.ROIstruct(count).ROIposition=ROIposition(ii,:);
                        SegmentData.ROIstruct(count).oi=temp.oi;
                        SegmentData.ROIstruct(count).vari=temp.vari;
                        SegmentData.ROIstruct(count).gi=temp.gi;
                        SegmentData.ROIstruct(count).Nphoton=Nphoton;
                    end
                end
            end
            out = SegmentData;
        end
        function out = drawROIs(SegmentData)
            imagesc(SegmentData.frame), axis image;
            colormap(gray), colorbar;
            hold on
            for ii=1:numel(SegmentData.ROIstruct)
                if ~isempty(SegmentData.ROIstruct(ii).ROIposition)
                    rectangle('Position',...
                        [SegmentData.ROIstruct(ii).ROIposition(1),...
                         SegmentData.ROIstruct(ii).ROIposition(2),...
                         SegmentData.ROIstruct(ii).ROIposition(3),...
                         SegmentData.ROIstruct(ii).ROIposition(4)],...
                              'EdgeColor','r','LineWidth',1);
                end
            end
            hold off
            out = SegmentData;
        end
        function out = plotROIs(SegmentData)
            figure('Name','ROIs');
            for ii = 1:numel(SegmentData.ROIstruct)
                subplot(ceil(numel(SegmentData.ROIstruct)/4),4,ii);
                imagesc(SegmentData.ROIstruct(ii).ROI), axis image;
                colormap(winter);
            end
            
            out = SegmentData;
        end
    end
end

%this function performs a uniform mean filtering on a frame
function out = uniformfilter(frame,q)
kernel_size=q;%kernel_size=6*PSF_half_width_measure+1; %PSF_half_width_measure is an integer
filtered_frame=zeros(size(frame));
if mod(q,2)==1 %if odd
    for ii=1:size(filtered_frame,1)
        for jj=1:size(filtered_frame,2)
            if (ii<(kernel_size-1)/2+1) || (ii>(size(frame,1)-(kernel_size-1)/2)) || (jj<(kernel_size-1)/2+1) || (jj>(size(frame,2)-(kernel_size-1)/2)) 
                filtered_frame(ii,jj)=0;
            else
                ii_temp=ii-(kernel_size-1)/2:ii+(kernel_size-1)/2;
                jj_temp=jj-(kernel_size-1)/2:jj+(kernel_size-1)/2;
                cropped_region=frame(ii_temp,jj_temp);
                filtered_frame(ii,jj)=sum(sum((cropped_region)));
            end
        end
    end
else
    for ii=1:size(filtered_frame,1)
        for jj=1:size(filtered_frame,2)
            if (ii<(kernel_size/2)) || (ii>(size(frame,1)-((kernel_size/2)+1))) || (jj<(kernel_size/2)) || (jj>(size(frame,2)-((kernel_size/2)+1)))
                filtered_frame(ii,jj)=0;
            else
                ii_temp=ii-((kernel_size/2)-1):ii+((kernel_size/2)+1);
                jj_temp=jj-((kernel_size/2)-1):jj+((kernel_size/2)+1);
                cropped_region=frame(ii_temp,jj_temp);
                filtered_frame(ii,jj)=sum(sum((cropped_region)));
            end
        end
    end
end
out = filtered_frame;
end

%this function performs a uniform mean filtering on a frame
function ret=max_filter(frame,q)
kernel_size=q;%kernel_size=6*PSF_half_width_measure+1; %PSF_half_width_measure is an integer
filtered_frame=zeros(size(frame));
if mod(q,2)==1 %if odd
    for ii=1:size(filtered_frame,1)
        for jj=1:size(filtered_frame,2)
            if (ii<(kernel_size-1)/2+1) || (ii>(size(frame,1)-(kernel_size-1)/2)) || (jj<(kernel_size-1)/2+1) || (jj>(size(frame,2)-(kernel_size-1)/2)) 
                filtered_frame(ii,jj)=0;
            else
                ii_temp=ii-(kernel_size-1)/2:ii+(kernel_size-1)/2;
                jj_temp=jj-(kernel_size-1)/2:jj+(kernel_size-1)/2;
                cropped_region=frame(ii_temp,jj_temp);
                filtered_frame(ii,jj)=max(max(cropped_region));
            end
        end
    end
else
    for ii=1:size(filtered_frame,1)
        for jj=1:size(filtered_frame,2)
            if (ii<(kernel_size/2)) || (ii>(size(frame,1)-((kernel_size/2)+1))) || (jj<(kernel_size/2)) || (jj>(size(frame,2)-((kernel_size/2)+1)))
                filtered_frame(ii,jj)=0;
            else
                ii_temp=ii-((kernel_size/2)-1):ii+((kernel_size/2)+1);
                jj_temp=jj-((kernel_size/2)-1):jj+((kernel_size/2)+1);
                cropped_region=frame(ii_temp,jj_temp);
                filtered_frame(ii,jj)=max(max(cropped_region));
            end
        end
    end
end
ret=filtered_frame;
end

%this function estimates the background in photon counts inside an ROI 
%by calculating the averate of the pixels directly bordering an ROI.
function out = processROI(fullROI,fulloi,fullvari,fullgi)
ROI=(fullROI-fulloi)./fullgi;
nelem=numel(fullROI)-numel(fullROI(2:end-1,2:end-1));
out = struct('ROI',fullROI(2:end-1,2:end-1),...
    'BACKGROUND',...
    abs((1/nelem)*...
    (sum(ROI(1,1:end))+...
     sum(ROI(end,1:end))+...
     sum(ROI(2:end-1,1))+...
     sum(ROI(2:end-1,end)))...
     *ones(size(fullROI,1)-2)),...
    'oi',abs(fulloi(2:end-1,2:end-1)),...
    'vari',abs(fullvari(2:end-1,2:end-1)),...
    'gi',abs(fullgi(2:end-1,2:end-1)));
end

%this function calculates the variance of an inner region, 1, of equal area 
%to an outer region, 2, see below:
%
%           2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
%           2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
%           2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
%           2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2
%           2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2
%           2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2
%           2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2
%           2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2
%           2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2
%           2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2
%           2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2
%           2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2
%           2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2
%           2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2
%           2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2
%           2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
%           2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
%           2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
%
%the outer region has a width, j, determined by the following formulae:
% 
%            j=K/(2*(sqrt(2)+2))=K/(2*(sqrt(2)+2)).
%
%The algorithm relies on the fact that variances are additive. The variance
%of region 1 and region 1 should add to the variance of the whole ROI,
%regardless of if a signal is present or not. If there is no signal then,
%because the regions have approximately the same number of elements,
%two times the variance region 1 should approximately equal the variance
%of the whole ROI. Therefor, with the scalar alpha, the condition that the 
%variance of region one is greater than alpha/2 times the variance of the
%whole region of interest is equivalent to the presence of a signal.
function out = r1r2vartest(ROI,gi,oi,BACKGROUND)
DATA=(ROI-oi)./gi-BACKGROUND;
K=size(DATA,1);
j=round(K/(2*(sqrt(2)+2)));
ROIvar=var(var(DATA));
r1var=var(var(DATA(j+1:K-j,j+1:K-j)));

alpha=1.2;
if r1var>(alpha/2)*ROIvar
    out = 1;
else
    out = 0;
end
end

%this function calculates whether the background subtracted photon count is less than or greater than a photon threshold
%if the photon number is less than the the threshold the function returns 0
%other wise if the photon number is greater than the threshold it returns 1
function out = photonthresholdtest(ROI,gi,oi,BACKGROUND,photonthreshold)
Nphoton=sum(sum(...
    (ROI-oi)./gi-BACKGROUND...
    ));

if numel(photonthreshold)==1
    if Nphoton>=photonthreshold
        out = 1;
    else
        out = 0;
    end
else
    if (Nphoton>photonthreshold(1)) && (Nphoton<photothreshold(2))
        out = 1;
    else
        out = 0;
    end
end
end