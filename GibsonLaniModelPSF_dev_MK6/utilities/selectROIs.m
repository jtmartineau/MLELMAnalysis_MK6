%this function identifies ROIs as well as estimates the background of each
%ROI and filters ROIs with a background-subtracted signal that is below
%photon_threshold.
function out = selectROIs(SegmentData,GLMPP)

%PSF_half_width_measure=2;
%half_width_multiplier=6;
%q=half_width_multiplier*PSF_half_width_measure+1;
%q_multiplier=2;

A1=uniformfilter(SegmentData.frame,SegmentData.oi,SegmentData.vari,SegmentData.gi,SegmentData.q)-...
    uniformfilter(SegmentData.frame,SegmentData.oi,SegmentData.vari,SegmentData.gi,2*SegmentData.q);
A2=max_filter(A1,5*SegmentData.PSF_half_width_measure);
A3=(A1==A2);
A3=A3(SegmentData.q_multiplier*SegmentData.q:end-SegmentData.q_multiplier*SegmentData.q,...
      SegmentData.q_multiplier*SegmentData.q:end-SegmentData.q_multiplier*SegmentData.q);
[I,J]=ind2sub(size(A3),find(A3));
Icorner=I-GLMPP.K/2+SegmentData.q_multiplier*q;
Jcorner=J-GLMPP.K/2+SegmentData.q_multiplier*q;
Iwidth=GLMPP.K*ones(size(Icorner));
Jwidth=GLMPP.K*ones(size(Jcorner));

ROIposition(:)=single([Jcorner,Icorner,Jwidth,Iwidth]);
warning('off','all');

for ii=1:size(ROIposition,1)
    
    fullROI=frame(ROIposition(ii,2)-1:ROIposition(ii,2)+ROIposition(ii,4),...
        ROIposition(ii,1)-1:ROIposition(ii,1)+ROIposition(ii,3));
    fulloi=    oi(ROIposition(ii,2)-1:ROIposition(ii,2)+ROIposition(ii,4),...
        ROIposition(ii,1)-1:ROIposition(ii,1)+ROIposition(ii,3));
    fullvari=vari(ROIposition(ii,2)-1:ROIposition(ii,2)+ROIposition(ii,4),...
        ROIposition(ii,1)-1:ROIposition(ii,1)+ROIposition(ii,3));
    fullgi=    gi(ROIposition(ii,2)-1:ROIposition(ii,2)+ROIposition(ii,4),...
        ROIposition(ii,1)-1:ROIposition(ii,1)+ROIposition(ii,3));
    
    temp=processROI(fullROI,fulloi,fullvari,fullgi,SegmentData.Kfull);
    
    r1r2varval=r1r2vartest(...
        temp.ROI,temp.gi,temp.oi,temp.BACKGROUND);
    photonthresholdval=photonthresholdtest(...
        ROI,gi,oi,BACKGROUND,photonthreshold);
    
    if (r1r2varval==1) && (photonthresholdval==1)
        SegmentData.ROIstruct(ii).ROI=temp.ROI;
        SegmentData.ROIstruct(ii).BACKGROUND=temp.BACKGROUND;
        SegmentData.ROIstruct(ii).ROIposition(ii,:)=ROIposition(1,:);
        SegmentData.ROIstruct(ii).oi=temp.oi;
        SegmentData.ROIstruct(ii).vari=temp.vari;
        SegmentData.ROIstruct(ii).gi=temp.gi;
        SegmentData.ROIstruct(ii).Nphoton=sum(sum(...
            (temp.ROI-temp.oi)./temp.gi-temp.BACKGROUND));
    else
        continue
    end
end
out = SegmentData;
end