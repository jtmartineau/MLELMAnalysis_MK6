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