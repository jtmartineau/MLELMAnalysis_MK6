% this function calculates the spectral center of mass of the obvserved
% spectrum

function out = SCOM(Spectrum,Wavelengths)
if (max(Spectrum)==0)
    out = median(Spectrum);
else
    out = (sum(Spectrum.*Wavelengths))/(sum(Spectrum));
end
end