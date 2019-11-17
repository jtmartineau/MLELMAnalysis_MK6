function ret=nPointSpectrum(n,delta_lambda,Wavelength,Spectrum,Filter)
Transmission=Spectrum.*Filter;

% calulate the index of the element in the transmitted spectrum with a
% value closest to the minimum lambda value. This block has to occur before
% the next one so the delta lambda is scalled correctly (in microns) to
% find the minimum lambda
SCOM_value=SCOM(Transmission,Wavelength);
min_lambda=SCOM_value-((n-1)/2)*delta_lambda;
min_lambda_vec=min_lambda.*ones(numel(Wavelength),1);
min_lambda_diff=abs(min_lambda_vec-Wavelength);
min_lambda_idx=find(min_lambda_diff==min(min_lambda_diff));


% spectral increments are in single nanometeres. delta_lambda is input in microns, but a value of delta_lambda=1 corresponds to a step of one nanometer. delta lambda is converted to nm and rounded to the nearest ingeral nanometer.
delta_lambda=round(1000*delta_lambda);

% select the n points
nPointWavelength=zeros(n,1);
nPointTransmission=zeros(n,1);
for ii=1:n
    nPointWavelength(ii)=Wavelength(min_lambda_idx+(ii-1)*delta_lambda);
    nPointTransmission(ii)=Transmission(min_lambda_idx+(ii-1)*delta_lambda);
end
ret=struct('Wavelength',nPointWavelength,...
           'Transmission',nPointTransmission);
end