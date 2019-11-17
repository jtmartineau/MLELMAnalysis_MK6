% this function caculates mu and the derivateves of mu for one PSF composed
% of five wavelengths sampled from a spectrum
function ret=mu_nPointSpectrum(K_X,K_Y,PIF,N_photon,NA,n,delta_lambda,Wavelength,Spectrum,Filter,x,x0,M,y,y0,rhotheta,Bstar,nw,z0,np,max_stroke,THETA_Cools)
Spectral=nPointSpectrum(n,delta_lambda,Wavelength,Spectrum,Filter);
mu_vec=zeros(K_Y,K_X,numel(Spectral.Wavelength));
for ii=1:numel(Spectral.Wavelength)
    lambda=Spectral.Wavelength(ii);
    mu_vec(:,:,ii)=Spectral.Transmission(ii)*mu_Cools_mex(K_X,K_Y,PIF,NA,N_photon,lambda,x,x0,M,y,y0,rhotheta,Bstar,nw,z0,np,max_stroke,THETA_Cools);
end                                          
mu=sum(mu_vec,3);
A=1/(sum(sum(mu)));
ret=N_photon*(A*mu);
end