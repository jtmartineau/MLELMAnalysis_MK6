function out = GibsonLaniPSFModel(x, y, x0, y0, zp, lambda,...
                            NA, M, ns, ni, DELTAt, zd, zdstar, np, na,...
                            rho, theta, Bstar, THETAplate, THETAZernike)
%#codegen
                        
% NA = numerical aperture
% M = Magnification
% x = horizontal (column-wise) coordinate in image plane
% x0 = center of PSF coordinate system in object space
% y = vertical (row-wise) coordinate in object space
% rho = radial polar coordinate in the pupil plane
% theta = angular polar coordinate in pupil plane
% lambda = wavlength
% ns = index of refraction of sample
% zp = z position in the sample plane of the emitter, object space
% ni = index of refcation of the cover slip
% DELTAt = distance between focal position and zero axial coordinate in 
%       the object plane
% zd = focal offset of the detector from the focal plane of the lens
%       system
% zdstar = the tube length of the imaging system before the detector
% np = the index of refraction of the phase plate
% THETAPlate = the actual topography of the phase plate in um
% ThetaZernike = the non-phase wrapped phase function of the super-position
%       of Zernike polynomial aberrations
% na = index of refraction of air

% this expression was converted to matlab code form mathematica code with
% the "ToMatlab" package in mathematica.


integrand = exp(1).^((sqrt(complex(-1))*2).*lambda.^(-1).*pi.*((-1).*DELTAt.*ni.*(complex(1+(...
-1).*NA.^2.*ni.^(-2).*rho.^2)).^(1/2)+(-1).*np.*THETAplate+(-1).*...
 na.*THETAZernike+(-1/2).*rho.*zd.^(-1).*zdstar.^(-1).*((-1).*zd+...
zdstar)+(-1).*ns.*(complex(1+(-1).*NA.^2.*ns.^(-2).*rho.^2)).^(1/2).*zp+...
M.^(-1).*NA.*rho.*(complex((x+(-1).*x0).^2+(y+(-1).*y0).^2)).^(1/2).*cos(...
theta+(-1).*atan2(x+(-1).*x0,y+(-1).*y0)))).*lambda.^(-1).*M.^(-2)...
.*NA.^2.*rho;

integrand(isinf(integrand))=0;
integrand(isnan(integrand))=0;

integral=sum(integrand.*complex(Bstar),3);

out = integral.*conj(integral); 
end
