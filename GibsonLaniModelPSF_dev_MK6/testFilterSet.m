% this script tests FilterSetClass
clc
clear

GLMPP=GLMPPClass;
GLMPP=GLMPP.init('GLMPP_temp.csv');
GLMPP=GLMPP.update('/media/strgA/July2019LiveCellImaging_exported_tiffs/CAL_20190911040227/GLMPP_XPP_20190911040227.csv');
PupilPlane=PupilPlaneClass;
PupilPlane=PupilPlane.init('PupilPlane_temp.csv');
PupilPlane=PupilPlane.update('/media/strgA/July2019LiveCellImaging_exported_tiffs/CAL_20190911040227/PupilPlane_XPP_20190911040227.csv');
SpectralStructurein = load('/media/strgA/July2019LiveCellImaging_exported_tiffs/SpectralStructure_SEMPNTB_AXF0647_AXF0660_AXF0680_AXF0700.mat');
filterwidthvecin=ones(1,10);
mincutoffin=0.650;
maxcutoffin=0.850;
Nphoton=1000;
x0in=0;
y0in=0;
zpin=0;

GibsonLaniPSF=GibsonLaniPSFClass;
GibsonLaniPSF=GibsonLaniPSF.init(GLMPP,PupilPlane);
GibsonLaniPSF=...
    GibsonLaniPSF.fitlersetstandard(filterwidthvecin,...
                    mincutoffin, maxcutoffin, x0in, y0in, zpin, ...
                    SpectralStructurein);

for ii=1:size(GibsonLaniPSF.filtersetpsfarray,3)
   subplot(size(GibsonLaniPSF.filtersetpsfarray,3),1,ii)
   imagesc(poissrnd(Nphoton*GibsonLaniPSF.filtersetpsfarray(:,:,ii)));
   axis image
   colormap gray
   title(sprintf('Nphoton = %0.3g',...
       sum(sum(...
       poissrnd(Nphoton*GibsonLaniPSF.filtersetpsfarray(:,:,ii))...
       ))));
end

FilterSet = FilterSetClass;
FilterSet = FilterSet.init(SpectralStructurein,...
                filterwidthvecin, mincutoffin, maxcutoffin);
FilterSet = FilterSet.create;

figure(1)
for ii = 1:numel(FilterSet.SpectralStructure)
    for jj = 1:numel(FilterSet.fields)
        subplot(3,ceil(numel(FilterSet.filterwidthvec)/3),ii)
        hold on
        plot(FilterSet.SpectralStructure(ii).(FilterSet.fields{jj}).Wavelength,...
            FilterSet.SpectralStructure(ii).(FilterSet.fields{jj}).Spectrum);
        plot(FilterSet.SpectralStructure(ii).(FilterSet.fields{jj}).Wavelength,...
            FilterSet.SpectralStructure(ii).(FilterSet.fields{jj}).Filter);
        title(sprintf('Filter %d',ii))
        xlim([0.6,0.9])
        hold off
    end
end

fnames=fieldnames(SpectralStructurein);

figure(2)
for jj = 1:numel(FilterSet.fields)
    psfarray=[];
    for ii = 1:numel(FilterSet.SpectralStructure)
        wavelengthsin=FilterSet.SpectralStructure(ii).(FilterSet.fields{jj}).Wavelength;
        spectrain=FilterSet.SpectralStructure(ii).(FilterSet.fields{jj}).Spectrum;
        filtersin=FilterSet.SpectralStructure(ii).(FilterSet.fields{jj}).Filter;
        photonfraction=FilterSet.SpectralStructure(ii).(FilterSet.fields{jj}).photonfraction;
        GibsonLaniPSF=GibsonLaniPSF.spectralintensstandard(x0in, y0in, zpin, wavelengthsin, spectrain, filtersin);
        psfarray=cat(2, psfarray, Nphoton*photonfraction*GibsonLaniPSF.spectralPSF); 
    end
    %psfarray(isnan(psfarray))=0;    
    subplot(numel(FilterSet.fields),1,jj)
    imagesc(poissrnd(psfarray+20*ones(size(psfarray)))), axis image, colormap gray; colorbar;
    %title(sprintf('sum of PSF array = %0.3g',sum(sum(psfarray))));
    title(sprintf(fnames{jj}));
end