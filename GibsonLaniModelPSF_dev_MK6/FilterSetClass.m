classdef FilterSetClass
    properties
        fields
        SpectralStructure
        mincutoff
        maxcutoff
        transmissionrange
        Wavelength
        numberlambda
        Filter
        specintegral
        filterwidthvec
        filteredges
        intstep
    end
    methods
        function out = init(FilterSet, SpectralStructurein,...
                filterwidthvecin, mincutoffin, maxcutoffin)
            FilterSet.fields = fieldnames(SpectralStructurein);
            FilterSet.SpectralStructure = SpectralStructurein;
            for ii = 1:numel(FilterSet.fields)
                cast(FilterSet.SpectralStructure.(FilterSet.fields{ii}).Wavelength,'single');
                cast(FilterSet.SpectralStructure.(FilterSet.fields{ii}).Spectrum,'single');
                cast(FilterSet.SpectralStructure.(FilterSet.fields{ii}).Filter,'single');
                FilterSet.SpectralStructure.(FilterSet.fields{ii}).Filter(:) =...
                    single(1);
            end
            FilterSet.mincutoff = single(mincutoffin);
            FilterSet.maxcutoff = single(maxcutoffin);
            FilterSet.transmissionrange = FilterSet.maxcutoff-...
                FilterSet.mincutoff;
            FilterSet.Wavelength = single(SpectralStructurein.(FilterSet.fields{1}).Wavelength);
            FilterSet.numberlambda = numel(FilterSet.Wavelength);
            %filterwidthvecin consists of numbers, the relative size of
            %which relates to the relative width of the spectral bins.
            %These widths are normalized such that their sum is one and
            %then the vector is multiplied by the transmission range so
            %that the elements of filterwidthvec give the width of the
            %spectral bins in microns
            FilterSet.filterwidthvec = single(FilterSet.transmissionrange.*...
                (filterwidthvecin./sum(filterwidthvecin))); 
            for ii = 1:numel(FilterSet.filterwidthvec)
                for jj = 1:numel(FilterSet.fields)
                    
                    FilterSet.SpectralStructure(ii).(FilterSet.fields{jj}).Filter = ...
                        single(padarray(SpectralStructurein.(FilterSet.fields{jj}).Wavelength,[200,0]));
                    FilterSet.SpectralStructure(ii).(FilterSet.fields{jj}).Wavelength = ...
                        cat(1,...
                        ...
                        single(linspace(...
                        min(SpectralStructurein.(FilterSet.fields{jj}).Wavelength)-0.200,...
                        min(SpectralStructurein.(FilterSet.fields{jj}).Wavelength)-.0001,...
                        200)'),...
                        ...
                        single(SpectralStructurein.(FilterSet.fields{jj}).Wavelength),...
                        ...
                        single(linspace(...
                        max(SpectralStructurein.(FilterSet.fields{jj}).Wavelength)+0.001,...
                        max(SpectralStructurein.(FilterSet.fields{jj}).Wavelength)+0.200,...
                        200)...
                            )');
                    FilterSet.SpectralStructure(ii).(FilterSet.fields{jj}).Spectrum = ...
                        single(padarray(SpectralStructurein.(FilterSet.fields{jj}).Spectrum,[200,0]));
                    FilterSet.SpectralStructure(ii).(FilterSet.fields{jj}).photonfraction = ... 
                        single(0);
                end
            end
            FilterSet.filteredges=single(zeros(2,numel(FilterSet.filterwidthvec)));
            for ii = 1:numel(FilterSet.filterwidthvec)
                if (ii==1)
                    FilterSet.filteredges(1,ii)=FilterSet.mincutoff;
                    FilterSet.filteredges(2,ii)=FilterSet.mincutoff+...
                        FilterSet.filterwidthvec(ii);
                else
                    FilterSet.filteredges(1,ii)=FilterSet.filteredges(2,ii-1);
                    FilterSet.filteredges(2,ii)=FilterSet.filteredges(2,ii-1)+FilterSet.filterwidthvec(ii);
                end
            end
            FilterSet.intstep=0.001;
            
            out = FilterSet;
        end
        function out = create(FilterSet)
            for ii = 1:numel(FilterSet.filterwidthvec)
                for jj = 1:numel(FilterSet.fields)
                    FilterSet = FilterSet.notch(ii);
                    FilterSet.SpectralStructure(ii).(FilterSet.fields{jj}).Filter = ...
                        FilterSet.Filter;
                    FilterSet.intstep = FilterSet.SpectralStructure(ii).(FilterSet.fields{jj}).Wavelength(2)-...
                    FilterSet.SpectralStructure(ii).(FilterSet.fields{jj}).Wavelength(1);
                    FilterSet.SpectralStructure(ii).(FilterSet.fields{jj}).photonfraction = ...
                        (FilterSet.intstep*trapz(...
                            FilterSet.SpectralStructure(ii).(FilterSet.fields{jj}).Filter.*...
                            FilterSet.SpectralStructure(ii).(FilterSet.fields{jj}).Spectrum...
                        ))/...
                        (FilterSet.intstep*trapz(...
                            FilterSet.SpectralStructure(ii).(FilterSet.fields{jj}).Spectrum...
                        ));
                end
            end
            
            out = FilterSet;
        end
        function out = notch(FilterSet,ii)
            FilterSet.Filter = heaviside((FilterSet.SpectralStructure(1).(FilterSet.fields{1}).Wavelength-FilterSet.filteredges(1,ii))).*...
                heaviside(-(FilterSet.SpectralStructure(1).(FilterSet.fields{1}).Wavelength-FilterSet.filteredges(2,ii)));
            out = FilterSet;
        end
    end
end