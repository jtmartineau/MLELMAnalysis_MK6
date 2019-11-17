%this function performs a uniform mean filtering on a frame
function out = uniformfilter(frame,oi,vari,gi,q)
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
                cropped_oi=oi(ii_temp,jj_temp);
                cropped_vari=vari(ii_temp,jj_temp);
                cropped_gi=gi(ii_temp,jj_temp);
                filtered_frame(ii,jj)=sum(sum(((cropped_region-cropped_oi)./(cropped_gi.*cropped_vari))))/sum(sum(1./cropped_vari));
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
                cropped_oi=oi(ii_temp,jj_temp);
                cropped_vari=vari(ii_temp,jj_temp);
                cropped_gi=gi(ii_temp,jj_temp);
                filtered_frame(ii,jj)=sum(sum(((cropped_region-cropped_oi)./(cropped_gi.*cropped_vari))))/sum(sum(1./cropped_vari));
            end
        end
    end
end
out = filtered_frame;
end