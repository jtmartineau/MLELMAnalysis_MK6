function [xcom,ycom] = COM(ROI,PixelSize,M)
    K=size(ROI,2);
    [xx,yy]=meshgrid(...
        linspace(-((K*PixelSize)/M)/2,((K*PixelSize)/M)/2,K),...
        linspace(-((K*PixelSize)/M)/2,((K*PixelSize)/M)/2,K));
    xcom=sum(sum(xx.*abs(ROI)))/sum(sum(abs(ROI)));
    ycom=sum(sum(yy.*abs(ROI)))/sum(sum(abs(ROI)));
end