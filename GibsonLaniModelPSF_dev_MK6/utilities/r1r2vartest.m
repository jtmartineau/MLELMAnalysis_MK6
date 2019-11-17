%this function calculates the variance of an inner region, 1, of equal area 
%to an outer region, 2, see below:
%
%           2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
%           2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
%           2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
%           2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2
%           2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2
%           2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2
%           2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2
%           2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2
%           2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2
%           2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2
%           2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2
%           2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2
%           2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2
%           2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2
%           2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2
%           2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
%           2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
%           2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
%
%the outer region has a width, j, determined by the following formulae:
% 
%            j=K/(2*(sqrt(2)+2))=K/(2*(sqrt(2)+2)).
%
%The algorithm relies on the fact that variances are additive. The variance
%of region 1 and region 1 should add to the variance of the whole ROI,
%regardless of if a signal is present or not. If there is no signal then,
%because the regions have approximately the same number of elements,
%two times the variance region 1 should approximately equal the variance
%of the whole ROI. Therefor, with the scalar alpha, the condition that the 
%variance of region one is greater than alpha/2 times the variance of the
%whole region of interest is equivalent to the presence of a signal.
function out = r1r2vartest(ROI,gi,oi,BACKGROUND)
DATA=(ROI-oi)./gi-BACKGROUND;
K=size(DATA,1);
j=round(K/(2*(sqrt(2)+2)));
ROIvar=var(var(DATA));
r1var=var(var(DATA(j+1:K-j,j+1:K-j)));

alpha=3;
if r1var>(alpha/2)*ROIvar
    out = 1;
else
    out = 0;
end
end