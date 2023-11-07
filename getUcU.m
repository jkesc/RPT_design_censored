% function y=getUcU(x,UcU1,UcU2,eta)
% %returns a distribution of reduced UcU for values of x between 0 and 1. x should be
% %the fraction of the total bladelength.
% %This version of the function is based on a cosine function, because I
% %don't know what else I should base it on. The real distribution should be
% %somewhat more steep in the middle, at least.
% baseline =[1 0.95 0.90 0.85 0.80  0.5 0.30 0.25 0.2 0.15 0.1 0.05 0];%[1 0.95 0.9 0.85 0.75 0.60 0.50 0.40 0.25 0.15 0.10 0.5 0];
% %[1 0.78 0.63 0.44 0.33 0.25 0.16 0.09 0.04 0.02 0]; is stolen from the
% %hydraulic compendium.
% fractionOfStreamline=linspace(0,1,length(baseline));
% yx=interparc(x,baseline,fractionOfStreamline);
% y=zeros(size(x));
% for i=1:1:size(yx)
%     y(i)=yx(i,1);
% end
% 
% %y=eta.* cos(x.*pi)./2+0.5; %I just made something up
% %Lets try with a Bezier curve.
% [xBez,yBez]=UcUBezier;
% y=spline(xBez,eta.*yBez,x);
% 
% end
function UcU= getUcU(FractionOfBladelength,UcU1,UcU2,iEllipse,jEllipse,P)
    UcU=zeros(iEllipse,jEllipse);
    for j=1:1:jEllipse
        UcU(:,j)=UcU2(j)+ppval(extractBezierByPoints(0,P),FractionOfBladelength(:,j)).*(UcU1(j)-UcU2(j));
    end
end

