function deltaXoutlet=getDeltaX2(D2,Dhub,jEllipse)
%This function is supposed to output a vector for x-values (or R-values for
%dividing a circular outlet into jEllipse-1 parts of equal area.

deltaAoutlet=pi*((D2./2).^2-(Dhub./2).^2)./(jEllipse-1);       %outlet area between two streamlines
deltaXoutlet=zeros(jEllipse,1);
deltaXoutlet(1)=Dhub/2;
deltaXoutlet(end)=D2/2;
for i=2:1:jEllipse-1
   deltaXoutlet(i)=sqrt(deltaXoutlet(i-1).^2+deltaAoutlet./pi);
end
end
