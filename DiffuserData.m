%DiffuserData
%chemical properties
mu=1.002e-3;
rho=999.9;
roughness=0;
g=9.82;%gravitational acceleration
%The diameter at volute outlet is 1000mm
%Assuming constant velocity through the volute:
%Case specific properties
Q=10;
c2u=sum(cu(1,:))/jEllipse;%mean tangential velocity at outlet
c2m=sum(cm(1,:))/jEllipse;%meridional velocity at outlet
u2=mean(ufun(R(1,:)));%runner velocity at outlet
alpha=deg2rad(0);


%Dimensions
LGuideVane=0.4;%Length of first part of guide vane
LGuideVane2=0.1;%Length of second part of guide vane
LGuideVane1=LGuideVane-LGuideVane2;%Length of guide vane
zle=17;%number of guide vanes
DOutlet=2;%diameter of the volute at the inlet to the diffuser
D2Outlet=3;%diameter of the diffuser at the outlet to the penstock
LDiffuser=5;%length of the diffuser
RGuideVanes=1;%location of the guide vane pin(?)
SpiralCasingRadii=[3.  , 3.15, 3.3 , 3.45, 3.6 , 3.75, 3.9 , 4.05, 4.2 , 4.35, 4.5 ,4.65, 4.8 , 4.95, 5.1 , 5.25, 5.4 , 5.55, 5.7 , 5.85, 6.];%List of Radii of the outer volute wall in m
SpiralCasingDiameter=SpiralCasingRadii-RGuideVanes;%Diameter of the volute
OuterSpiralCasingLengths=[ 900.        ,  921.05263158,  942.10526316,  963.15789474, 984.21052632, 1005.26315789, 1026.31578947, 1047.36842105,1068.42105263, 1089.47368421, 1110.52631579, 1131.57894737,1152.63157895, 1173.68421053, 1194.73684211, 1215.78947368,1236.84210526, 1257.89473684, 1278.94736842, 1300.        ];%list with length of each outer volute element

b3=1;%height guide vane inlet
b4=1;%height guide vane outlet

%Calculated dimenions
% assume the guide vanes to be flat
r3=sqrt((RGuideVanes-cos(alpha)*LGuideVane1)^2+(sin(alpha)*LGuideVane1)^2);
r4=sqrt((RGuideVanes+cos(alpha)*LGuideVane2)^2+(sin(alpha)*LGuideVane2)^2);


a3=2*r3*sin(2*pi/zle/2);%width guide vane inlet
a4=2*r4*sin(2*pi/zle/2);%width guide vane outlet

c_volute=4*Q/(pi*DOutlet^2);%assuming constant velocity in the entire volute

B=(SpiralCasingDiameter(1:end-1)./2+RGuideVanes);%Length from turbine centre to volute centre
C=(SpiralCasingDiameter(2:end)./2+RGuideVanes);
a=OuterSpiralCasingLengths;%Renaming to make expression easier
b=SpiralCasingRadii(1:end-1);%same
c=SpiralCasingRadii(2:end);%same
A=sqrt(B.^2+C.^2+B.*C./b./c.*(a.^2-b.^2-c.^2));%finding length at the centre of volute, using law of cosineson two triangles with one shared angle

meanSpiralCasingLengths=A;%Renaming to make sense
meanSpiralCasingDiameters=(SpiralCasingDiameter(1:end-1)+SpiralCasingDiameter(2:end))/2;%finding mean diameter of each element

%Losses
HDiffuser=DiffuserLossesPump(c2u,c2m,Q,a3,b3,a4,b4,LGuideVane,mean(u2),zle,g);%calculating losses in the guide vanes
% HStayvanes=%Should calculate losses in stay vanes as well.
%%
AR=D2Outlet^2/DOutlet^2;%AR of volute outlet
fprintf('\nGülich, figure 1.19, L/R1= %.2f, AR-1=%.2f\n',2*LDiffuser/DOutlet,AR-1)
cp=0.5;%cp of volute outlet
cx=c_volute;
HVolute=voluteLossesPump(roughness,meanSpiralCasingLengths,c_volute,mu,rho,meanSpiralCasingDiameters,Q,u2,cx,cp,AR);%calculating losses in the volute
