%Common problems: 
%The ellipse on TE and/or LE is sometimes problematic. may need to remove
%some points at start and end of the curve, not to get blade trimming error
%or something. Can be done by changing curveEnd and curveStart.

%Currently, we design the blade assuming no losses in the QH-curve because
%there are large uncertainties coupled to both blade blockage, slip and
%efficiency.

%Lower values of beta result in longer blades.
%TO DO:

%High priority:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Change gamma in such a way that the slip angle corresponds with the
%numerically found slip angle.

% Increase the outlet blade angle or find a way to reduce slip.

%Try to define the betaAlongBlade by two different variations, so the blade
%can have a decent ucu variation at high span, where it needs it the most.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Medium priority
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Currently assume that the booster pump delivers  m of head at BEP.
%Should be changed at off BEP calculations.

%PT and Helene said that the head increases by appx. m at Q= m^3/s. This
%does not seem to fit with the measured values at all? what is going on?

%Write a file with all the actual values given, so that this may be
%redacted from the script to be added in the appendices.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Low priority
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make the plots pretty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%TurbinDesign2
clear
close all
clc
figno=1;
% %% Calculating approximate required Head. This only works on the school PC.
% VS=load('N:\Documents\NTNU\V2020\TEP4900\SiraKvinaInfo\Vasstandsdata.mat');
% DeltaH=VS.Roskreppfjorden-VS.Oyarvatn;
% Day=VS.Day((DeltaH>0));
% Hour=VS.Hour((DeltaH>0));
% month=VS.month((DeltaH>0));
% Year=VS.year((DeltaH>0));
% DeltaH=DeltaH(DeltaH>0);
% Day=Day(DeltaH<100);
% Hour=Hour(DeltaH<100);
% month=month(DeltaH<)100;
% Year=Year(DeltaH<100);
% DeltaH=DeltaH(DeltaH<100);
% figure(figno)
% figname.TotalHead=figno;
% figno=figno+1;
% hold on
% plot(month,DeltaH,'o')
% avg=zeros(12,1);
% for i=1:1:12
%    avg(i)=mean(DeltaH(month==i));
% end
% plot(avg)
% Head=mean(avg); %=100 may be more realistic
% title({'Head difference between higher and lower reservoir throughout the year',['Mean value:  ',num2str(Head),' m']})
% xlabel('Month')
% ylabel('Head [m]')
% legend('all values','average per month')
% %% mean difference between oyarvatn and the turbine
% figure(figno)
% figname.HeadDownstream=figno;
% figno=figno+1;
% hold on
% H_turbine=100;
% DeltaH=VS.Oyarvatn-H_turbine;
% Day=VS.Day((abs(DeltaH)<20));
% Hour=VS.Hour((abs(DeltaH)<20));
% month=VS.month((abs(DeltaH)<20));
% Year=VS.year((abs(DeltaH)<20));
% DeltaH=DeltaH((abs(DeltaH)<20));
% Day=Day(abs(DeltaH)<20);
% Hour=Hour(abs(DeltaH)<20);
% month=month(abs(DeltaH)<20);
% Year=Year(abs(DeltaH)<20);
% DeltaH=DeltaH(DeltaH<20);
% figure(figname.HeadDownstream)
% hold on
% plot(month,DeltaH,'o')
% avg=zeros(12,1);
% for i=1:1:12
%    avg(i)=mean(DeltaH(month==i));
% end
% plot(avg)
% %Assuming there will be most pumping between mid july and early september,
% %so months 7, 8 and 9:
% %Head=mean(avg(7:9))=100
% %but it may be stupid to not allow for pumping all year around, so
% HeadDiffSubmergence=mean(avg);%Approximately 8m
% title({'Difference between Øyarvatn and the turbine',['Mean value:  ',num2str(HeadDiffSubmergence),' m']})
% ylabel('Head [m]')
% xlabel('Month')
%%
dir='AnsysFiles';
if~exist(dir,'dir')
    mkdir(dir)
end
global g rho mu nu patm pvapour Values
%Choices%%%%%%%%%%%%%%%%
Choices;
%%%%%%%%%%%%%%%%%%%%%%
%Physical quantities
g=9.82;
patm=101325;                                    %Atmospheric pressure in Pa
rho=999.9;                                        %density of water
mu=1.519e-3;
nu=mu/rho;
pvapour=872.1;                                  %vapourization pressure of water at T = 5 [deg]
%Geometry values
getValues;
%Arbitrarily chosen temperature: 5 degrees.
HeadDiffSubmergence=Values.HeadDiffSubmergence;      %The height difference between the average lower reservoir height and the RPT
R1GuideVanes=Values.R1GuideVanes;                              %Radius of the guideVanes                                         %g in norway
Head=Values.Head;                                  %The design head
Q=Values.Q;%;                                           %From efficiencies measured in 2005, Q_BEP is approximately 60. Assume this to be the case for the new one as well.
 D1=Values.D1;%                       %Current diameter at the guide vanes
%D1=2*R1GuideVanes./1.05;                       D1 \approxeq 1.05D_guidevanes
b=Values.b;
Dhub=Values.Dhub;%;                                     %(Quite arbitrary chosen) hub diameter
A=pi.*D1.*b;                                    %Area at 1, no blade blockage
%A1=A;
% D2=sqrt(4.*D1.*b/acceleration+Dhub.^2);        %Diameter imposed by
%assumed acceleration
D2=Values.D2;                                       %Diameter imposed by existing geometry
A2=pi*(D2.^2-Dhub.^2)./4;
acceleration=A/A2;
D2m=D2;                                         %THIS SHOULD BE CALCULATED BY sqrt(0.5*(D2min^2+D2max^2), Geometric diameter
e=0.3*D1;                                     %blade thickness [X,X]*D1 Should be in the higher range for high heads.
n= Values.n;                                         %current rotational velocity increased by one syncronous vel.
%b=A2./(pi.*D1);                                %b at 1

iEllipse=Values.iEllipse;
jEllipse=Values.jEllipse;
ufun=@(R) (2.*pi.*R.*n)./(60);
R2=getDeltaX2(D2,Dhub,jEllipse);
runnerheight=Values.runnerheight;%; %From draft tube to top of guide vane.
eta=getEta_h(Q,Q,n,Head);                                                   %Don't I change this in QH-calculations anyway? By use of equations from Gülich
Routlet(1)=Values.Routlet(1);                                                         %Values from the supplied figure.
Zoutlet(1)=Values.Zoutlet(1);
Routlet(2)=Values.Routlet(2);
Zoutlet(2)=Values.Zoutlet(2);
RLe1=Values.RLe1;                                                                  %The target values for The leading edge
ZLe1=Values.ZLe1;
RLe2=Values.RLe2;
ZLe2=Values.ZLe2;
%% starting point for streamlines 
aEllipse=R1GuideVanes;                                                       %The half length of the ellipsis in the x-direction
bEllipse=runnerheight;                                                      %The half length of the ellipsis in the y-direction
switch(Geometry)
    case 'Ellipse'
        xEllipse=linspace(0,-aEllipse);                                        %to -D1 etc. because are in the 2nd quadrant of ellipsis.
        yEllipse= bEllipse.*sqrt(1-(xEllipse./aEllipse).^2);
        xyEllipse=interparc(iEllipse,xEllipse,yEllipse,'spline');
        xEllipse=xyEllipse(:,1);
        xEllipse=R1GuideVanes+xEllipse;
        yEllipse=xyEllipse(:,2);
    case 'Bezier'
        P1=[D1/2,runnerheight];
        %axis([0 R1GuideVanes 0 runnerheight])
        P2=Values.P2;%ginput();%[0,runnerheight];
        P3=[0,0];
        bezierpoints=   [P1
                        P2
                        P3];
       BezierCurve=extractBezierByPoints(0,bezierpoints);
       BezierCurve=ppval(BezierCurve,linspace(0,R1GuideVanes));
%        BezierCurve=Bezier(0,-aEllipse,bEllipse,'Plotnt');                %Creating an arbitrary chosen bezier-curve to connect the outlet and inlet
       xyBezier=interparc(iEllipse,linspace(0,2),BezierCurve,'spline'); %Dividing it into iEllipse equally spaced parts
       xBezier=xyBezier(:,1);                                                   %Extracting data for the points as two vectors
       yBezier=xyBezier(:,2);
       xEllipse=xBezier(end:-1:1);                                              %and putting them in the same form as would have been used for the Ellipsis
       yEllipse=yBezier(end:-1:1);
    case 'BezierDefByPoints'
        P1=[D1/2,runnerheight];                                             %Defining start and stop of hub curve
        P3=[0,0];       
        distLim=0.001;                                                       %Criterium for when the curve is close enough to the required point at the inlet of the current geometry
        dR=0;                                                               %initial value for iterative value for R
        dZ=0.1*runnerheight;                                                %iterative change in Z
        rPoint=Routlet(1)./2;                                                   %Radius for the control point for the curve                                                    
        zPoint=runnerheight;                                                %z-value for the bezier control point
        if Recording
            OuterStreamlineObj=VideoWriter('outer.avi');
            OuterStreamlineObj.FrameRate=10;
            open(OuterStreamlineObj);
        end
        while true                                                          %do-while loop in matlab
            rPoint=rPoint+dR;                                               %change the radius
            P2=[rPoint,zPoint];                                             %Define the control point
            bezierpoints=   [P1                                             %vector with all the control points
                            P2
                            P3];
           BezierCurve=extractBezierByPoints(0,bezierpoints);
           BezierCurve=ppval(BezierCurve,linspace(0,R1GuideVanes));                 
           xyBezier=interparc(iEllipse,linspace(0,1.6),BezierCurve,'spline'); %Dividing it into iEllipse equally spaced parts
           xBezier=xyBezier(:,1);                                                   %Extracting data for the points as two vectors
           yBezier=xyBezier(:,2);
           xEllipse=xBezier(end:-1:1);                                              %and putting them in the same form as would have been used for the Ellipsis
           yEllipse=yBezier(end:-1:1);
           [pt,dist,t]=distance2curve([xEllipse,yEllipse],[Routlet(1),Zoutlet(1)],'pchip');%find the distance to the hub line and the point on the line closest to the required point
           if dist<distLim                                                  %If the distance is small enough, the solution has converged
               if yEllipse(2)>yEllipse(1)                                   %don't want the stream to flow upwards again.
                   zPoint=zPoint-dZ;                                        %If it does, reduce Z of the point
               else
               break                                                        %else the solution has converged. Exit the loop
               end
           end
           if pt(1)>Routlet(1)
               dR=-1*dist;
           else
               dR=1*dist;
           end
           figure(figno)
           hold off
           plot(xEllipse./Values.DimensionlessRadius,yEllipse./Values.DimensionlessRadius)
           hold on
           xlabel('R/R_1')
           ylabel('Z/R_1')
           title('Streamline along the hub')
           figname.hubStreamline=figno;
           
           
          
           plot([P1(1),P2(1),P3(1)]./Values.DimensionlessRadius,[P1(2),P2(2),P3(2)]./Values.DimensionlessRadius,'bo')
           plot(pt(1)./Values.DimensionlessRadius,pt(2)./Values.DimensionlessRadius,'rd')
           plot(Routlet(1)./Values.DimensionlessRadius,Zoutlet(1)./Values.DimensionlessRadius,'ks')
           line([P1(1),P2(1),P3(1)]./Values.DimensionlessRadius,[P1(2),P2(2),P3(2)]./Values.DimensionlessRadius)
           drawnow
           if Recording
           currFrame=getframe(gcf);
           writeVideo(OuterStreamlineObj,currFrame);
           end
        end
        if Recording
            close(OuterStreamlineObj)
        end
        figno=figno+1;
    otherwise
        error('ERROR: The variable "Geometry" has no matching value.')
end
%% Streamlines
% x2R=@(x) R1GuideVanes+x;                                                           %Translation from x-coordinate to R-coordinate.
deltaYinlet=b/(jEllipse-1);                                                 %inlet value spacing for all the other streamlines
deltaXoutlet=getDeltaX2(D2,Dhub,jEllipse);                                  %outlet value placing for all the other streamlines calculated with same annular area at outlet
R=zeros(iEllipse,jEllipse);                                                 %initializing matrix of R-values
Z=R;                                                                        %initializing matrix of Z-values
R(:,1)=xEllipse;                                                       %setting values along the ring
Z(:,1)=yEllipse;                                                            %setting values along the ring
Z(iEllipse,:)=0;                                                            %Z-values at the outlet
R(iEllipse,:)=deltaXoutlet(1:1:end);
Z(1,:)=Z(1,1)-[0:deltaYinlet:b];                                            %I get red lines, but it seems to work anyhow.
R(1,:)=D1./2;                                                               %R-values at the inlet
alphaij=@(iv,jv,Rv,Zv) abs(atan((Zv(iv-1,jv)-Zv(iv+1,jv))./((Rv(iv-1,jv)-Rv(iv+1,jv)))));              %Function for the angle alphaij
dA1=pi*D1*b/(jEllipse-1);
dA2=pi./4*(D2.^2-Dhub.^2)./(jEllipse-1);                                                            %The meridional area is constant for a RPT If acceleration is introduced, this would have to change be a function of the distance along the blade.
bij=@(iv,jv,Rv,Zv) abs((Rv(iv,jv)-Rv(iv,jv+1))/(sin(alphaij(iv,jv,Rv,Zv))));
figure(figno)
figname.meridionalProjection=figno;
figno=figno+1;
hold on
axis equal
GGeom=zeros(size(R));
%plot(R(1,:)./Values.DimensionlessRadius,Z(1,:)./Values.DimensionlessRadius,'o')
%plot(R(:,1)./Values.DimensionlessRadius,Z(:,1)./Values.DimensionlessRadius,'o')
switch Geometry
    case 'BezierDefByPoints'
        dist=0;                                                             %initial distance between the blades
        ax=Values.ax;                                                            %initial x value to define the acceleration along the blade
        ay=Values.ay;                                                             %same for y
        dx=Values.dx;                                                            %Value to change ax and ay with for each iteration
        dy=Values.dy;
        if Recording
            shroudVideoObject=VideoWriter('shroud.avi');
            shroudVideoObject.FrameRate=10;
            open(shroudVideoObject)
        end
        while true                                                          %do-while
            P.Aij=   [                                                       %Define Bezier control point which defines the change in area along the blade
                    ax    ay
                    ax    ay
                    ];                                                      %At the moment it is defined by two equal points, because the algorithm didn't converge with only one point.
               % plot(R(:,1)./Values.DimensionlessRadius,Z(:,1)./Values.DimensionlessRadius,'-bo')
            for j=1:1:(jEllipse-1)                                          %jEllipse-1 that's where we find the streamline nr jEllipse.
                 for i=2:1:iEllipse
                     GGeom(i,j)=GGeom(i-1,j)+sqrt((R(i-1,j)-R(i,j)).^2+(Z(i-1,j)-Z(i,j)).^2);
                 end
                 for i=2:1:(iEllipse-1)
                   % plot(R(i,j)./Values.DimensionlessRadius,Z(i,j)./Values.DimensionlessRadius,'-ro')
                    R(i,j+1)=sqrt(R(i,j)^2+((getdAij(GGeom(i,j),dA1,dA2,P.Aij,max(GGeom(:,j)))*sin(alphaij(i,j,R,Z))/pi)));
                    Z(i,j+1)=Z(i,j)-bij(i,j,R,Z)*cos(alphaij(i,j,R,Z));
                    
                  %plot([R(i-1,j)./Values.DimensionlessRadius,R(i+1,j)./Values.DimensionlessRadius],[Z(i-1,j)./Values.DimensionlessRadius,Z(i+1,j)./Values.DimensionlessRadius],'-ro')
                  %plot([R(i,j)./Values.DimensionlessRadius,R(i,j+1)./Values.DimensionlessRadius],[Z(i,j)./Values.DimensionlessRadius,Z(i,j+1)./Values.DimensionlessRadius],'-ro')

                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                RZInterp=interparc(iEllipse,R(:,j+1),Z(:,j+1),'Spline');
                R(:,j+1)=RZInterp(:,1);% To NOT distribute the points for every iteration, remove these three lines.
                Z(:,j+1)=RZInterp(:,2);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 %plot(R(:,j+1)./Values.DimensionlessRadius,Z(:,j+1)./Values.DimensionlessRadius,'-b') %debugging
            end
            ptLast=pt;                                                      %To check wether the line has changed the side of the point
            [pt,dist,t]=distance2curve([R(:,jEllipse),Z(:,jEllipse)],[Routlet(2),Zoutlet(2)],'pchip');%nearest point on the line, distance between the points and parametric variable
            if sign(ptLast(2)-Zoutlet(2))~=sign(pt(2)-Zoutlet(2))             %If the line switches which side of the point it is on
                dy=.5*dy;                                                   %Reduce the iterative increase by a factor of 0.5, to ensure convergence
            end
            if dist<distLim                                                 %If the line has converged sufficiently
                break                                                       %Exit the while loop
            elseif Z(1,jEllipse)<Z(2,jEllipse)                              %if the streamline goes upward at the outlet (1)
                ax=ax+dx;                                                   %increase the acceleration at the outlet (1)
            elseif pt(2)>Zoutlet(2)                                          %If the point is greater than the curve
                ay=ay-dy;                                                   %Need less steep curve for acceleration in the start
            else                                                            %If the point is below the curve
                ay=ay+dy;                                                   %Need more acceleration in the start, to contract curve
            end
            if max(((getdAij(GGeom(:,jEllipse-1),dA1,dA2,P.Aij,max(GGeom(:,jEllipse-1))))-dA1)./(dA2-dA1))>1 %if the area has a local maximum which is not at the inlet
                ay=ay-dy;                                                   %reduce the area
                ax=ax-2*dx;                                                 %change conditions and try again
            end
           figure(figno)
           subplot(2,1,1)
           hold off
           plot(R./Values.DimensionlessRadius,Z./Values.DimensionlessRadius)
           hold on
           xlabel('R/R_1')
           ylabel('Z/R_1')
           %plot([P1(1),P2(1),P3(1)]./Values.DimensionlessRadius,[P1(2),P2(2),P3(2)]./Values.DimensionlessRadius,'bo')
           plot(pt(1)./Values.DimensionlessRadius,pt(2)./Values.DimensionlessRadius,'rd')
           plot(Routlet(2)./Values.DimensionlessRadius,Zoutlet(2)./Values.DimensionlessRadius,'ks')
           subplot(2,1,2)
           hold off
           plot(GGeom(:,jEllipse-1)./max(GGeom(:,jEllipse-1)),((getdAij(GGeom(:,jEllipse-1),dA1,dA2,P.Aij,max(GGeom(:,jEllipse-1))))-dA1)./(dA2-dA1))
           hold on
           xlabel('G/G_{max}')
           ylabel('dA')
           
           plot(ax,ay,'dk')
           drawnow
           if Recording
               currFrame=getframe(gcf);
               writeVideo(shroudVideoObject,currFrame)
           end
        end
        if Recording
            close(shroudVideoObject)
        end
        figname.shroudStreamline=figno;
        figno=figno+1;
    otherwise
    P.Aij=Values.P.Aij;
    for j=1:1:(jEllipse-1)%jEllipse-1 that's where we find the streamline nr jEllipse.
        for i=2:1:iEllipse
         GGeom(i,j)=GGeom(i-1,j)+sqrt((R(i-1,j)-R(i,j)).^2+(Z(i-1,j)-Z(i,j)).^2);
        end
        for i=2:1:(iEllipse-1)
          %  plot(R(i,j)./Values.DimensionlessRadius,Z(i,j)./Values.DimensionlessRadius,'o')
            R(i,j+1)=sqrt(R(i,j)^2+((getdAij(GGeom(i,j),dA1,dA2,P.Aij,max(GGeom(:,j)))*sin(alphaij(i,j,R,Z))/pi)));
            Z(i,j+1)=Z(i,j)-bij(i,j,R,Z)*cos(alphaij(i,j,R,Z));

        %      plot([R(i-1,j),R(i+1,j)]./Values.DimensionlessRadius,[Z(i-1,j),Z(i+1,j)]./Values.DimensionlessRadius,'-bo')
        %      plot([R(i,j),R(i,j+1)]./Values.DimensionlessRadius,[Z(i,j),Z(i,j+1)]./Values.DimensionlessRadius,'-ro')

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        RZInterp=interparc(iEllipse,R(:,j+1),Z(:,j+1),'Spline');
        R(:,j+1)=RZInterp(:,1);% To NOT distribute the points for every iteration, remove these three lines.
        Z(:,j+1)=RZInterp(:,2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %plot(R(:,j+1)./Values.DimensionlessRadius,Z(:,j+1)./Values.DimensionlessRadius,'o-') %debugging
    end
end
dAij=getdAij(GGeom(i,j),dA1,dA2,P.Aij,max(GGeom(:,j)));

    GGeom(:,end)=GGeom(:,end-1)+sqrt((R(:,end)-R(:,end-1)).^2+(Z(:,end)-Z(:,end-1)).^2);
    figure(figname.meridionalProjection)
    hold on
    plot(R./Values.DimensionlessRadius,Z./Values.DimensionlessRadius)
    axis equal
    xlabel('R/R_1')
    ylabel('Z/R_1')
    title('Streamlines of the turbine in the R-Z plane')
    axis equal
% plot(P2(:,1)./Values.DimensionlessRadius,P2(:,2)./Values.DimensionlessRadius,'d')
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main dims and streamlines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% above
%% QH calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Angles below
%The entire design collapses if the trailing edge is too far down the
%meridional view, so don't do that.
%Cutting off the blades
%columns like this: [xvalue,yvalue;...]. decreasing x-values, increasing 
%y-values (unsure if the last part matters).
switch Geometry
    case 'BezierDefByPoints'
P.LEp=Values.P.LEp;%The points which generate the Leading edge by use of Bezier curves
P.TEpmid=Values.P.TEpmid;
if isempty(P.TEpmid)
    P.TEpmid=[Routlet(2)+(Routlet(1)-Routlet(2))/2 Zoutlet(2)+(Zoutlet(1)-Zoutlet(2))/2];
end

P.thetaTE1=atan((abs(Routlet(2)-P.TEpmid(1)))./((abs(Zoutlet(2)-P.TEpmid(2)))));
P.TEp1=[P.TEpmid(1)+1.05*(sqrt((Routlet(2)-P.TEpmid(1))^2+(Zoutlet(2)-P.TEpmid(2))^2))*sin(P.thetaTE1) ...
    P.TEpmid(2)-1.05*(sqrt((Routlet(2)-P.TEpmid(1))^2+(Zoutlet(2)-P.TEpmid(2))^2))*cos(P.thetaTE1)];

P.thetaTE2=atan((abs(Routlet(1)-P.TEpmid(1)))./((abs(Zoutlet(1)-P.TEpmid(2)))));
P.TEp2=[P.TEpmid(1)-1.05*(sqrt((Routlet(1)-P.TEpmid(1))^2+(Zoutlet(1)-P.TEpmid(2))^2))*sin(P.thetaTE2) ...
    P.TEpmid(2)+1.05*(sqrt((Routlet(1)-P.TEpmid(1))^2+(Zoutlet(1)-P.TEpmid(2))^2))*cos(P.thetaTE2)];

P.TEp=[ 
    P.TEp1
    P.TEpmid
    P.TEp2
    ];%The points which generate the Trailing edge by use of Bezier curves
    otherwise
        P.LEp=Values.P.LEp;%The points which generate the Leading edge by use of Bezier curves
        P.TEp=Values.P.TEp;
end
RGeom=R;
ZGeom=Z;
[R,Z]=cutOffBlade(RGeom,ZGeom,P.TEp,P.LEp,iEllipse,jEllipse);
figure(figname.meridionalProjection)
plot(R./Values.DimensionlessRadius,Z./Values.DimensionlessRadius,'o')
P.TERunner=Values.P.TERunner;
P.LERunner=Values.P.LERunner;
[RRunner,ZRunner]=cutOffBlade(RGeom,ZGeom,P.TERunner,P.LERunner,iEllipse,jEllipse);

%% G-streamline
G=zeros(size(Z)); %Initializing the G matrix
for i=2:iEllipse
    G(i,:)=G(i-1,:)+sqrt((R(i-1,:)-R(i,:)).^2+(Z(i-1,:)-Z(i,:)).^2);
end
DeltaG=G(2:end,:)-G(1:end-1,:);

%Now we want to find beta1 and beta2 for all values along LE and TE
%To find beta1 along TE, we assume free vortex theory: cu*r=const, where
%the reference value is found at D1. 
GToAdd=getGGeomAtG0(R(1,:),Z(1,:),RGeom,ZGeom);                             %Finds the value of GGeom where G of the blade itself is 0.
QHCalculations
switch FindBetaBlade
    case 'freeVortex'
        beta1=getBeta1B(gamma,cm1,tau_1,D1,n,beta1B,R);                             %finds values for beta1 responding to free vortex theory
        beta2=zeros(size(beta1));                                                   %Initializing vector
        for j=1:1:jEllipse
            beta2(j)=atan((Q./(jEllipse-1)./getdAij((G(end,j)+GToAdd(j)),dA1,dA2,P.Aij,GGeom(end,j)))./ufun(R(end,j)));%finding beta2 according to the assumed cm at the point in the runner.
        end
        % beta2=atan((4.*Q)./(ufun(R(iEllipse,:)).*pi.*(D2.^2-Dhub.^2)));           %values for beta2 if the blades would end at the end of the geometry.
        beta1=beta1+(beta1<0).*pi;                                                  %if the angle gets larger than 90 degrees, trigonometry leads to it being 180deg smaller than intended. therefore we add 180deg (or pi radians) to fix this.
    case 'streamlines'
        beta2=zeros(size(beta1));                                                   %Initializing vector
        for j=1:1:jEllipse
            beta2(j)=atan((Q./(jEllipse-1)./getdAij((G(end,j)+GToAdd(j)),dA1,dA2,P.Aij,GGeom(end,j)))./ufun(R(end,j)));%finding beta2 according to the assumed cm at the point in the runner.
        end
    otherwise
        error('ERROR, no matching value for variable FindBetaBlade')
end

%% Energy distribution and finding beta along the blade
switch bladeLoadingMethod
    case 'UcU'
        velAngular=n.*2.*pi./60;%is a scalar, angular velocity of the turbine.
        uAlongBlade=velAngular.*R;% get one value for each node, peripheral velocity of the runner
        FractionOfBladelength=G./G(end,:);%fraction of the streamline along the blade
        cmAlongBlade=zeros(iEllipse,jEllipse);
        for j=1:1:length(cmAlongBlade)
            cmAlongBlade(:,j)=Q./(jEllipse-1)./getdAij((G(:,j)+GToAdd(:,j)),dA1,dA2,P.Aij,GGeom(end,j));
        end
        UcU1=uAlongBlade(1,:).*(uAlongBlade(1,:)-(cmAlongBlade(1,:))./(tan(beta1)));
        UcU2=uAlongBlade(end,:).*(uAlongBlade(end,:)-(cmAlongBlade(end,:))./(tan(beta2)));%Should be 0 (or very close to 0.
        P.UcU=Values.P.UcU;
        UcUAlongBlade=getUcU(FractionOfBladelength,UcU1,UcU2,iEllipse,jEllipse,P.UcU); %Should return a distribution for UcU along the blade.
        cuAlongBlade=UcUAlongBlade./uAlongBlade;
        betaAlongBlade=atan(cmAlongBlade./(uAlongBlade-cuAlongBlade)); %finding beta by use of velocity triangles along the blade.
        %Will probably have to design the UcU distribution in such a way that the
        %different streamlines have different distributions. i.e. smaller values
        %close to the hub.
        betaAlongBlade=betaAlongBlade+(betaAlongBlade<=0).*pi;%so my theory is that the un-smoothness of the blade originates from atan switching to the wrong quadrant. adding pi to all negative values should solve this.
        betaAlongBladeW=betaAlongBlade;%To make it compatible with later updates, I think...
    case 'BetaDistribution'
        %This part is supposed to sketch a smooth curve for beta, based on
        %a guess.
        
        P.BetaDist=Values.P.BetaDist;
        P.Spanwise=Values.P.Spanwise;%To control how the blade changes in the spanwise direction
        switch betaEval
            case 'lineSegment' %To evaluate values of beta at line segment between two nodes, and thus get more accurate values of theta.
                                %The exception is at start and end, because
                                %beta1 and 2 MUST be right.
                        FractionOfBladelength=G./G(end,:);
                        betaAlongBladeW=getBetaDist(beta1,beta2,iEllipse,jEllipse,FractionOfBladelength,'Node',P.BetaDist,P.Spanwise);%So the distribution seems to be equally spaced, and not

                        FractionOfBladelength=(G(1:end-1,:)+DeltaG./2)./G(end,:);
                        betaAlongBlade=getBetaDist(beta1,beta2,iEllipse,jEllipse,FractionOfBladelength,betaEval,P.BetaDist,P.Spanwise);%to monitor relative velocity
            case 'Node'
                        FractionOfBladelength=G./G(end,:);

                        betaAlongBlade=getBetaDist(beta1,beta2,iEllipse,jEllipse,FractionOfBladelength,betaEval,P.BetaDist,P.Spanwise);%So the distribution seems to be equally spaced, and not according to the fraction of G along the blade? Also seems like 'Bezier' is worthless? can change to Spline which does absolutely nothing. Who wrote this code? (I fear it was me :'-( )
                        betaAlongBladeW=betaAlongBlade;%to monitor relative velocity, w, later on.
            otherwise
                error('ERROR: no fitting value for betaEval')
        end
        
    otherwise
        error('ERROR: The variable "bladeLoadingMethod" has no matching value.')
end
%%Finding the velocities along the blade
cm=zeros(iEllipse,jEllipse);
for j=1:1:jEllipse
    cm(:,j)=(Q./(jEllipse-1)./getdAij((G(:,j)+GToAdd(:,j)),dA1,dA2,P.Aij,GGeom(end,j)));
end
w=sqrt(cm.^2+(cm./tan(betaAlongBladeW)).^2);
c=sqrt(cm.^2+(ufun(R)-cm./tan(betaAlongBladeW)).^2);
% c=sqrt(ufun(R).^2+w.^2-2.*ufun(R).*w.*cos(betaAlongBladeW))
cu=ufun(R)-cm./tan(betaAlongBladeW);
%plotting c, w, and UCu
% ucu
figure(figno)
figname.UcU=figno;
figno=figno+1;
plot(G./max(G),(ufun(R).*cu)./max(ufun(R).*cu))
xlabel('G/G_{max}')
ylabel('u \cdot c_u/max(u \cdot c_u)')
title('u \cdot c_u distribution')
%w
figure(figno)
figname.w=figno;
figno=figno+1;
plot(G./max(G),w./Values.DimensionlessVelocity)
xlabel('G/G_{max}')
ylabel('w/u_1')
title('w distribution')
%c
figure(figno)
figname.c=figno;
figno=figno+1;
plot(G./max(G),c./Values.DimensionlessVelocity)
xlabel('G/G_{max}')
ylabel('c/u_1')
title('c distribution')
%cm
figure(figno)
figname.cm=figno;
figno=figno+1;
plot(G./max(G),cm./Values.DimensionlessVelocity)
xlabel('G/G_{max}')
ylabel('c_m/u_1')
title('c_m distribution')
%u
figure(figno)
figname.u=figno;
figno=figno+1;
plot(G./max(G),ufun(R)./Values.DimensionlessVelocity)
xlabel('G/G_{max}')
ylabel('u/u_1')
title('u distribution')
%% H-streamline
H=G;
switch betaEval
    case 'Node'
DeltaH=DeltaG./tan(betaAlongBlade(1:end-1,:));%Not including the last beta, because I get trouble with the sizes of the matrices. maybe i should take the mean values of beta?
    case 'lineSegment'
        DeltaH=DeltaG./tan(betaAlongBlade);
end
for i=2:iEllipse
    H(i,:)=H(i-1,:)+DeltaH(i-1,:);
end
%% G-H plane
figure(figno)
figname.GHplane=figno;
figno=figno+1;
% for j=1:1:jEllipse %Debugging
%     plot(H(:,j),G(:,j),'-o')
% end
plot(H./max(max(H)),G./max(max(G)),'o-')
hold on
xlabel('H/H_{max}')
ylabel('G/G_{max}')
title('Plot of the G-H plane')
%% Polar coordinates
 %deltaTheta=(DeltaH./R(1:end-1,:));
deltaTheta=2.*asin((DeltaH./2./R(2:end,:)));%More accurate in my opinion
%deltaTheta=getDeltaTheta(R,DeltaH);
theta=zeros(size(R));%Initializing theta
for i=2:iEllipse
    theta(i,:)=theta(i-1,:)+deltaTheta(i-1,:);%finding theta for the domain by setting the initial value of theta equal to zero
end
figure(figno)
figname.RThetaPlane=figno;
figno=figno+1;
polarplot(theta,R./Values.DimensionlessRadius)%Plotting theta in polar coordinates
title('Theta-R/R_{max} in polar coordinates')
%% Back to cartesian coordinates
BladesPlotted='all';
figure(figno)
figname.streamPlane=figno;
figno=figno+1;
hold on
grid on
switch BladesPlotted
    case 'one'
        X=R.*cos(theta);
        Y=R.*sin(theta);
        mesh(X./Values.DimensionlessRadius,Y./Values.DimensionlessRadius,Z./Values.DimensionlessRadius)%(rad2deg(betaAlongBlade)>=90).*2)%or surf. Whatever you prefer
    case 'all'
        for i=0:1:(zla-1)
        X=R.*cos(theta+i*2*pi/zla);
        Y=R.*sin(theta+i*2*pi/zla);
        surf(X./Values.DimensionlessRadius,Y./Values.DimensionlessRadius,Z./Values.DimensionlessRadius)%(rad2deg(betaAlongBlade)>=90).*2)%or surf. Whatever you prefer
        end
        [XH,YH,ZH]=cylinder(RGeom(:,1));%Coords for hub
        [XS,YS,ZS]=cylinder(RGeom(:,jEllipse));%Coords for shroud
        surf(XH./Values.DimensionlessRadius,YH./Values.DimensionlessRadius,ZGeom(:,1).*ones(1,21)./Values.DimensionlessRadius)
        surf(XS./Values.DimensionlessRadius,YS./Values.DimensionlessRadius,ZGeom(:,jEllipse).*ones(1,21)./Values.DimensionlessRadius)
    otherwise
        error('ERROR: unrecognised value of BladesPlotted')
end
xlabel('X/R_1')
ylabel('Y/R_1')
zlabel('Z/R_1')
title('3D-plot of the streamlines')
%% Plot of the UcU dist. and the beta dist.

figure(figno)
figname.betaDist=figno;
figno=figno+1;
hold on
% yyaxis left
switch betaEval
    case 'Node'
for j=1:jEllipse
plot(G(:,j)./G(iEllipse,j),rad2deg(betaAlongBlade(:,j)));
end
    case 'lineSegment'
        for j=1:jEllipse
            plot((G(1:end-1,j)+DeltaG(:,j)./2)./G(iEllipse,j),rad2deg(betaAlongBlade(:,j)));
        end
    otherwise
        error('ERROR: no matching betaEval')
end
% yyaxis right
% plot(UcUAlongBlade)
%title('Plot of the UcU distribution and the beta distribution')
xlabel('G/G_{max}')
ylabel('\beta')
title('Beta distribution')

%Write to readable file for DM:
%txtForDesignModeler(X,Y,Z,R,'HubAndShroud')

% M=G(iEllipse,1)-G(:,1);
% writeAnglesToDm(M,betaAlongBlade(:,1),'BetaShroud')
% 
% M=G(iEllipse,jEllipse)-G(:,jEllipse);
% writeAnglesToDm(M,betaAlongBlade(:,jEllipse),'BetaHub')
% 
% M=G(iEllipse,ceil(jEllipse/2))-G(:,ceil(jEllipse/2));
% writeAnglesToDm(M,betaAlongBlade(:,ceil(jEllipse/2)),'BetaMiddle')


% M=G(iEllipse,1)-G(:,1);
%writeAnglesToDm(M,theta(:,1),'ThetaShroud',dir)
%writeThicknessToDM(M,e,'Shroud',dir)

% rM=G(iEllipse,jEllipse)-G(:,jEllipse);
%writeAnglesToDm(M,theta(:,jEllipse),'ThetaHub',dir)
%writeThicknessToDM(M,e,'Hub',dir)


% M=G(iEllipse,ceil(jEllipse/2))-G(:,ceil(jEllipse/2));
% writeAnglesToDm(M,theta(:,ceil(jEllipse/2)),'ThetaMiddle',dir)
% writeThicknessToDM(M,e,'Middle',dir)

fclose('all');
% The shapes of the trailing and leading edges.
LEshape
TEshape
figname.LETE=figno;
figno=figno+1;
%BladeCurveX=zeros(jEllipse,((iEllipse-4)+13)*2+1);
curveStart=3;
curveEnd=2;
BladeCurveX=Xup(curveStart:1:end-curveEnd,:);
for i=[1,29:-1:25,3,30:1:34,2]
    BladeCurveX=[BladeCurveX;TE(:,1,i)'];
end
BladeCurveX=[BladeCurveX;Xdown(end-curveEnd:-1:curveStart,:)];
for i=[2,34:-1:30,3,25:1:29,1]%[2,34:-1:30,3,25:1:29,1]
    BladeCurveX=[BladeCurveX;LE(:,1,i)'];
end
BladeCurveX=[BladeCurveX;Xup(curveStart,:)];

LeIndex=iEllipse-(curveStart-1)-curveEnd+1+5+1;%number of elements in a str.line-elements removed at start and end+elements up to, and including the LE
TeIndex=LeIndex+5+1+iEllipse-curveEnd-(curveStart-1)+1+5+1;%Elements to the LE, +elements until streamline starts again,-elements not included in streamline+elements at TE, including TE.

BladeCurveY=Yup(curveStart:1:end-curveEnd,:);
for i=[1,29:-1:25,3,30:1:34,2]
    BladeCurveY=[BladeCurveY;TE(:,2,i)'];
end
BladeCurveY=[BladeCurveY;Ydown(end-curveEnd:-1:curveStart,:)];
for i=[2,34:-1:30,3,25:1:29,1]%[2,34:-1:30,3,25:1:29,1]
    BladeCurveY=[BladeCurveY;LE(:,2,i)'];
end
BladeCurveY=[BladeCurveY;Yup(curveStart,:)];

BladeCurveZ=Zup(curveStart:1:end-curveEnd,:);
for i=[1,29:-1:25,3,30:1:34,2]
    BladeCurveZ=[BladeCurveZ;TE(:,3,i)'];
end
BladeCurveZ=[BladeCurveZ;Zdown(end-curveEnd:-1:curveStart,:)];
for i=[2,34:-1:30,3,25:1:29,1]%[2,34:-1:30,3,25:1:29,1]
    BladeCurveZ=[BladeCurveZ;LE(:,3,i)'];
end
BladeCurveZ=[BladeCurveZ;Zup(curveStart,:)];

figure(figno)
figname.curveSurf=figno;
figno=figno+1;
surf(BladeCurveX./Values.DimensionlessRadius,BladeCurveY./Values.DimensionlessRadius,BladeCurveZ./Values.DimensionlessRadius)%,'b')
hold on
plot3(BladeCurveX(LeIndex,:)./Values.DimensionlessRadius,BladeCurveY(LeIndex,:)./Values.DimensionlessRadius,BladeCurveZ(LeIndex,:)./Values.DimensionlessRadius,'ro')
plot3(BladeCurveX(TeIndex,:)./Values.DimensionlessRadius,BladeCurveY(TeIndex,:)./Values.DimensionlessRadius,BladeCurveZ(TeIndex,:)./Values.DimensionlessRadius,'rx')
cb=colorbar;
cb.Label.String='Z/D_1';
xlabel('X/R_1')
ylabel('Y/R_1')
title('The blade')
axis equal
view(2)
%% Should maybe be commented out, to avoid changing the geometry used in Ansys, unintentionally.
   makeCurveFile
   getVolume
%%
% figure(figno)
% figname.Edges=figno
% figno=figno+1;
% hold on
% plot3(LE3.X,LE3.Y,LE3.Z,'bo')
% plot3(LE.X,LE.Y,LE.Z,'ro')
% xlabel('x')
% ylabel('y')
% hold off
%% BC for cfx and other output:
clc
fprintf('Pressure at 1: \t \t\t\t\t\t\t\t%.2f \t\t[kPa]\n',rho*g*(Head+hf(Q)-hBooster(Q))/1000)
fprintf('Pressure at 2: \t\t\t\t\t\t\t\t%.2f \t\t[kPa] \n',(hBooster(Q)+HeadDiffSubmergence-hfDraftTube(Q)).*rho.*g./1000)
fprintf('Volume flow rate: \t \t\t\t\t\t\t%.2f \t\t[m^3/s]\n',Q)
fprintf('Meridional velocity at 1 boundary: \t\t\t%.2f \t\t[m/s]\n',Q/(pi*2*R1GuideVanes*b))
fprintf('Meridional velocity at 2 boundary: \t\t\t%.2f\t\t[m/s]\n',Q/(pi*(D2.^2)./4))
fprintf('Mass flow rate: \t\t\t\t\t\t\t%.2f \t[kg/s]\n',Q*rho)
fprintf('Greatest value for c: \t\t\t\t\t\t%.2f \t\t[m/s]\n',max(max(c)))
fprintf('Greatest value for w: \t\t\t\t\t\t%.2f \t\t[m/s]\n',max(max(w)))
fprintf('Required head:\t\t\t\t\t\t\t\t%.2f \t\t[m]\n',Hreq)
fprintf('Time for one passing is:\t\t\t\t\t%.4f \t\t[s]\n',60/n/zla)
% fprintf('Velocity at RGuidevanes=%.2f m:\t\t\t\t%.2f\t\t[m/s]\n',RGuidevanes,Q/(pi*(RGuidevanes(end)+RGeom(1,1))*2*(ZGeom(1,1)-ZGeom(1,jEllipse))))
if exist('V','var')
    fprintf('\nApproximate volume of the geometry: \t\t%.2f\t\t[m^3]\n',V)
    fprintf('Time for flow to go from inlet to outlet:\t%.2f\t\t[s]\n',V/Q)
end

deHaller=(Q/(pi*2*R1GuideVanes*b)*sin(deg2rad(betaDim)))/(Q/(pi*(D2.^2)./4));%(DIXON page 85&86)
if deHaller<0.72
    fprintf('\nWARNING: The DeHaller criterion is not fulfilled: (c1/c2) = %.2f <0.72\n',deHaller)
    beep()%Very important. makes a beeping sound.
else
    fprintf('\nThe DeHaller criterion is fulfilled: (c1/c2) = %.2f >=0.72\n',deHaller)
end
fprintf('Smallest near wall grid size dS wrt. y+:\t%.2e \t[m]\n',min(min(getGridSize(200,rho,nu,max(max(G)),w))))
switch addDraftTube
    case 'yes'
        AIn=pi.*(max(ZDraft).*tan(draftTubeAngle)+RGeom(iEllipse,jEllipse)).^2;
        VIn=Q./AIn;
        PtotIn=patm+rho*g*(HeadDiffSubmergence+hBooster(Q)-hfDraftTube(Q)+LDraftTube);
        REIn=rho*VIn*sqrt(AIn/pi)*2/mu;%Reynolds number at inlet to draft tube
        if REIn>2300
            PstatIn=PtotIn-Values.TurbulentEnergyCorrectionFactor.*rho.*VIn.^2./2;
        else
            PstatIn=PtotIn-Values.LaminarEnergyCorrectionFactor.*rho.*VIn.^2./2;
        end
        fprintf('Static pressure at inlet:\t\t\t\t\t%.2f\t\t[kPa]\n',PstatIn./1000)
end
save('P.mat','P')%saves the struct with the points for later reference. Should be moved to individual folder if you decide to run a simulation with this design.
NPSHA=getNPSHA(hBooster(Q),patm,pvapour,HeadDiffSubmergence,'pump');
NPSHR=getNPSHR(cm,ufun(R),'pump');
if NPSHA<=max(max(NPSHR))
    fprintf('Warning, cavitation will occur at some point\nNPSH_A= %.2f < %.2f =max(NPSH_R)\n',NPSHA,max(max(NPSHR)))
else
    fprintf('NPSHA= %.2f > %.2f =max(NPSHR)\n',NPSHA,max(max(NPSHR)))
end
fprintf('Smallest boundary layer thickness:\t\t\t%f\t[m]',min(min(getBLThickness(w,(max(G)-G)))))