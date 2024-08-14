%makeGuideVaneCurves
%% AngularTransformFunction
%first we make the shape. zero at sharp edge
P.GuideVaneSs=zeros(19,2);
P.GuideVaneSs(:,1)=linspace(0,540,19)./1000;
P.GuideVaneSs(:,2)=[];
P.GuideVanePs=P.GuideVaneSs(end-1:-1:1,:);
P.GuideVanePs(:,2)=[];
%finding the geometry of th LE
Circle.radius=;
Circle.Centre=-Circle.radius;
Circle.SSx=-;
Circle.SSy=;
Circle.PSS=sqrt(Circle.SSx^2+Circle.SSy^2);
Circle.tangentSS=sqrt(Circle.PSS^2-Circle.radius^2);
Circle.theta1SS=atan(Circle.SSy/Circle.SSx)+pi;
Circle.theta2SS=acos(Circle.radius/Circle.PSS);
Circle.SSStart=Circle.radius.*cos(Circle.theta1SS-Circle.theta2SS)+Circle.Centre;
P.TipSs(:,1)=linspace(Circle.SSStart,,10);
P.TipSs(:,2)=sqrt(Circle.radius.^2-(P.TipSs(:,1)-Circle.Centre).^2);

Circle.PSx=;
Circle.PSy=;
Circle.PPS=sqrt(Circle.PSx^2+Circle.PSy^2);
Circle.tangentPS=sqrt(Circle.PPS^2-Circle.radius^2);
Circle.theta1PS=atan(Circle.PSy/Circle.PSx)+pi;
Circle.theta2PS=acos(Circle.radius/Circle.PPS);
Circle.PSStart=Circle.radius.*cos(Circle.theta1PS-Circle.theta2PS)+Circle.Centre;
P.TipPs(:,1)=linspace(Circle.PSStart,,10);
P.TipPs(:,2)=-sqrt(Circle.radius.^2-(P.TipPs(:,1)-Circle.Centre).^2);

%P.GuideVane=[P.GuideVanePs(end:-1:1,:);P.TipPs;P.TipSs(end:-1:1,:);P.GuideVaneSs(end-1:-1:1,:)];
P.GuideVane=real([P.GuideVaneSs(1:end-1,:);P.TipSs;P.TipPs(end:-1:1,:);P.GuideVanePs]);
plot(P.GuideVane(:,1),P.GuideVane(:,2))
axis equal
C(1)=;
C(2)=C(1)*;
%Transforming to cylinder coordinates
[P.GuideVane(:,1),P.GuideVane(:,2)]=angularTransform(P.GuideVane(:,1),real(P.GuideVane(:,2)),C(1),C(2));
C(1)=;%placing the tap where it should be
alpha=pi/3;%atan((Q/A)./(ufun(D1/2)-(Q/A)/tan(beta1Slip)));%rotation around the tap in radians
Zvalue=[];%at the hub, in the middle and at the shroud
%Transforming to cartesianCoordinates
[P.GuideVane(:,1),P.GuideVane(:,2)]=cartesianTransform(P.GuideVane(:,1),P.GuideVane(:,2)+alpha,C(1),C(2));
plot(P.GuideVane(:,1),P.GuideVane(:,2))
axis equal
GuideVanePoints=real(P.GuideVane);
P.TipMin=C(1)-;
%% blade.crv
% Format: X, Y and Z. leading edge and trailing edge have a fourth column
% indexed by "le1" and "te1",respectively. Blade curves are separated by a
% single line with a #.
%The curves have to start at the hub and go towards the shroud
fid=fopen([dir,'/guideVane_blade.crv'],'w');
for j=1:length(Zvalue)
    fprintf(fid,'#\n');
    for i=1:1:length(GuideVanePoints)
       fprintf(fid,'%f\t%f\t%f',GuideVanePoints(i,1),GuideVanePoints(i,2),Zvalue(j));
       if i==(length(P.GuideVaneSs)+length(P.TipSs)-1)
         fprintf(fid,'\tte1');
       elseif i==length(GuideVanePoints)
           fprintf(fid,'\\tle1');
       end
       fprintf(fid,'\n');
    end
end
fclose(fid);
% %% hub.crv
% % format: The axial projection of the shroud and hub, respectively(in
% % R,theta,Z coordinates, but theta may remain 0)
% fid=fopen([dir,'/guideVane_hub.crv'],'w');
% fprintf(fid,'%f\t%f\t%f\n',C(1)-0.5,0,Zvalue(1));
% fprintf(fid,'%f\t%f\t%f\n',C(1)+0.5,0,Zvalue(1));
% fclose(fid);
% %% shroud.crv %Mostly the same as for hub.
% fid=fopen([dir,'/guideVane_shroud.crv'],'w');
% fprintf(fid,'%f\t%f\t%f\n',C(1)-0.5,0,Zvalue(3));
% fprintf(fid,'%f\t%f\t%f\n',C(1)+0.5,0,Zvalue(3));
% fclose(fid);