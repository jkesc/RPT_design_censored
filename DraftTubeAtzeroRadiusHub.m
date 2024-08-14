%HubToZeroRadius
fid=fopen([dir,'/hubZeroRadius.crv'],'w');
j=1;                                                                        %initialize j
i=iEllipse;  
switch addDraftTube
    case 'yes'
        LDraftTube=RGeom(iEllipse,jEllipse)*10;
        ZDraft=linspace(LDraftTube,LDraftTube/20,20);
        draftTubeAngle=deg2rad(4);
        for iDraft=1:length(ZDraft)
            fprintf(fid,'%f\t%f\t%f\n',0,0,ZGeom(i,j)-ZDraft(iDraft));
            plot(0,ZGeom(i,j)-ZDraft(iDraft),'-x')
        end
        fprintf(fid,'%f\t%f\t%f\n',R(iEllipse,j),0,ZGeom(i,j)-ZDraft(iDraft));
            plot(R(iEllipse,j),ZGeom(i,j)-ZDraft(iDraft),'-x')
    otherwise
end
while(ZGeom(i,j)<Z(iEllipse,j))                                             %while the hub has not reached the blade yet
    fprintf(fid,'%f\t%f\t%f\n',R(iEllipse,j),0,ZGeom(i,j));                    %write the coordinates of the hub
    plot(R(iEllipse,j),ZGeom(i,j),'-x')
    i=i-1;
end
for i=iEllipse:-1:1                                                         %Then write all the elements along the blade(from inlet to outlet)
       fprintf(fid,'%f\t%f\t%f\n',R(i,j),0,Z(i,j));
       plot(R(i,j),Z(i,j),'-x')
end
iloop=iEllipse;                                                            %index to find when the hub is not on the blade anymore
while(RGeom(iloop,j)<R(1,j))                                               %and corresponding loop
    iloop=iloop-1;
end
for i=iloop:-1:1                                                            %Print the remainder of the hub
       fprintf(fid,'%f\t%f\t%f\n',RGeom(i,j),0,ZGeom(i,j));
       plot(RGeom(i,j),ZGeom(i,j),'-x')
end
switch addGuideVanes
    case 'yes'
        RGuidevanes=[0.5,1,1.5,2.0,2.5];
        AGuidevanes=2*pi*RGeom(1,j)*(Z(1,1)-Z(1,jEllipse));%want to keep the area after guide vanes constant to avoid separation
        DeltaZGuidevanes=AGuidevanes./(2.*pi.*(RGeom(1,j)))-AGuidevanes./(2.*pi.*(RGeom(1,j)+RGuidevanes));%The height of the shroud when the area is constant.
        for iGuidevanes=1:length(RGuidevanes)
            fprintf(fid,'%f\t%f\t%f\n',RGeom(i,j)+RGuidevanes(iGuidevanes),0,ZGeom(i,j)-DeltaZGuidevanes(iGuidevanes)./2);
            plot(RGeom(i,j)+RGuidevanes(iGuidevanes),ZGeom(i,j)-DeltaZGuidevanes(iGuidevanes)./2,'-x');
        end
    otherwise
end
% fprintf(fid,'%f\t%f\t%f\n',RGeom(i,j)+1,0,ZGeom(i,j));
% plot(RGeom(i,j)+1,ZGeom(i,j),'o');
fclose(fid);