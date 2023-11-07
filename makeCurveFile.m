%makeCurveFile
%(this is the one used in BladeDesign2 written by Karl E. Not to be confused
%with makeCurveFiles written by Helene N.D.
%% blade.crv
% Format: X, Y and Z. leading edge and trailing edge have a fourth column
% indexed by "le1" and "te1",respectively. Blade curves are separated by a
% single line with a #.
%The curves have to start at the hub and go towards the shroud
switch addGuideVaneGeometry
            case 'yes'
                makeGuideVaneCurves
    otherwise
end
fidcsv=fopen([dir,'/blade.csv'],'w');
fidBlade=fopen([dir,'/blade.crv'],'w');                                     %Opening a txt file for the blade
for j=(1:1:size(BladeCurveX,2))                                             %For number of streamlinse
    fprintf(fidBlade,'#\n');                                                %# to seperate blades
    for i=1:1:length(BladeCurveX)                                           %For number of elements in each curve
        if i==LeIndex                                                       %If this is the curve element corresponding to le
            string="le1";                                                   %write this in the file
        elseif i==TeIndex                                                   %same for te
                string="te1";
        else
            string="";                                                      %otherwise, just don't write anything
        end
        fprintf(fidBlade,'%f\t%f\t%f\t%s\n',BladeCurveX(i,j),BladeCurveY(i,j),BladeCurveZ(i,j),string);%print out the corresponding X Y and Z values
        fprintf(fidcsv,'%f,%f,%f\n',BladeCurveX(i,j),BladeCurveY(i,j),BladeCurveZ(i,j));%Trying to create a polysurface for use in cfd Post
    end
end
fclose(fidBlade);                                                           %Close the file
fclose(fidcsv);
%% hub.crv
% format: The axial projection of the shroud and hub, respectively(in
% R,theta,Z coordinates, but theta may remain 0)
figure(figno)                                                               %just to see if it does what it should
figname.HubShroud=figno;
figno=figno+1;
hold on
fid=fopen([dir,'/hub.crv'],'w');                                            %open a hub file
j=1;                                                                        %initialize j
i=iEllipse;                                                                 %i starts at the end, so at the circular outlet
switch addDraftTube
    case 'yes'
        LDraftTube=RGeom(iEllipse,jEllipse)*10;
        ZDraft=linspace(LDraftTube,LDraftTube/20,20);
        draftTubeAngle=deg2rad(4);
        for iDraft=1:length(ZDraft)
            fprintf(fid,'%f\t%f\t%f\n',R(iEllipse,j),0,ZGeom(i,j)-ZDraft(iDraft));
            plot(R(iEllipse,j)./Values.DimensionlessRadius,(ZGeom(i,j)-ZDraft(iDraft))./Values.DimensionlessRadius,'o')
        end
    otherwise
end
while(ZGeom(i,j)<Z(iEllipse,j))                                             %while the hub has not reached the blade yet
    fprintf(fid,'%f\t%f\t%f\n',R(iEllipse,j),0,ZGeom(i,j));                    %write the coordinates of the hub
    plot(R(iEllipse,j)./Values.DimensionlessRadius,ZGeom(i,j)./Values.DimensionlessRadius,'o')
    i=i-1;
end
for i=iEllipse:-1:1                                                         %Then write all the elements along the blade(from inlet to outlet)
       fprintf(fid,'%f\t%f\t%f\n',R(i,j),0,Z(i,j));
       plot(R(i,j)./Values.DimensionlessRadius,Z(i,j)./Values.DimensionlessRadius,'o')
end
iloop=iEllipse;                                                            %index to find when the hub is not on the blade anymore
while(RGeom(iloop,j)<R(1,j))                                               %and corresponding loop
    iloop=iloop-1;
end
for i=iloop:-1:1                                                            %Print the remainder of the hub
    switch addGuideVanes
        case 'yes'
        switch addGuideVaneGeometry
            case 'yes'
                if RGeom(i,j)<P.TipMin
               fprintf(fid,'%f\t%f\t%f\n',RGeom(i,j),0,ZGeom(i,j));
               plot(RGeom(i,j)./Values.DimensionlessRadius,ZGeom(i,j)./Values.DimensionlessRadius,'o')
                end
            otherwise
               fprintf(fid,'%f\t%f\t%f\n',RGeom(i,j),0,ZGeom(i,j));
               plot(RGeom(i,j)./Values.DimensionlessRadius,ZGeom(i,j)./Values.DimensionlessRadius,'o')
        end
        otherwise
               fprintf(fid,'%f\t%f\t%f\n',RGeom(i,j),0,ZGeom(i,j));
               plot(RGeom(i,j)./Values.DimensionlessRadius,ZGeom(i,j)./Values.DimensionlessRadius,'o')
     end
end
switch addGuideVanes
    case 'yes'
        switch addGuideVaneGeometry
            case 'yes'
                fidGV=fopen([dir,'/hubGV.crv'],'w');
                RGuidevanes=[P.TipMin,P.TipMin+0.54,P.TipMin+1];
                fprintf(fid,'%f\t%f\t%f\n',RGuidevanes(1),0,ZGeom(i,j));
                plot(RGuidevanes(iGuidevanes)./Values.DimensionlessRadius,ZGeom(i,j)./Values.DimensionlessRadius,'o');
                for iGuidevanes=1:length(RGuidevanes)
                    fprintf(fidGV,'%f\t%f\t%f\n',RGuidevanes(iGuidevanes),0,ZGeom(i,j));
                    plot(RGuidevanes(iGuidevanes)./Values.DimensionlessRadius,ZGeom(i,j)./Values.DimensionlessRadius,'o');
                end
                fclose(fidGV);
            otherwise
        RGuidevanes=Values.RGuidevanes;
        AGuidevanes=2*pi*RGeom(1,j)*(Z(1,1)-Z(1,jEllipse));%want to keep the area after guide vanes constant to avoid separation
        DeltaZGuidevanes=AGuidevanes./(2.*pi.*(RGeom(1,j)))-AGuidevanes./(2.*pi.*(RGeom(1,j)+RGuidevanes)).*(1-RGuidevanes./max(RGuidevanes).*0.1);%Making the area smaller with increasing radius to avoid recirculation
        for iGuidevanes=1:length(RGuidevanes)
            fprintf(fid,'%f\t%f\t%f\n',RGeom(i,j)+RGuidevanes(iGuidevanes),0,ZGeom(i,j)-DeltaZGuidevanes(iGuidevanes)./2);
            plot((RGeom(i,j)+RGuidevanes(iGuidevanes))./Values.DimensionlessRadius,(ZGeom(i,j)-DeltaZGuidevanes(iGuidevanes)./2)./Values.DimensionlessRadius,'o');
        end
        end
    otherwise
end
% fprintf(fid,'%f\t%f\t%f\n',RGeom(i,j)+1,0,ZGeom(i,j));
% plot(RGeom(i,j)+1,ZGeom(i,j),'o');
fclose(fid);
%% shroud.crv %Mostly the same as for hub.
fid=fopen([dir,'/shroud.crv'],'w');
j=jEllipse;
% fprintf(fidShroud,'%f\t%f\t%f\n',R(iEllipse,j),0,Z(iEllipse,j)-0.5);
% for i=iEllipse:-1:1
%         fprintf(fidHub,'%f\t%f\t%f\n',R(i,j),0,Z(i,j));
% end
% fprintf(fidHub,'%f\t%f\t%f\n',R(1,j)+0.5,0,Z(1,j));
i=iEllipse;
switch addDraftTube
    case 'yes'
        for iDraft=1:length(ZDraft)
            RDraft=ZDraft(iDraft).*tan(draftTubeAngle);
            fprintf(fid,'%f\t%f\t%f\n',RGeom(i,j)+RDraft,0,ZGeom(i,j)-ZDraft(iDraft));
            plot((RGeom(i,j)+RDraft)./Values.DimensionlessRadius,(ZGeom(i,j)-ZDraft(iDraft))./Values.DimensionlessRadius,'o')
        end
    otherwise
end
while(ZGeom(i,j)<Z(iEllipse,j))
    fprintf(fid,'%f\t%f\t%f\n',RGeom(i,j),0,ZGeom(i,j));
    plot(RGeom(i,j)./Values.DimensionlessRadius,ZGeom(i,j)./Values.DimensionlessRadius,'o')
    i=i-1;
end
for i=iEllipse:-1:1
       fprintf(fid,'%f\t%f\t%f\n',R(i,j),0,Z(i,j));
       plot(R(i,j)./Values.DimensionlessRadius,Z(i,j)./Values.DimensionlessRadius,'o')
end
iloop=iEllipse;
while(RGeom(iloop,j)<R(1,j))
    iloop=iloop-1;
end
switch addGuideVaneGeometry
    case 'yes'
    fidGV=fopen([dir,'/shroudGV.crv'],'w');
    otherwise
end
for i=iloop:-1:1
    switch addGuideVanes
        case 'yes'
         switch addGuideVaneGeometry
            case 'yes'
                if i==iloop
                    fprintf(fid,'%f\t%f\t%f\n',RGeom(i,j),0,ZGeom(i,j));
                    plot(RGeom(i,j)./Values.DimensionlessRadius,ZGeom(i,j)./Values.DimensionlessRadius,'o')
                end
                fprintf(fidGV,'%f\t%f\t%f\n',RGeom(i,j),0,ZGeom(i,j));
                plot(RGeom(i,j)./Values.DimensionlessRadius,ZGeom(i,j)./Values.DimensionlessRadius,'o')
            otherwise
            fprintf(fid,'%f\t%f\t%f\n',RGeom(i,j),0,ZGeom(i,j));
            plot(RGeom(i,j)./Values.DimensionlessRadius,ZGeom(i,j)./Values.DimensionlessRadius,'o')
         end
        otherwise
        fprintf(fid,'%f\t%f\t%f\n',RGeom(i,j),0,ZGeom(i,j));
        plot(RGeom(i,j)./Values.DimensionlessRadius,ZGeom(i,j)./Values.DimensionlessRadius,'o')
    end
end
switch addGuideVanes
    case 'yes'
         switch addGuideVaneGeometry
            case 'yes'
%                 fprintf(fid,'%f\t%f\t%f\n',RGuidevanes(2),0,ZGeom(i,j));
%                 plot(RGuidevanes(iGuidevanes)./Values.DimensionlessRadius,ZGeom(i,j)./Values.DimensionlessRadius,'o');
                for iGuidevanes=2:length(RGuidevanes)
                    fprintf(fidGV,'%f\t%f\t%f\n',RGuidevanes(iGuidevanes),0,ZGeom(i,j));
                    plot(RGuidevanes(iGuidevanes)./Values.DimensionlessRadius,ZGeom(i,j)./Values.DimensionlessRadius,'o');
                end
                fclose(fidGV);
            otherwise
                for iGuidevanes=1:length(RGuidevanes)
                    fprintf(fid,'%f\t%f\t%f\n',RGeom(i,j)+RGuidevanes(iGuidevanes),0,ZGeom(i,j)+DeltaZGuidevanes(iGuidevanes)./2);
                    plot((RGeom(i,j)+RGuidevanes(iGuidevanes))./Values.DimensionlessRadius,(ZGeom(i,j)+DeltaZGuidevanes(iGuidevanes)./2)./Values.DimensionlessRadius,'o');
                end
         end
    otherwise
end
axis equal
% fprintf(fid,'%f\t%f\t%f\n',RGeom(i,j)+1,0,ZGeom(i,j));
% plot(RGeom(i,j)+1,ZGeom(i,j),'o');
fclose(fid);
