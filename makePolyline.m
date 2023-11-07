%makePolyLine
%(this is the one used in BladeDesign2 written by Karl E. Not to be confused
%with makeCurveFiles written by Helene N.D.
%% blade.txt
% Format: X, Y and Z. leading edge and trailing edge have a fourth column
% indexed by "le1" and "te1",respectively. Blade curves are separated by a
% single line with a #.
%The curves have to start at the hub and go towards the shroud
fidtxt=fopen([dir,'/bladePolyline.txt'],'w');
for j=(1:1:size(BladeCurveX,2))                                             %For number of streamlinse
    fprintf(fidtxt,'[Name]\n');
    fprintf(fidtxt,'Polyline %i\n',j);
    fprintf(fidtxt,'[Data]\n');
    fprintf(fidtxt,'X [ m ], Y [ m ], Z [ m ], Area [ m^2 ], Density [ kg m^-3 ]\n');
    for i=1:1:length(BladeCurveX)                                           %For number of elements in each curve
        fprintf(fidtxt,'%f,%f,%f,0,0\n',BladeCurveX(i,j),BladeCurveY(i,j),BladeCurveZ(i,j));%Trying to create a polysurface for use in cfd Post
    end
    fprintf(fidtxt,'[Lines]\n');
    for i=1:1:length(BladeCurveX)-1
        fprintf(fidtxt,'%i,%i\n',i-1,i);
    end
    fprintf(fidtxt,'%i,%i\n',i,0);
end
fclose(fidtxt);