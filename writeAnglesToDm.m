function writeAnglesToDm(M,beta,curve,dir)
fid=fopen([dir,'/AngleCurve',curve,'DM.ha'],'w');
%try
fprintf(fid,'%i\n',length(beta));
for i=length(M):-1:1
    fprintf(fid,'%1.6f\t%1.6f\n',M(i),beta(i));
end
%DM requires a .ha-file with number of inputs as first line, and then tab
%separated values of M' and beta
%M' is the length along the blade nondimensionalized to the radius of the
%blade.
%M is the length of the blade (integrate(sqrt(dR^2+dZ^2)).
%M begins on Le and goes to Te. In Ansys-terms, this means from 2 to 1 in
%this script...
%Remember to change the input stuff to M in blade modeler.
%catch
    fclose('all');
 %   error('ERROR')
%end
end