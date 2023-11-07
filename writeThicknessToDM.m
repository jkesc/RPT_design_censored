function writeThicknessToDM(M,e,curve,dir)
fid=fopen([dir,'/Thickness',curve,'DM.ht'],'w');
fprintf(fid,'2\n');
fprintf(fid,'%f\t%f\n',M(1),e);
fprintf(fid,'%f\t%f\n',M(end),e);
fclose(fid);
end