function txtForDesignModeler(X,Y,Z,R,type)
%type='HubAndShroud';
switch type
    case 'Streamlines'
    fid=fopen('writeToDM.txt','w');
   % try
        fprintf(fid,'#------------------------------------------------------------\n');
        fprintf(fid,'# List of Point Coordinates\n');
        fprintf(fid,'# Format is: integer Group, integer ID, X, Y, Z\n');
        fprintf(fid,'# all delimited by spaces, with nothing after the Z value.\n');
        for j=1:size(X,2)
            fprintf(fid,'# Group %i\n',j);
            for i=1:size(X,1)
            fprintf(fid,'%i %i %2.3f %2.3f %2.3f\n',j,i,X(i,j),Y(i,j),Z(i,j));
            end
            fprintf(fid,'\n');
        end
        fprintf(fid,'#------------------------------------------------------------');
        fclose(fid);
   % catch
    fclose(fopen('all'));
        error('ERROR')
    %end
    case 'HubAndShroud'
         fid=fopen('writeShroudAndHubToDM.txt','w');
   % try
        fprintf(fid,'#------------------------------------------------------------\n');
        fprintf(fid,'# List of Point Coordinates\n');
        fprintf(fid,'# Format is: integer Group, integer ID, X, Y, Z\n');
        fprintf(fid,'# all delimited by spaces, with nothing after the Z value.\n');
        for j=[1,size(R,2)]
            fprintf(fid,'%i %i %2.3f 0 %2.3f\n',j,1,R(1,j)+0.5,Z(1,j));
            for i=1:size(R,1)
            fprintf(fid,'%i %i %2.3f 0 %2.3f\n',j,i+1,R(i,j),Z(i,j));
            end
            fprintf(fid,'%i %i %2.3f 0 %2.3f\n',j,(numel(R)./length(R))+2,R(end,j),Z(end,j)-0.5);
            fprintf(fid,'\n');
        end
        fprintf(fid,'#trailing edge and leading edge:\n');
        for i=[1,size(R,1)]
            for j=1:1:size(R,2)
                fprintf(fid,'%i %i %2.3f 0 %2.3f\n',length(R)+i,j,R(i,j),Z(i,j));
            end
                fprintf(fid,'\n');
        end
        fprintf(fid,'#------------------------------------------------------------');
        fclose(fid);
    %catch
   % fclose(fopen('all'));
%    error('ERROR')
    %end
end
end