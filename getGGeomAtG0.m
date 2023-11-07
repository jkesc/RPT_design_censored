function GToAdd = getGGeomAtG0(Rintersection,Zintersection,RGeom,ZGeom)
%gets the value of GGeom where G starts.
GToAdd=zeros(1,size(RGeom,2));                                                                  %initializing GToAdd
for j=1:size(GToAdd,2)                                                                           %Iterating through the lines
    i=2;                                                                                    %Starting at i=2 to include i=1 at i-1
    while RGeom(i,j)>Rintersection                                                              %While we are away from the starting point
        GToAdd(1,j)=GToAdd(1,j)+sqrt((RGeom(i-1,j)-RGeom(i,j)).^2+(ZGeom(i-1,j)-ZGeom(i,j)).^2);                      %adding the length to the value of G
        i=i+1;                                                                              %incrementing index
    end
    GToAdd(1,j)=GToAdd(1,j)+sqrt((RGeom(i-1,j)-Rintersection(1,j)).^2+(ZGeom(i-1,j)-Zintersection(1,j)).^2);  %Adding the rest of the line
end
end
