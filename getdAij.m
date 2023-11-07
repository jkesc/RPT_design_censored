function dAij=getdAij(G,dA1,dA2,Pmid,Gmax)
if nargin==4%If G is a matrix/vector, we can use that to find the maximal value of G
    Gfrac=G./max(G);
else 
    Gfrac=G/Gmax;%Else we need a reference input to find G
end
P=[
    0         0
    Pmid
    1         1
    ];%Some points by which the acceleration is defined in this script. May be changed if necessary.
dAijChip=extractBezierByPoints(0,P);%Input these points to get a bezier curve as pchip.
dAij=dA1+(ppval(dAijChip,Gfrac)).*(dA2-dA1);%evaluate the curve at the points along the Length of the streamline.
end