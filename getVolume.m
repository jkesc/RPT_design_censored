%getVolume.m
%Finds the volume of the geometry, to use in transient simulation.
switch addDraftTube     %Approximation of volume from leading edge of the blade to the point added at the "inlet" in CFX
    case 'yes'
        VDraft=sum((Z(end,1:end-1)-(ZGeom(end,1:end-1)-max(ZDraft))).*pi.*(R(end,2:end).^2-R(end,1:end-1).^2));
    case 'no'
        VDraft=0;
    otherwise
        error('ERROR: No valid value for addDraftTube')
end
switch addGuideVanes    %Approximation of volume of the part after the blade and up to the "outlet" in CFX
    case 'yes'
        VGuide=sum(((RGeom(1,1:end-1)+RGuidevanes(end)).^2-R(1,1:end-1).^2).*(Z(1,1:end-1)-Z(1,2:end)+ZGeom(1,1:end-1)-ZGeom(1,2:end))./2.*pi);
    case 'no'
        VGuide=0;
    otherwise
        error('ERROR: No valid value for addGuideVanes')
end
dAijV=zeros(iEllipse,jEllipse-1);   %The cross sectional area at almost all streamlines
for j=1:1:jEllipse-1                
dAijV(:,j)=getdAij((G(:,j)+GToAdd(:,j)),dA1,dA2,P.Aij,GGeom(end,j));
end
Vmain=sum(sum((G(2:end,1:end-1)-G(1:end-1,1:end-1)).*getdAij((G(2:end,1:end-1)+GToAdd(:,1:end-1)),dA1,dA2,P.Aij,GGeom(end,1:end-1))));%multiplying the length of an element with its area to get the volume
V=Vmain+VDraft+VGuide;  %Adding all approximate volumes.
