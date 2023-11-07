
function [Rout,Zout]=cutOffBlade(R,Z,TE,LE,iEllipse,jEllipse)
TEpchip=extractBezierByPoints(0,TE);
Rout=zeros(size(R));
Zout=zeros(size(Z));
    for j=1:1:jEllipse
        %Cut off at TE, first
        Zpchip=pchip(R(:,j),Z(:,j));                                 %Piecewise polynomial of streamline wrt r
        rr=linspace(max(R(end,j),min(TE(:,1))),min(R(1,j),max(TE(:,1))));   %Interval where both streamline and TE are defined
        zDiff=ppval(Zpchip,rr)-ppval(TEpchip,rr);
        zDiffpchip=pchip(rr,zDiff);
        Rcut=unique(fnzeros(zDiffpchip));
        rz=interparc(iEllipse,linspace(R(end,j),Rcut),ppval(Zpchip,linspace(R(end,j),Rcut)),'pchip');
        Rout(:,j)=rz(:,1);
        Zout(:,j)=rz(:,2);
    end
    R=Rout(end:-1:1,:);
    Z=Zout(end:-1:1,:);
    %Switch axis and cut off LE. Had problems with this one due to vertical
    %lines at the outlet.
    LEpchip=extractBezierByPoints(0,LE(:,end:-1:1));
    hold on
        for j=1:1:jEllipse
            Rpchip=pchip(Z(:,j),R(:,j));                                 %Piecewise polynomial of streamline wrt r
            zz=linspace(max(Z(end,j),min(LE(end,2),LE(1,2))),min(Z(1,j),max(LE(1,2),LE(end,2))));   %Interval where both streamline and TE are defined
            rDiff=ppval(Rpchip,zz)-ppval(LEpchip,zz);
            rDiffpchip=pchip(zz,rDiff);
            Zcut=unique(fnzeros(rDiffpchip));
            zr=interparc(iEllipse,linspace(Zcut,Z(1,j)),ppval(Rpchip,linspace(Zcut,Z(1,j))),'pchip');
            Rout(:,j)=zr(:,2);
            Zout(:,j)=zr(:,1);
        end
    Rout=Rout(end:-1:1,:);
    Zout=Zout(end:-1:1,:);
end