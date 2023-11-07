%LEshape.m 
% The creation of the thickness will be as in the powerpoint supplied by
% Rainpower.
% Finding the normal to the surface created by the blades:

[Nx,Ny,Nz]=surfnorm(X,Y,Z);
Xup=X-Nx.*e./2;
Xdown=X+Nx.*e./2;
Yup=Y-Ny.*e./2;
Ydown=Y+Ny.*e./2;
Zup=Z-Nz.*e./2;
Zdown=Z+Nz.*e./2;
% Defining the ratios of the ellipses to be constructed later
bEllipse=1*e;
aEllipse=e/2;


LE=zeros(jEllipse,3,34);%will obtain coordinates for jEllipse streamlines, 3 dimensions and 37 points pr streamline.
aAbs=sqrt((Xup(1,:)-Xup(2,:)).^2+(Yup(1,:)-Yup(2,:)).^2+(Zup(1,:)-Zup(2,:)).^2);
bAbs=bEllipse;
LE(:,1,1)=-bAbs./aAbs.*(Xup(1,:)-Xup(2,:))+Xup(1,:);
LE(:,2,1)=-bAbs./aAbs.*(Yup(1,:)-Yup(2,:))+Yup(1,:);
LE(:,3,1)=-bAbs./aAbs.*(Zup(1,:)-Zup(2,:))+Zup(1,:);

aAbs=sqrt((Xdown(1,:)-Xdown(2,:)).^2+(Ydown(1,:)-Ydown(2,:)).^2+(Zdown(1,:)-Zdown(2,:)).^2);
bAbs=bEllipse;
LE(:,1,2)=-bAbs./aAbs.*(Xdown(1,:)-Xdown(2,:))+Xdown(1,:);
LE(:,2,2)=-bAbs./aAbs.*(Ydown(1,:)-Ydown(2,:))+Ydown(1,:);
LE(:,3,2)=-bAbs./aAbs.*(Zdown(1,:)-Zdown(2,:))+Zdown(1,:);

aAbs=sqrt((Xup(1,:)-Xdown(1,:)).^2+(Yup(1,:)-Ydown(1,:)).^2+(Zup(1,:)-Zdown(1,:)).^2);
bAbs=e/2;
LE(:,1,3)=-bAbs./aAbs.*(Xup(1,:)-Xdown(1,:))+Xup(1,:);
LE(:,2,3)=-bAbs./aAbs.*(Yup(1,:)-Ydown(1,:))+Yup(1,:);
LE(:,3,3)=-bAbs./aAbs.*(Zup(1,:)-Zdown(1,:))+Zup(1,:);

aAbs=sqrt((LE(:,1,1)-LE(:,1,2)).^2+(LE(:,2,1)-LE(:,2,2)).^2+(LE(:,3,1)-LE(:,3,2)).^2);
bAbs=e/2;
LE(:,:,4)=-bAbs./aAbs.*(LE(:,:,1)-LE(:,:,2))+LE(:,:,1);

% Should be able to construct both sides from knowledge of LE 1-4.
%upside:
bAbs=e/12;
LE(:,:,9)=LE(:,:,4);

%Points away from the tip
for i=5:1:9
    aAbs=sqrt((LE(:,1,1)-LE(:,1,i-1)).^2+(LE(:,2,1)-LE(:,2,i-1)).^2+(LE(:,3,1)-LE(:,3,i-1)).^2);
    LE(:,:,i)=bAbs./aAbs.*(LE(:,:,1)-LE(:,:,i-1))+LE(:,:,i-1);
    LE(:,:,i+5)=bAbs./aAbs.*(LE(:,:,2)-LE(:,:,(i+5)-1))+LE(:,:,(i+5)-1);
end

%Points on the tip, but flat.
%upside
bAbs=e/12;
aAbs=e/2;
LE(:,1,15)=bAbs./aAbs.*(Xup(1,:)'-LE(:,1,3))+LE(:,1,3);
LE(:,2,15)=bAbs./aAbs.*(Yup(1,:)'-LE(:,2,3))+LE(:,2,3);
LE(:,3,15)=bAbs./aAbs.*(Zup(1,:)'-LE(:,3,3))+LE(:,3,3);

for i=16:1:19
    aAbs=sqrt((Xup(1,:)'-LE(:,1,i-1)).^2+(Yup(1,:)'-LE(:,2,i-1)).^2+(Zup(1,:)'-LE(:,3,i-1)).^2);
    LE(:,1,i)=bAbs./aAbs.*(Xup(1,:)'-LE(:,1,i-1))+LE(:,1,i-1);
    LE(:,2,i)=bAbs./aAbs.*(Yup(1,:)'-LE(:,2,i-1))+LE(:,2,i-1);
    LE(:,3,i)=bAbs./aAbs.*(Zup(1,:)'-LE(:,3,i-1))+LE(:,3,i-1);
end
%downside
bAbs=e/12;
aAbs=e/2;
LE(:,1,20)=bAbs./aAbs.*(Xdown(1,:)'-LE(:,1,3))+LE(:,1,3);
LE(:,2,20)=bAbs./aAbs.*(Ydown(1,:)'-LE(:,2,3))+LE(:,2,3);
LE(:,3,20)=bAbs./aAbs.*(Zdown(1,:)'-LE(:,3,3))+LE(:,3,3);

for i=21:1:24
    aAbs=sqrt((Xdown(1,:)'-LE(:,1,i-1)).^2+(Ydown(1,:)'-LE(:,2,i-1)).^2+(Zdown(1,:)'-LE(:,3,i-1)).^2);
    LE(:,1,i)=bAbs./aAbs.*(Xdown(1,:)'-LE(:,1,i-1))+LE(:,1,i-1);
    LE(:,2,i)=bAbs./aAbs.*(Ydown(1,:)'-LE(:,2,i-1))+LE(:,2,i-1);
    LE(:,3,i)=bAbs./aAbs.*(Zdown(1,:)'-LE(:,3,i-1))+LE(:,3,i-1);
end

%Continue from here with the ellipsis-form
j=0;
aAbs=bEllipse;
for i=25:1:29
    j=j+1;
    bAbs=sqrt(1-((j./6.*e./2)./(aEllipse)).^2).*bEllipse;
    LE(:,:,i)=bAbs./aAbs.*(LE(:,:,i-10)-LE(:,:,i-20))+LE(:,:,i-20);
end
j=0;
for i=30:1:34
    j=j+1;
    bAbs=sqrt(1-((j./6.*e./2)./(aEllipse)).^2).*bEllipse;
    LE(:,:,i)=bAbs./aAbs.*(LE(:,:,i-10)-LE(:,:,i-20))+LE(:,:,i-20);
end

figure(figno)
hold on
mesh(Xup,Yup,Zup)
mesh(Xdown,Ydown,Zdown)
for i=25:1:34 
    plot3(LE(:,1,i),LE(:,2,i),LE(:,3,i),'o')
end
