%Diffusor losses pump
function H = DiffuserLossesPump(c2u,c2m,Q,a3,b3,a4,b4,LGuideVane,u2,zle,g)
    c2=sqrt(c2u.^2+c2m.^2);
    c3q=Q./zle./a3./b3;
    AR=a4*b4/a3/b3;
    fprintf('Read the value of c_p for L/h=%f and A_R-1=%f in Gülich figure 1.18',LGuideVane.*sqrt(pi./a3./b3),AR-1);
    cp=0.4;%input('\nc_p: ');
    zetaLe=(c3q./u2).^2.*(0.3.*(c2./c3q-1).^2+1-cp-(1)./(AR.^2));
    H=zetaLe.*u2.^2./2./g;
end
