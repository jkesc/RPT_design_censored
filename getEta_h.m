function eta_h=getEta_h(Q,QBEP,n,HBEP)
nq=n.*sqrt(QBEP)./(HBEP.^0.75);
m=0.08*0.5.*(1/QBEP).^0.15.*(45./nq).^0.06;
eta_hMax=1-0.055*(1/QBEP).^m-0.09.*(log(nq/45)).^2.5;       %Gülich page 142
qRed=Q/QBEP;
eta_h=eta_hMax.*(1-0.6.*(qRed-0.9).^2-0.25.*(qRed-0.9).^3); %as in Gülich page 166
end