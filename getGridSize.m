function deltaS=getGridSize(yPlusReq,rho,nu,L,U_inf)
Re=U_inf.*L./nu;
Cf=0.370.*(log10(Re)).^(-2.584);%Schultz-Grunov, https://www.cfd-online.com/Wiki/Skin_friction_coefficient
wallStress=0.5.*Cf.*rho.*U_inf.^2;
u_wall=sqrt(wallStress/rho);
deltaS=yPlusReq.*nu./u_wall;
end