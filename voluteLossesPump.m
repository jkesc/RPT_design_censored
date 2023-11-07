function H = voluteLossesPump(roughness,L,c,mu,rho,D,Q,u2,cx,cp,AR)
    Area=pi.*D.*L;
    zeta_spR=1/(Q.*u2.^2).*(sum((getFrictionFactor(roughness,L,c,mu,rho)+0.0015).*c.^3.*Area));
    zeta_spD=cx^2/u2^2.*(1-cp-1/AR^2);
    H=zeta_spD+zeta_spR;
end
