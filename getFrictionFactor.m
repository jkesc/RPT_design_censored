function cf=getFrictionFactor(roughness,L,c,mu,rho)
    REL=rho.*c.*L./mu;
    cf=0.136./((-log(0.2.*roughness./L+12.5./REL)).^2.15);
end