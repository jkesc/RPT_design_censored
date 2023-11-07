function gamma= getGammaEmpirical(u1,cu1,cu1_inf)
    gamma=1-(cu1_inf-cu1)./(u1);
end