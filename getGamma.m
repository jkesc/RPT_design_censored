function gamma=getGamma(beta1B,zla,D2m,D1,f)
     %following eqns in gulich table 3.2 page 132
    epsLim=exp(-8.16.*sin(beta1B)./zla);
    D2mred=D2m/D1;
    if D2mred>epsLim
        kw=1;
    else
        kw=1-((D2mred-epsLim)./(1-epsLim)).^3;
    end
    gamma=f.*(1-(sqrt(sin(beta1B)))./(zla.*0.7)).*kw;
end