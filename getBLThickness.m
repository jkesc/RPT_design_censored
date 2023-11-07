function delta=getBLThickness(U,x)
    global mu rho
    delta=5.5./sqrt(rho.*U.*x./mu).*x;
end