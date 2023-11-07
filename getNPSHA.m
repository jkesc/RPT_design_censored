function NPSHA=getNPSHA(hBooster,p_atm,p_vapour,headDiffSubmergence,turbomachine)
    global rho g
    h_atm=p_atm./rho./g;
    h_vapour=p_vapour./rho./g;
    switch upper(turbomachine)
        case 'TURBINE'
            NPSHA=h_atm+hBooster+headDiffSubmergence-h_vapour;
        case 'PUMP'
            NPSHA=h_atm+hBooster+headDiffSubmergence-h_vapour;
        otherwise
            error('ERROR: no matchig value for variable "turbomachine"')
    end