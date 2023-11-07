function NPSHR=getNPSHR(cm2,u2,Turbomachine)
global g
Turbomachine=upper(Turbomachine);
switch Turbomachine
    case 'PUMP'
        a=2;b=0.25;%Strictest
%         a=1.6;b=0.2;%Most relaxed
    case 'TURBINE'
        a=1.15;b=0.15;%Strictest
%         a=1.05;b=0.05;%Most relaxed
end
NPSHR=a.*cm2.^2/(2.*g)+b.*u2.^2/(2.*g);
end