function deltaTheta=getDeltaTheta(R,deltaH)
R1=R(1:end-1,:);
R2=R(2:end,:);
deltaR=(deltaH.^2-R2.^2+R1.^2)./(2.*R1);
deltaTheta=acos((R1-deltaR)./(R2));
end