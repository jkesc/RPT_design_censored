function beta1B=getBeta1B(gamma,c1m,tau1,D1m,n,beta1Bref,R)
%Code based on Johan Gülichs centrifugal pumps page 77
u1ref=pi*D1m*n/60;%Finding reference u at D1m
c1uRef=(gamma-(c1m.*tau1)/(u1ref.*tan(beta1Bref))).*u1ref;%finding flow velocity at D1m
I=c1uRef.*(D1m./2);%assuming free vortex (cu/r=const)
c1u=I./R(1,:);%finding cu1 at all of the inlet nodes
u1=2.*pi.*R(1,:).*n/60;%finding u1 at all inlet nodes
beta1B=atan((c1m.*tau1)./(u1.*gamma-c1u));%finding beta1B for all of the inlet nodes
end