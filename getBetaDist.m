function betaAlongBlade=getBetaDist(beta1B,beta2B,iEllipse,jEllipse,Fraction,option,P,PSpanwise)
if nargin==7
    PSpanwise=[0 1;1 1];
end
switch option
    case 'Node'
        betaAlongBlade=zeros(iEllipse,jEllipse);
    case 'lineSegment'
        betaAlongBlade=zeros(iEllipse-1,jEllipse);
    otherwise
        error('ERROR: "Option" not valid')
end
%                 for j=1:1:jEllipse            
%             [p]=[0 beta1B(j);0.7 1.0.*beta1B(j);0.7 beta2B(j);1 beta2B(j)];%input('Enter co-odinates of points within brackets ->[x1 y1;x2 y2;x3 y3;...;xn yn] ');
%             n=size(p,1);
%              n1=n-1;
%             sigma=zeros(1,n);
%                 i=0:1:n1;
%         sigma(i+1)=factorial(n1)./(factorial(i).*factorial(n1-i));
%     %         for    i=0:1:n1
%     %             sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  % for calculating (x!/(y!(x-y)!)) values 
%     %         end
%             l=[];
%             UB=zeros(1,n);
%             for u=0:0.002:1
%                 for d=1:n
%                     UB(d)=sigma(d)*((1-u)^(n-d))*(u^(d-1));
%                 end
%                 l=cat(1,l,UB);                                      %catenation 
%             end
%             P=l*p;
%     %        figure(8)
%     %         line(P(:,1),P(:,2))
%     %         line(p(:,1),p(:,2))
%         foo=interparc(iEllipse,P(:,1),P(:,2));
%         betaAlongBlade(:,j)=foo(:,2);
%         end


TransversalChangeFactor=ppval(extractBezierByPoints(0,PSpanwise),[0:1:jEllipse-1]./(jEllipse-1));%To change the blade by a factor along its span
Pold=P;
for j=1:1:jEllipse
    P(2:end-1,2)=Pold(2:end-1,2).*TransversalChangeFactor(j);
    betaAlongBlade(:,j)=ppval(extractBezierByPoints(0,P),Fraction(:,j)).*(beta1B(j)-beta2B(j))+beta2B(j);
    %betaAlongBlade(2:end-1,j)=betaAlongBlade(2:end-1,j).*TransversalChangeFactor(j);
end
betaAlongBlade(1,:)=beta1B;
betaAlongBlade(end,:)=beta2B;
end