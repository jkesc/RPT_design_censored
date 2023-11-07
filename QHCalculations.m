% Comments
%CHECK OUT GÜLICH PAGE 140, THINK THE LOSSES ARE THERE.

%% 5. Impeller blade number (This one should be sat outside of the switch statement Also must find out of number of guide vanes)        
        %zla number of blades
        %zle number of vanes
        resonance=true;
        zle=Values.zle;
        zla=Values.zla;
        nu1=[1,2,3].*[1;1;1];
        nu2=nu1';
        m=zeros(3,3);
        while resonance==true
            m=abs(nu2.*zla-nu1.*zle);
            if sum(sum((m==0)+(m==1)))==0
                resonance=false;
            else
                zla=zla-1;
            end
        end
%% Initializing variables
switch FindBetaBlade
    case 'freeVortex'
        D2m=sqrt(0.5*((2*R(iEllipse,jEllipse))^2+(2*R(iEllipse,1))^2));%Mean diameter at 2 (Suction side)
        D1=(R(1,jEllipse)+R(1,1));%Assuming D1 is the arithmetic average value.=sqrt(0.5*((2*R(1,jEllipse))^2+(2*R(1,1))^2));        %assuming that D1 ia supposed to be the geometric average as well.
        u=(pi.*D1.*n)./(60);                            %tangential velocity at 1
        switch nondimensional
            case 'yes'
                Values.DimensionlessH=u^2/g;%u_1^2/2g (head coefficient)
        end
        A2=pi.*(R(iEllipse,jEllipse).^2-R(iEllipse,1).^2);                      %Area at 2
        Hmax=(u).^2./g;                                 %Highest theoretical shutoff height
        lambdaLa=pi./2;                                 %Angle between blades and plates
        alpha2=pi./2;                                   %Angle between circumferential direction and absolute velocity
        %Shuld probably make this like in gulich page 133, eqn 3.3.7 
        %Hfun=@(Q,beta,u,A,eta_h) (u.^2./g-Q.*u./(A.*g.*tan(beta))).*eta_h;
        %PS the last term divides by tan(alpha), so if alpha=90deg(no prerotation),
        %the last term becomes 0.
        Hfun=@(Q,beta1B,u,A1,A2,eta_h,gamma,tau_1,D2m,D1,alpha2,f) (eta_h*u.^2)/(g).*(gamma-(Q)./(f*A1*u*tan(beta1B)).*(tau_1+(A1*(D2m./D1)*tan(beta1B))/(A2*tan(alpha2))));%p.133 Gülich
        Hideal=@(u,Q,A,beta) u./g.*(u-Q./(A.*tan(beta)));
        hf=Values.hf;                         %Friction losses of the existing system
        hfDraftTube=Values.hfDraftTube;                %FrictionLosses from draft tube to the lower reservoir.
        hBooster=Values.hBooster;%Assuming ideal slipless characteristic for boosterpump, working in the same range of Q as the RPT.
        %Beta with slip:
        beta1B=deg2rad(10);                             %Setting initial value for blade angle
        deltaBeta=1;
        bladeBlockage= @(beta1B,lambdaLa) 1./(1-(zla.*e)./(pi.*D1.*sin(beta1B).*sin(lambdaLa)));%Gülich page XXX
        while true                                          %Lifehack for "do while" loop in matlab
            beta1B=beta1B+deg2rad(deltaBeta);
            f=Values.f;                                         %Different for semi axial turbine. not quite sure wether the pump is radial or semi-axial (It's not PURELY radial at least).
            cm1=Q./A;
            c_1uB=u-(cm1)./(tan(beta1B));        
        switch empiricalGamma
            case 'yes'
                cu1CFD=Values.cu1CFD(1);%Not sure how I'll extract the values from CFX yet...
                gamma=getGammaEmpirical(u,cu1CFD,c_1uB);
            case 'no'
                gamma=Values.gamma;%frequently between 0.7 and 0.8
            otherwise
                error('ERROR: no valid value for variable EmpiricalGamma')
        end
            beta1Slip=atan(cm1/(u-(c_1uB-(1-gamma).*u)));
            tau_1=bladeBlockage(beta1B,lambdaLa);
            eta_h=getEta_h(Q,Q,n,Head);%finding the hydraulic efficiency at Q=QBEP. therefore same Q twice.
            %eta_h=.94;
            %if (Head+hf(Q))< Hideal(u,Q,A,beta1B)
            Hreq=Hfun(Q,beta1B,u,A,A2,eta_h,gamma,tau_1,D2m,D1,alpha2,f);
            if (Head+hf(Q)-hBooster(Q)) <Hreq
                break;
            end
        end
        % betaDim=atan((Q)./(A.*(u.^2-Head.*g)));            %Dimensioning beta1'(So with slip). This is however still without losses. Should rather make QH-curves first, and then decide.
        Qfun=@(Hfun) tan(beta1Slip).*A.*(u.^2-Hfun.*g);     %function to find Q
        % DeltaQ=Qfun(53)-Qfun(Head);                        %Difference between Q at highest and lowest difference, no losses accounted for
        % RelativeIncreaseOfQ=DeltaQ./Qfun(Head);            %Relative increase of Q
        % Qmax=Qfun(0);

        steps=97;%number of Q-values to plot the head for.
        figure(figno)%plotting the pump characteristic 
        figname.PumpChar=figno;
        figno=figno+1;
        hold on
        yyaxis left
        plot(Hfun(0:1:steps,beta1B,u,A,A2,getEta_h(0:1:steps,Q,n,Head),gamma,tau_1,D2m,D1,alpha2,f)./Values.DimensionlessH); %getGamma(beta1B,zla,D2m,D1,f),tau_1,D2m,D1,alpha2,f));           %beta1B,u,A,getEta_h(0:1:105,Q,n,Head)));
        %plot(hBooster(0:1:steps)+Hfun(0:1:steps,beta1B,u,A,A2,getEta_h(0:1:steps,Q,n,Head),getGamma(beta1B,zla,D2m,D1,f),tau_1,D2m,D1,alpha2,f));           %beta1B,u,A,getEta_h(0:1:105,Q,n,Head)));

        %plot(Hfun(0:1:steps,beta1B,u,A,A2,eta_h,1,tau_1,D2m,D1,alpha2,f));
        % legend('Head with losses')
        plot((Head+hf(0:1:steps)-hBooster(0:1:steps))./Values.DimensionlessH)
        plot((Hideal(u,0:1:steps,A,beta1B))./Values.DimensionlessH)
        plot((Hideal(u,0:1:steps,A,beta1Slip))./Values.DimensionlessH)
        % legend(['System curve at H=' Head])
        switch nondimensional
            case 'yes'
                ylabel('\Psi')
            case 'no'
                ylabel('Head')
        end
        yyaxis right
        plot(getEta_h(0:1:steps,Q,n,Head),'-x')
        legend('Head with losses',sprintf('System curve at H=%.2f m',Head./Values.DimensionlessH),'characteristic without losses','lossless characteristic with slip',sprintf('Efficiency for Q_{BEP}=%.2fm^3/s',Q),'Location','southwest')
        ylabel('Efficiency')
        xlabel('Q')
        %title('Pump characteristic')
        betaDim=rad2deg(beta1B);
    case 'streamlines'
        beta1=zeros(1,jEllipse);
        Qold=Q;
        for j=1:1:jEllipse
            D2m=2*R(iEllipse,j);% decided to use the raduius of the streamlines as reference, and to use average data from half the interval above nde half below for calculations of volume flow rate etc.
            D1=2*R(1,j);%Assuming D1 is the arithmetic average value.=sqrt(0.5*((2*R(1,jEllipse))^2+(2*R(1,1))^2));        %assuming that D1 ia supposed to be the geometric average as well.
            u=(pi.*D1.*n)./(60);                            %tangential velocity at 1
            if j==1
                A2=0.5.*(getdAij((G(iEllipse,j)+GToAdd(:,j)),dA1,dA2,P.Aij,GGeom(end,j)));% using only half of the area for first and last line
                A1=0.5.*(getdAij((G(1,j)+GToAdd(:,j)),dA1,dA2,P.Aij,GGeom(end,j)));
                Q=Qold./(jEllipse-1)./2;
            elseif j==jEllipse
                A2=0.5.*(getdAij((G(iEllipse,j-1)+GToAdd(:,j-1)),dA1,dA2,P.Aij,GGeom(end,j-1)));
                A1=0.5.*(getdAij((G(1,j-1)+GToAdd(:,j-1)),dA1,dA2,P.Aij,GGeom(end,j-1)));
                Q=Qold./(jEllipse-1)./2;
            else
                A2=0.5.*(getdAij((G(iEllipse,j)+GToAdd(:,j)),dA1,dA2,P.Aij,GGeom(end,j))+getdAij((G(iEllipse,j-1)+GToAdd(:,j-1)),dA1,dA2,P.Aij,GGeom(end,j-1)));%As the code uses dAij to find the area between streamline j and j+1, the area for each streamline is assumed to be the average of the area for j and j+1
                A1=0.5.*(getdAij((G(1,j)+GToAdd(:,j)),dA1,dA2,P.Aij,GGeom(end,j))+getdAij((G(1,j-1)+GToAdd(:,j-1)),dA1,dA2,P.Aij,GGeom(end,j-1)));
                Q=Qold./(jEllipse-1);
            end
            Hmax=(u).^2./g;                                 %Highest theoretical shutoff height
            lambdaLa=pi./2;                                 %Angle between blades and plates
            alpha2=pi./2;                                   %Angle between circumferential direction and absolute velocity
            %Shuld probably make this like in gulich page 133, eqn 3.3.7 
            %Hfun=@(Q,beta,u,A1,eta_h) (u.^2./g-Q.*u./(A1.*g.*tan(beta))).*eta_h;
            %PS the last term divides by tan(alpha), so if alpha=90deg(no prerotation),
            %the last term becomes 0.
            Hfun=@(Q,beta1B,u,A1,A2,eta_h,gamma,tau_1,D2m,D1,alpha2,f) (eta_h*u.^2)/(g).*(gamma-(Q)./(f*A1*u*tan(beta1B)).*(tau_1+(A1*(D2m./D1)*tan(beta1B))/(A2*tan(alpha2))));%p.133 Gülich
            Hideal=@(u,Q,A1,beta) u./g.*(u-Q./(A1.*tan(beta)));
            hf=Values.hf;                         %Friction losses of the existing system
            hfDraftTube=Values.hfDraftTube;                %FrictionLosses from draft tube to the lower reservoir.
            hBooster=Values.hBooster; %Assuming ideal slipless characteristic for boosterpump, working in the same range of Q as the RPT.
            %Beta with slip:
            beta1B=deg2rad(10);                             %Setting initial value for blade angle
            deltaBeta=1;
            bladeBlockage= @(beta1B,lambdaLa) 1./(1-(zla.*e)./(pi.*D1.*sin(beta1B).*sin(lambdaLa)));%Gülich page XXX
            while true                                          %Lifehack for "do while" loop in matlab
                beta1B=beta1B+deg2rad(deltaBeta);
                f=Values.f;                                         %Different for semi axial turbine. not quite sure wether the pump is radial or semi-axial (It's not PURELY radial at least).
            %    gamma=getGamma(beta1B,zla,D2m,D1,f);
                    cm1=Q./A1;
                c_1uB=u-(cm1)./(tan(beta1B));
                 switch empiricalGamma
                    case 'yes'
                        %cu1CFD=cu1CFD(1);%Not sure how I'll extract the values from CFX yet...
                        gamma=getGammaEmpirical(u,cu1CFD(j),c_1uB);
                    case 'no'
                        gamma=Values.gamma;
                     otherwise
                         error('ERROR: no matching values for variable empiricalGamma')
                 end
                beta1Slip=atan(cm1/(u-(c_1uB-(1-gamma).*u)));
                tau_1=bladeBlockage(beta1B,lambdaLa);
                eta_h=getEta_h(Q,Q,n,Head);%finding the hydraulic efficiency at Q=QBEP. therefore same Q twice.
                %eta_h=.94;
                %if (Head+hf(Q))< Hideal(u,Q,A1,beta1B)
                Hreq=Hfun(Q,beta1B,u,A1,A2,eta_h,gamma,tau_1,D2m,D1,alpha2,f);
                if (Head+hf(Q)-hBooster(Q)) <Hreq
                    break;
                end
            end
            % betaDim=atan((Q)./(A1.*(u.^2-Head.*g)));            %Dimensioning beta1'(So with slip). This is however still without losses. Should rather make QH-curves first, and then decide.
            % DeltaQ=Qfun(53)-Qfun(Head);                        %Difference between Q at highest and lowest difference, no losses accounted for
            % RelativeIncreaseOfQ=DeltaQ./Qfun(Head);            %Relative increase of Q
            Qmax=Q*1.5;
            Qval=linspace(0,Qmax);
            steps=97;%number of Q-values to plot the head for.
            figure(figno)%plotting the pump characteristic 
            figname.PumpChar=figno;
            figno=figno+1;
            hold on
            yyaxis left
            plot(Qval./Values.DimensionlessQ,Hfun(Qval,beta1B,u,A1,A2,eta_h,gamma,tau_1,D2m,D1,alpha2,f)./Values.DimensionlessH); %getGamma(beta1B,zla,D2m,D1,f),tau_1,D2m,D1,alpha2,f));           %beta1B,u,A1,getEta_h(0:1:105,Q,n,Head)));
            %plot(hBooster(0:1:steps)+Hfun(0:1:steps,beta1B,u,A1,A2,getEta_h(0:1:steps,Q,n,Head),getGamma(beta1B,zla,D2m,D1,f),tau_1,D2m,D1,alpha2,f));           %beta1B,u,A1,getEta_h(0:1:105,Q,n,Head)));

            %plot(Hfun(0:1:steps,beta1B,u,A1,A2,eta_h,1,tau_1,D2m,D1,alpha2,f));
            % legend('Head with losses')
            plot(Qval./Values.DimensionlessQ,(Head+hf(Qval)-hBooster(Qval))./Values.DimensionlessH)
            plot(Qval./Values.DimensionlessQ,Hideal(u,Qval,A1,beta1B)./Values.DimensionlessH)
            plot(Qval./Values.DimensionlessQ,Hideal(u,Qval,A1,beta1Slip)./Values.DimensionlessH)
            legend('with slip','system curve','\eta=1 without slip','\eta=1 with slip')
            % legend(['System curve at H=' Head])
            ylabel('Head')
            yyaxis right
            %plot(getEta_h(Qval,Q,n,Head),'-x')
            %legend('Head with losses',sprintf('System curve at H=%.2f m',Head),'characteristic without losses','lossless characteristic with slip',sprintf('Efficiency for Q_{BEP}=%.2fm^3/s',Q),'Location','southwest')
            ylabel('Efficiency')
            xlabel('Q')
            %title('Pump characteristic')
            betaDim=rad2deg(beta1B);
            beta1(j)=beta1B;
        end
        Q=Qold;
    otherwise
        error('ERROR, no matching value for variable FindBetaBlade')
end