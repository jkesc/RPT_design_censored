%Values
global Values g rho mu nu patm pvapour
Values.HeadDiffSubmergence=;      %The height difference between the average lower reservoir height and the RPT
Values.R1GuideVanes=;                              %Radius of the guideVanes
Values.D1=;                        %Current diameter at the guide vanes
Values.Head=;                                   %The design head
Values.Q=;                                           %From efficiencies measured in 2005, Q_BEP is approximately 60. Assume this to be the case for the new one as well.
Values.b=;
Values.Dhub=;                                     %(Quite arbitrary chosen) hub diameter
Values.D2=;                                       %Diameter imposed by existing geometry
Values.n= ;                                         %current rotational velocity increased by one syncronous vel.
Values.iEllipse=15;
Values.jEllipse=15;
Values.runnerheight=; %From draft tube to top of guide vane.
Values.Routlet(1)=;                                                         %Values from the supplied figure.
Values.Zoutlet(1)=;
Values.Routlet(2)=;
Values.Zoutlet(2)=;
Values.RLe1=;                                                                  %The target values for The leading edge
Values.ZLe1=;
Values.RLe2=;
Values.ZLe2=;
Values.P2=[    Values.runnerheight];%ginput();%[0,runnerheight];
Values.ax=0.25;                                                            %initial x value to define the acceleration along the blade
Values.ay=0.75;                                                             %same for y
Values.dx=0.01;                                                            %Value to change ax and ay with for each iteration
Values.dy=0.01;
Values.P.Aij=[
        
    ];
switch Geometry
	case 'BezierDefByPoints'
	Values.P.LEp=[                                                                     %P for points
        
    ];%The points which generate the Leading edge by use of Bezier curves
	Values.P.TEpmid=[
    %1.2406    0.8690
    %1.2639    0.5361
    ];
	otherwise
	Values.P.LEp=[
            
    ];%The points which generate the Leading edge by use of Bezier curves
	Values.P.TEp=[ 
            1.4071    0.0910
            0.7877    0.2961
            0.3077    0.6910
        %     1.5539    0.4903
        %     1.2406    0.8690
        %     1.2223    1.3040
            ];
end
Values.P.TERunner=[%I don't recall making this...
   
    ];
Values.P.LERunner=[
   
    ];
Values.zle=;
Values.zla=8;      %Has to be smaller or equal to 8.
Values.hf=@(Q) *Q.^2;
Values.hfDraftTube=@(Q) *Q.^2;                %FrictionLosses from draft tube to the lower reservoir.
Values.hBooster=@(Q) ;%Assuming ideal slipless characteristic for boosterpump, working in the same range of Q as the RPT.

Values.f=0.98;                                         %Different for semi axial turbine. not quite sure wether the pump is radial or semi-axial (It's not PURELY radial at least).
Values.gamma=0.8;
Values.P.UcU=[0 1; 0.3 1 ;0.5 0; 1 0];
Values.P.BetaDist=[0 1;0.7 0.8;0.7 0.1;1 0];%[0 1;0.1 0.1; 1 0];
Values.P.Spanwise=[0 1;0.9 1;1 2];%To control how the blade changes in the spanwise direction
Values.RGuidevanes=[0.5,1,1.5,2.0,2.5];
%Values for cu1 from CFD:
Values.cu1CFD=[];
Values.TurbulentEnergyCorrectionFactor=1.075;%Energy correction factor for turbulent flow (approximate)
Values.LaminarEnergyCorrectionFactor=2;%Energy correction factor for laminar flow
switch nondimensional
    case 'yes'
        Values.DimensionlessRadius=Values.R1GuideVanes;%Value for nondimenzionalization of lengths
        Values.DimensionlessVelocity=Values.DimensionlessRadius*Values.n/60*2*pi;%Value for nondimenzionalization of velocities
        Values.DimensionlessQ=Values.Q;%DesignFlowRate
        Values.DimensionlessH=(Values.R1GuideVanes*Values.n/60*pi*2)^2/g;%u_1^2/2g (head coefficient)
    case 'no'
        Values.DimensionlessRadius=1;%Value for nondimenzionalization of lengths
        Values.DimensionlessVelocity=1;%Value for nondimenzionalization of velocities
        Values.DimensionlessQ=1;%DesignFlowRate
        Values.DimensionlessH=1;
    otherwise
        error('ERROR: no matching value for variable "nondimensional". Must be either yes or no.')
end