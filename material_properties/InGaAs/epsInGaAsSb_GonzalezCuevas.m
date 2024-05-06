function [eps,Eg] = epsInGaAsSb_GonzalezCuevas(x,y,w,T)
%EPSALGAASSB_GONZALEZCUEVAS Computes the permittivity of In(x)Ga(1-x)As(y)Sb(1-y), as a function of x and y
% data from Gonzalez-Cuevas, Journal of Applied Physics (2017)

    arguments
        x (1,1) double {mustBeInRange(x,0,1)} % fraction of In (-)
        y (1,1) double {mustBeInRange(y,0,1)} % fraction of As (-)
        w (:,1) double {mustBePositive}       % angular frequency (rad.s^-1)
        T (1,1) double {mustBePositive}       % temperature (K)
    end

    global hb e %#ok<GVMIS>
    
    %% Section 1: Loading data
    
    %GaAs
    dGaAs0=[1.52,1.85,3.04,3.27,5.13,1.82,5.638,5.52,2.89,21.32,1e-3,1e-6,0.6631,0.037;...
            5.5e-4,3.5e-4,7.2e-4,7.2e-4,6.6e-4,6.1e-4,3.88e-5,0,0,0,0,6.6e-5,4.1e-4,0;...
            225,225,205,205,43,204,0,0,0,0,0,0,0,0];
    dGaAsT=dGaAs0(1,:);
    dGaAsT(1:6)=dGaAs0(1,1:6)-dGaAs0(2,1:6)*T^2./(T+dGaAs0(3,1:6));
    dGaAsT(7)=dGaAs0(1,7)+T*dGaAs0(2,7);
    dGaAsT(12:13)=dGaAs0(1,12:13)+T*dGaAs0(2,12:13);
    
    %InAs
    dInAs0=[0.42,0.79,2.61,2.88,4.74,1.13,6.0518,19.95,1.75,21.511,8.66,0.25e-3,0.403,0.35;...
            2.8e-4,3.4e-4,5e-4,5e-4,5.6e-4,2.8e-4,2.74e-5,0,0,0,0,3e-4,3e-4,0;...
            93,248,0,0,0,93,0,0,0,0,0,0,0,0];
    dInAsT=dInAs0(1,:);
    dInAsT(1:6)=dInAs0(1,1:6)-dInAs0(2,1:6)*T^2./(T+dInAs0(3,1:6));
    dInAsT(7)=dInAs0(1,7)+T*dInAs0(2,7);
    dInAsT(12:13)=dInAs0(1,12:13)+T*dInAs0(2,12:13);
    
    %GaSb
    dGaSb0=[0.81,1.57,2.19,2.62,4.32,0.88,6.086,1.102,3.34,4.93,1e-3,1e-6,0.6804,0.02;...
            4.2e-4,4.2e-4,6.8e-4,6.7e-4,9e-4,6e-4,4.72e-5,0,0,0,0,1.24e-4,4.8e-4,0;...
            140,140,147,176,376,140,0,0,0,0,0,0,0,0];
    dGaSbT=dGaSb0(1,:);
    dGaSbT(1:6)=dGaSb0(1,1:6)-dInAs0(2,1:6)*T^2./(T+dGaSb0(3,1:6));
    dGaSbT(7)=dGaSb0(1,7)+T*dGaSb0(2,7);
    dGaSbT(12:13)=dGaSb0(1,12:13)+T*dGaSb0(2,12:13);
    
    %InSb
    dInSb0=[0.24,1.2,2,2.49,4.24,0.93,6.4696,16.83,2.64,38.83,4.2,1.5e-4,0.6339,0.56;...
            3.2e-4,3.2e-4,6.8e-4,6.5e-4,5.4e-4,3.2e-4,3.48e-5,0,0,0,0,2.8e-4,4.9e-4,0;
            170,170,132,170,0,170,0,0,0,0,0,0,0,0];
    dInSbT=dInSb0(1,:);
    dInSbT(1:6)=dInSb0(1,1:6)-dInSb0(2,1:6)*T^2./(T+dInSb0(3,1:6));
    dInSbT(7)=dInSb0(1,7)+T*dInSb0(2,7);
    dInSbT(12:13)=dInSb0(1,12:13)+T*dInSb0(2,12:13);
    %dInSbT(9)=dInSbT(9)*2; % correction temporaire pour mieux matcher
    %dInSbT(13)=dInSbT(13)*1.6;
    
    %% Section 2: Interpolation
    
    % GaSb
    if x==0 && y==0
        dInGaAsSbT=dGaSbT;
    
    % GaAs
    elseif x==0 && y==1
        dInGaAsSbT=dGaAsT;
       
    % InSb
    elseif x==1 && y==0
        dInGaAsSbT=dInSbT;
        
    % InAs
    elseif x==1 && y==1
        dInGaAsSbT=dInAsT;
        
    else
    
        % Interpolation : energy levels
        C=[0.61 0.20 0.50 0.49 0.27 0.72];
        dInGaAsT=x*dInAsT(1:6)+(1-x)*dGaAsT(1:6)-x*(1-x)*C;
        C=[0.58 1.20 0.55 0.40 0.60 0.57];
        dInAsSbT=y*dInAsT(1:6)+(1-y)*dInSbT(1:6)-y*(1-y)*C;
        C=[1.20 0.61 0.00 0.00 0.00 1.09];
        dGaAsSbT=y*dGaAsT(1:6)+(1-y)*dGaSbT(1:6)-y*(1-y)*C;
        C=[0.42 0.10 0.38 0.28 0.24 0.38];
        dInGaSbT=x*dInSbT(1:6)+(1-x)*dGaSbT(1:6)-x*(1-x)*C;
        dInGaAsSbT=1/(x*(1-x)+y*(1-y))*(x*(1-x)*(y*dInGaAsT+(1-y)*dInGaSbT)+...
                                        y*(1-y)*(x*dInAsSbT+(1-x)*dGaAsSbT));

        % Interpolation : other parameters
        C=[-2.813 12.389 -1.648 17.938 -5.017 0.231 0.103 -0.079;...
           -0.016 0.001 0.002 0.001 -0.002 -0.046 0.058 -0.014];
        dInGaAsSbT=[dInGaAsSbT x*(y*dInAsT(7:14)+(1-y)*dInSbT(7:14))+...
                           (1-x)*(y*dGaAsT(7:14)+(1-y)*dGaSbT(7:14))+...
                          C(1,:)*x*(1-x)+C(2,:)*y*(1-y)];
                  
    end
    
    %% Section 3: Computation of the dielectric function
    
    E=hb*w/e;
    Eg=dInGaAsSbT(1)*e;
    
    X0=(E+1i*dInGaAsSbT(11))/dInGaAsSbT(1);
    Xs0=(E+1i*dInGaAsSbT(11))/dInGaAsSbT(2);
    X1=(E+1i*dInGaAsSbT(12))/dInGaAsSbT(3);
    Xs1=(E+1i*dInGaAsSbT(12))/dInGaAsSbT(4);
    X2=(E+1i*dInGaAsSbT(13))/dInGaAsSbT(5);
    Xg=(E+1i*dInGaAsSbT(14))/dInGaAsSbT(6);
    fX0=X0.^(-2).*(2-sqrt(1+X0)-sqrt(1-X0));
    fXs0=Xs0.^(-2).*(2-sqrt(1+Xs0)-sqrt(1-Xs0));
    
    B1=44*(dInGaAsSbT(3)+1/3*(dInGaAsSbT(4)-dInGaAsSbT(3)))/(dInGaAsSbT(7)*dInGaAsSbT(3)^2);
    B2=44*(dInGaAsSbT(3)+2/3*(dInGaAsSbT(4)-dInGaAsSbT(3)))/(dInGaAsSbT(7)*dInGaAsSbT(4)^2);
    
    eps1=dInGaAsSbT(8)*dInGaAsSbT(1)^(-1.5)*(fX0+0.5*(dInGaAsSbT(1)/dInGaAsSbT(2))^1.5*fXs0);
    eps2=-B1.*X1.^(-2).*log(1-X1.^2)-B2.*Xs1.^(-2).*log(1-Xs1.^2);
    eps3=dInGaAsSbT(9)./(1-real(X2).*X2);
    eps4=2*dInGaAsSbT(10)/pi*(-Xg.^(-2)*log(dInGaAsSbT(3)/dInGaAsSbT(6))+...
                          1/2*(1+Xg.^(-1)).^2.*log((Xg+dInGaAsSbT(3)/dInGaAsSbT(6))./(Xg+1))+...
                          1/2*(1-Xg.^(-1)).^2.*log((Xg-dInGaAsSbT(3)/dInGaAsSbT(6))./(Xg-1)));
    eps=eps1+eps2+eps3+eps4;

end

