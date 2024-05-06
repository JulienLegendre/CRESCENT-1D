function [eps,Eg] = epsAlGaAsSb_GonzalezCuevas(x,y,w,T)
%EPSALGAASSB_GONZALEZCUEVAS Computes the permittivity of Al(x)Ga(1-x)As(y)Sb(1-y), as a function of x and y
% data from Gonzalez-Cuevas, Journal of Applied Physics (2017)

    arguments
        x (1,1) double {mustBeInRange(x,0,1)} % fraction of Al (-)
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
    
    %AlAs
    dAlAs0=[3.09,3.37,3.98,4.18,4.86,2.46,5.6513,23.74,2.24,0.705,1e-3,6e-6,0.7421,0.49;...
            8.8e-4,8.8e-4,6.7e-4,6.7e-4,3.2e-4,6.1e-4,2.9e-5,0,0,0,0,2.9e-5,1.64e-4,0;...
            530,530,0,0,0,204,0,0,0,0,0,0,0,0];
    dAlAsT=dAlAs0(1,:);
    dAlAsT(1:6)=dAlAs0(1,1:6)-dAlAs0(2,1:6)*T^2./(T+dAlAs0(3,1:6));
    dAlAsT(7)=dAlAs0(1,7)+T*dAlAs0(2,7);
    dAlAsT(12:13)=dAlAs0(1,12:13)+T*dAlAs0(2,12:13);
    
    %GaSb
    dGaSb0=[0.81,1.57,2.19,2.62,4.32,0.88,6.086,1.102,3.34,4.93,1e-3,1e-6,0.6804,0.02;...
            4.2e-4,4.2e-4,6.8e-4,6.7e-4,9e-4,6e-4,4.72e-5,0,0,0,0,1.24e-4,4.8e-4,0;...
            140,140,147,176,376,140,0,0,0,0,0,0,0,0];
    dGaSbT=dGaSb0(1,:);
    dGaSbT(1:6)=dGaSb0(1,1:6)-dAlAs0(2,1:6)*T^2./(T+dGaSb0(3,1:6));
    dGaSbT(7)=dGaSb0(1,7)+T*dGaSb0(2,7);
    dGaSbT(12:13)=dGaSb0(1,12:13)+T*dGaSb0(2,12:13);
    
    %AlSb
    dAlSb0=[2.39,3.06,2.94,3.43,4.18,1.7,6.122,36.58,1.6,1.19,1e-3,1e-6,0.693,1.94;...
            4.2e-4,4.2e-4,4.7e-4,4.3e-4,4.7e-4,3.9e-4,2.6e-5,0,0,0,0,5.7e-5,3.06e-4,0;
            140,140,0,0,0,140,0,0,0,0,0,0,0,0];
    dAlSbT=dAlSb0(1,:);
    dAlSbT(1:6)=dAlSb0(1,1:6)-dAlSb0(2,1:6)*T^2./(T+dAlSb0(3,1:6));
    dAlSbT(7)=dAlSb0(1,7)+T*dAlSb0(2,7);
    dAlSbT(12:13)=dAlSb0(1,12:13)+T*dAlSb0(2,12:13);
    
    %% Section 2: Interpolation
    
    % GaSb
    if x==0 && y==0
        dAlGaAsSbT=dGaSbT;
    
    % GaAs
    elseif x==0 && y==1
        dAlGaAsSbT=dGaAsT;
       
    % AlSb
    elseif x==1 && y==0
        dAlGaAsSbT=dAlSbT;
        
    % AlAs
    elseif x==1 && y==1
        dAlGaAsSbT=dAlAsT;
        
    else
    
        % Interpolation : energy levels
        C=[0.37 0.07 0.45 0.00 0.02 0.06];
        dAlGaAsT=x*dAlAsT(1:6)+(1-x)*dGaAsT(1:6)-x*(1-x)*C;
        C=[0.72 0.15 0.00 0.00 0.00 0.28];
        dAlAsSbT=y*dAlAsT(1:6)+(1-y)*dAlSbT(1:6)-y*(1-y)*C;
        C=[1.20 0.61 0.00 0.00 0.00 1.09];
        dGaAsSbT=y*dGaAsT(1:6)+(1-y)*dGaSbT(1:6)-y*(1-y)*C;
        C=[0.69 0.30 0.28 0.32 0.00 0.55];
        dAlGaSbT=x*dAlSbT(1:6)+(1-x)*dGaSbT(1:6)-x*(1-x)*C;
        dAlGaAsSbT=1/(x*(1-x)+y*(1-y))*(x*(1-x)*(y*dAlGaAsT+(1-y)*dAlGaSbT)+...
                                        y*(1-y)*(x*dAlAsSbT+(1-x)*dGaAsSbT));

        % Interpolation : other parameters
        C=[0.002 0.012 0.019 0.010 -0.004 0.033 -0.031 0.008;...
           0.010 0.010 0.010 0.010  0.010 0.009  0.009 0.010];
        dAlGaAsSbT=[dAlGaAsSbT x*(y*dAlAsT(7:14)+(1-y)*dAlSbT(7:14))+...
                           (1-x)*(y*dGaAsT(7:14)+(1-y)*dGaSbT(7:14))+...
                          C(1,:)*x*(1-x)+C(2,:)*y*(1-y)];
                  
    end
    
    %% Section 3: Computation of the dielectric function
    
    E=hb*w/e;
    Eg=dAlGaAsSbT(1)*e;
    
    X0=(E+1i*dAlGaAsSbT(11))/dAlGaAsSbT(1);
    Xs0=(E+1i*dAlGaAsSbT(11))/dAlGaAsSbT(2);
    X1=(E+1i*dAlGaAsSbT(12))/dAlGaAsSbT(3);
    Xs1=(E+1i*dAlGaAsSbT(12))/dAlGaAsSbT(4);
    Xg=(E+1i*dAlGaAsSbT(14))/dAlGaAsSbT(6);
    fX0=X0.^(-2).*(2-sqrt(1+X0)-sqrt(1-X0));
    fXs0=Xs0.^(-2).*(2-sqrt(1+Xs0)-sqrt(1-Xs0));
    
    B1=44*(dAlGaAsSbT(3)+1/3*(dAlGaAsSbT(4)-dAlGaAsSbT(3)))/(dAlGaAsSbT(7)*dAlGaAsSbT(3)^2);
    B2=44*(dAlGaAsSbT(3)+2/3*(dAlGaAsSbT(4)-dAlGaAsSbT(3)))/(dAlGaAsSbT(7)*dAlGaAsSbT(4)^2);
    
    eps1=dAlGaAsSbT(8)*dAlGaAsSbT(1)^(-1.5)*(fX0+0.5*(dAlGaAsSbT(1)/dAlGaAsSbT(2))^1.5*fXs0);
    eps2=-B1.*X1.^(-2).*log(1-X1.^2)-B2.*Xs1.^(-2).*log(1-Xs1.^2);
    eps3=dAlGaAsSbT(9)*dAlGaAsSbT(5)^2./(dAlGaAsSbT(5)^2-E.^2-1i*E*dAlGaAsSbT(13));
    eps4=2*dAlGaAsSbT(10)/pi*(-Xg.^(-2)*log(dAlGaAsSbT(3)/dAlGaAsSbT(6))+...
                          1/2*(1+Xg.^(-1)).^2.*log((Xg+dAlGaAsSbT(3)/dAlGaAsSbT(6))./(Xg+1))+...
                          1/2*(1-Xg.^(-1)).^2.*log((Xg-dAlGaAsSbT(3)/dAlGaAsSbT(6))./(Xg-1)));
    eps=eps1+eps2+eps3+eps4;

end

