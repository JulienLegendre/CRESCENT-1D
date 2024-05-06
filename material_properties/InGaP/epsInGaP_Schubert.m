function [eps,Eg] = epsInGaP_Schubert(x,w,T)
%EPSINGAP_SCHUBERT Returns the interband contribution to dielectric
%function of InGaP.
% data from Schubert et al, Isotropic dielectric functions of highly
% disordered AlGaInP lattice-matched to GaAs (1999).

    arguments
        x (1,1) double {mustBeInRange(x,0,1)} % fraction of In (-)
        w (:,1) double {mustBePositive}       % angular frequency (rad.s^-1)
        T (1,1) double {mustBePositive}       % temperature (K)
    end

    global hb e %#ok<GVMIS>

    % extrapolation of the interband contribution for T!=300 K and x!=0.51
    % all the transition energy are supposed to vary as Eg
    % the temperature dependence is derived from the variation of Eg close
    % to room temperature for x=0.51.
    % the composition dependence is derived from the interpolation of Eg
    % for x=0, 0.51 and 0.71
    fEg=@(x,T) 1.34+0.8023*x+0.576*x^2-0.43e-3*(T-300);
    deltaE=fEg(x,T)-fEg(0.51,300);
    
    % Input data
    epsInf=0.52;
    E0=1.899+deltaE;
    Eg=E0*e;
    A0=10.44;
    G0=0.003;
    E1=3.224+deltaE;
    A1=5.27;
    G1=0.334;
    A1x=1.79;
    G1x=0.295;
    s1x=0.42;
    p1x=-0.76;
    E2=4.832+deltaE;
    A2=2.19;
    G2=0.743;
    s2=0.88;
    p2=0.16;

    % additional parameters
    E=hb*w/e;
    X0=(E+1i*G0)/E0;
    X1=(E+1i*G1)/E1;

    % calculation of the dielectric function
    eps0=A0*E0^(-1.5)*X0.^(-2).*(2-sqrt(1+X0)-sqrt(1-X0));
    eps1=-A1*X1.^(-2).*log(1-X1.^2);
    eps1x=A1x*exp(1i*p1x)./(E1-E-1i*G1x*exp(-s1x*(E1-E).^2));
    eps2=A2*exp(1i*p2)./(1-(E/E2).^2-1i*(E/E2)*(G2/E2).*exp(-s2*(E-E2).^2));
    eps=epsInf+eps0+eps1+eps1x+eps2;

end