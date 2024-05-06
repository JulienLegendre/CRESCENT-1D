function d = Properties_AlGaAs(d,w)
%PROPERTIES_ALGAAS Gives back AlGaAs properties

    arguments
        d (1,1) struct                  % layer whose properties are to be computed
        w (:,1) double{mustBePositive}  % angular frequency (rad.s^-1)
    end
    
    global hb kb m0 e %#ok<GVMIS>
    
    % Semiconductor properties
    % data from Adachi, GaAs and related materials (1994), 
    % Levinshtein et al, Handbook series on semiconductor parameters, Vol. 2 (1999)
    d.epsr=12.9-2.84*d.x;                                                           % dielectric constant (static)
    epsinf=10.9-2.41*d.x;                                                           % dielectric constant (high frequency)
    d.mx=[0.063+0.083*d.x ((0.51+0.25*d.x)^1.5+(0.082+0.068*d.x)^1.5)^(2/3)]*m0;    % electron/hole effective mass (kg)
    d.chi=(4.07-1.1*d.x)*e;                                                         % electron affinity (eV)
    d.Ar=d.mx*kb^2*e/2/pi^2/hb^3;                                                   % Richardson constant
    
    % Electron/hole mobility
    % data from Sotoodeh et al, Journal of Applied Physics (2000)
    d.CTmodel.mu_min=[ 500*(0.063*m0/d.mx(1))^1.5*(1/10.89-1/12.9)/(1/epsinf-1/d.epsr)  57.1*d.x.^2-67.1*d.x+20];
    d.CTmodel.mu_max=[9400*(0.063*m0/d.mx(1))^1.5*(1/10.89-1/12.9)/(1/epsinf-1/d.epsr)  781*d.x.^2-1073*d.x+491.5];
    d.CTmodel.Nref=  [(6e16)^(1-d.x)*(5.46e17)^d.x                                    10^(1.40*d.x.^2-0.988*d.x+17.17)];
    d.CTmodel.alpha= [0.394+(1-0.394)*d.x                                             0.421*d.x.^2-0.313*d.x+0.38];
    d.CTmodel.theta1=((1-d.x)*[2.1 2.2]+d.x*[2.1 2.24])./(1+d.x.*(1-d.x));
    d.CTmodel.theta2=[3 3];
    mu=mobilityCTstruct(d.Ndop'/1e6,d.T,d.CTmodel)/1e4;
    d.D=kb*d.T/e*mu;
    
    % Material properties
    if ~isfield(d,'Eu')
        d.Eu=NaN;
    end
    [d.eps,d.Eg,d.iwm,d.fracInterband]=epsAlGaAs(d.x,w,d.T,d.Ndop,d.DopSign,d.mx,sum(mu'.*[1-d.DopSign;1+d.DopSign]/2),d.Eu);  %dielectric function (-)

    % Semiconductor computed properties
    d.Nc=2*(d.mx(1)*kb*d.T/(2*pi*hb^2))^1.5;   % effective density of states in the conduction band (m^-3)
    d.Nv=2*(d.mx(2)*kb*d.T/(2*pi*hb^2))^1.5;   % effective density of states in the valence band (m^-3)
    d.ni=sqrt(d.Nc*d.Nv)*exp(-d.Eg/(2*kb*d.T));% intrinsic carrier concentration (m^-3)

    % Electron/hole lifetime
    % data from Sadi et al, IEEE Transactions on Electron Devices (2019)
    d.tauP=1/(3e5);
    d.tauN=1/(3e5);
    d.Brad=2e-10*1e-6; % normally not used in the solver
    d.Cp=1e-30*1e-12;
    d.Cn=1e-30*1e-12;
    
    % saving material used
    if d.x==0
        d.mat="GaAs";
    elseif d.x==1
        d.mat="AlAs";
    else 
        d.mat="AlGaAs";
    end

end