function d = Properties_InGaAs(d,w)
%PROPERTIES_INGAAS Gives back In_{x}Ga_{1-x}As properties
    
    arguments
        d (1,1) struct                  % layer whose properties are to be computed
        w (:,1) double {mustBePositive} % angular frequency (rad.s^-1)
    end

    global hb kb m0 e %#ok<GVMIS>
    
    % Semiconductor properties
    % Levinshtein et al, Handbook series on semiconductor parameters, Vol. 2 (1999)
    d.epsr=12.9+1.53*d.x+0.67*d.x^2;                                                        % dielectric constant (static)
    d.mx=[0.063-0.043*d.x+0.003*d.x^2 ((0.51-0.1*d.x)^1.5+(0.082-0.056*d.x)^1.5)^(2/3)]*m0; % electron/hole effective mass (kg)
    d.chi=(4.07+0.83*d.x)*e;                                                                % electron affinity (eV)
    d.Ar=d.mx*kb^2*e/2/pi^2/hb^3;                                                           % Richardson constant
    
    % Electron/hole mobility
    % data from Sotoodeh et al, Journal of Applied Physics (2000)
    d.CTmodel.mu_min=[polyval(polyfit([0 0.53 1],[500 300 400],2),d.x)  polyval(polyfit([0 0.53 1],[20 10 20],2),d.x)];
    d.CTmodel.mu_max=[polyval(polyfit([0 0.53 1],[9400 14000 34000],2),d.x)  polyval(polyfit([0 0.53 1],[491.5 320 530],2),d.x)];
    d.CTmodel.Nref=  10.^[polyval(polyfit([0 0.53 1],log10([6e16 1.3e17 1.1e18]),2),d.x)  polyval(polyfit([0 0.53 1],log10([1.48e17 4.9e17 1.1e17]),2),d.x)];
    d.CTmodel.alpha= [polyval(polyfit([0 0.53 1],[0.394 0.48 0.32],2),d.x)  polyval(polyfit([0 0.53 1],[0.38 0.403 0.46],2),d.x)];
    d.CTmodel.theta1=[polyval(polyfit([0 0.53 1],[2.1 1.59 1.57],2),d.x)  polyval(polyfit([0 0.53 1],[2.2 1.59 2.3],2),d.x)];
    d.CTmodel.theta2=[polyval(polyfit([0 0.53 1],[3 3.68 3],2),d.x) 3];
    mu=mobilityCTstruct(d.Ndop'/1e6,d.T,d.CTmodel)/1e4;
    d.D=kb*d.T/e*mu;
    
    % Material properties
    if ~isfield(d,'Eu')
        d.Eu=NaN;
    end
    [d.eps,d.Eg,d.iwm,d.fracInterband]=epsInGaAs(d.x,w,d.T,d.Ndop,d.DopSign,d.mx,sum(mu'.*[1-d.DopSign;1+d.DopSign]/2),d.Eu);  %dielectric function (-)
    
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
        d.mat="InAs";
    else 
        d.mat="InGaAs";
    end

end