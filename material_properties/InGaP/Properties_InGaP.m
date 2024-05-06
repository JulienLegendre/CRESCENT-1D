function d = Properties_InGaP(d,w)
%PROPERTIES_INGAP Gives back InGaP properties, x being the proportion of Ga
    
    arguments
        d (1,1) struct                  % layer whose properties are to be computed
        w (:,1) double {mustBePositive} % angular frequency (rad.s^-1)
    end

    global hb kb m0 e %#ok<GVMIS>
    
    % Semiconductor properties
    % Levinshtein et al, Handbook series on semiconductor parameters, Vol. 2 (1999)
    d.epsr=12.5-1.4*d.x;                                                            % dielectric constant (static)
    d.mx=[0.08+0.008*d.x/0.51 ((0.6+0.19*d.x)^1.5+(0.09+0.05*d.x)^1.5)^(2/3)]*m0;   % electron/hole effective mass (kg)
    d.chi=(4.38-0.58*d.x)*e;                                                        % electron affinity (J)
    d.Ar=d.mx*kb^2*e/2/pi^2/hb^3;                                                   % Richardson constant
    
    % Electron/hole mobility
    % data from Sotoodeh et al, Journal of Applied Physics (2000)
    d.CTmodel.mu_min=[ 400*(1-d.x/0.51)+ 400*d.x/0.51       sum(polyfit([0 0.51 1],[10 15 10],2).*[d.x^2 d.x 1])];
    d.CTmodel.mu_max=[5200*(1-d.x/0.51)+4300*d.x/0.51       sum(polyfit([0 0.51 1],[170 150 147],2).*[d.x^2 d.x 1])];
    d.CTmodel.Nref=  [10^sum(polyfit([0 0.51 1],[log10(3e17) log10(2e16) log10(4.4e18)],2).*[d.x^2 d.x 1])  10^sum(polyfit([0 0.51 1],[log10(4.87e17) log10(1.5e17) log10(1e18)],2).*[d.x^2 d.x 1])];
    d.CTmodel.alpha= [sum(polyfit([0 0.51 1],[0.47 0.7 0.8],2).*[d.x^2 d.x 1])                              sum(polyfit([0 0.51 1],[0.62 0.8 0.85],2).*[d.x^2 d.x 1])];
    d.CTmodel.theta1=[sum(polyfit([0 0.51 1],[2 1.66 1.6],2).*[d.x^2 d.x 1]) 2];
    d.CTmodel.theta2=(1-d.x)*[3.25 3]+d.x*[0.71 0];
    mu=mobilityCTstruct(d.Ndop'/1e6,d.T,d.CTmodel)/1e4;
    d.D=kb*d.T/e*mu;
    
    % Material properties
    [d.eps,d.Eg,d.iwm,d.fracInterband]=epsInGaP(d.x,w,d.T,d.Ndop,d.DopSign,d.mx,sum(mu'.*[1-d.DopSign;1+d.DopSign]/2),d.Eu);  %dielectric function (-)

    % Semiconductor computed properties
    d.Nc=2*(d.mx(1)*kb*d.T/(2*pi*hb^2))^1.5;   % effective density of states in the conduction band (m^-3)
    d.Nv=2*(d.mx(2)*kb*d.T/(2*pi*hb^2))^1.5;   % effective density of states in the valence band (m^-3)
    d.ni=sqrt(d.Nc*d.Nv)*exp(-d.Eg/(2*kb*d.T));% intrinsic carrier concentration (m^-3)

    % Electron/hole lifetime - tau and Brad are GaAs data
    % data from Sadi et al, IEEE Transactions on Electron Devices (2019),
    % Levinshtein et al, Handbook series on semiconductor parameters, Vol. 2 (1999)
    d.tauP=1/(3e5);
    d.tauN=1/(3e5);
    d.Brad=1e-10*1e-6;
    d.Cp=3e-30*1e-12;
    d.Cn=3e-30*1e-12;

    % saving material used
    d.mat="InGaP";

end