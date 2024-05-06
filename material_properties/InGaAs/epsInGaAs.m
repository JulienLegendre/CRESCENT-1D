function [eps,Eg,varargout]=epsInGaAs(x,w,T,N,Nsign,mx,muMaj,Eu)
%EPSINGAAS Computes the dielectric function of In_{x}Ga_{1-x}As
    arguments
        x (1,1) double {mustBeInRange(x,0,1)}             % fraction of Al (-)
        w (:,1) double {mustBePositive}                     % frequency (rad.s^-1)
        T (1,1) double {mustBePositive}                     % temperature (K)
        N (1,:) double {mustBePositive}                     % doping level per layer (m^-3)
        Nsign (1,:) double {mustBeMember(Nsign,[-1;1])}     % doping type (-1 -> n, 1 -> p) 
        mx (1,2) double {mustBePositive}                    % electron hole effective mass [m_electrons m_holes] (kg)
        muMaj (1,:) double {mustBePositive}                 % majority carrier mobility (m^2.V^-1.s^-1)
        Eu (1,1) double = nan                               % Urbach energy (eV)
    end

    global c hb eps0 e %#ok<GVMIS>
    
    %% Drude-Lorentz
    
    % Lorentz-Drude
    weq=w/(2*pi*100*c)*ones(1,length(N));
    % Lorentz
    % data from Adachi, GaAs and related materials (1994), 
    % Adachi, Optical constants of crystalline and amorphous semiconductors (1999),
    % Levinshtein et al, Handbook series on semiconductor parameters, Vol. 2 (1999)
    epsInf=10.9+1.4*x; % epsInf for GaAs comes from Adachi1994, while temperature dependence comes from Levinshtein1999
    wLOG=292.2;
    wTOG=269.2;
    wLOI=242.5;
    wTOI=218.5;
    gG=3.3;
    gI=4.9;
    % Drude
    mx=mx(1)*(1-Nsign)/2+mx(2)*(Nsign+1)/2;
    wp=ones(length(weq),1)*sqrt(N*e^2./(mx*epsInf*eps0));
    gammap=ones(length(weq),1)*(e./(mx.*muMaj));
    epsIntra=epsInf*((1-x)*(wLOG^2-wTOG^2)./(wTOG^2-weq.^2-1i*weq*gG)+x*(wLOI^2-wTOI^2)./(wTOI^2-weq.^2-1i*weq*gI)-wp.^2./(w.*(w+1i*gammap)));
    
    %% Interband
    
    [epsInter,Eg]=epsInGaAsSb_GonzalezCuevas(x,1,w,T);

    % Correction with Urbach tail
    if Eu>0
        wc=sort(unique([w;linspace(min(w(1),(Eg-0.3*e)/hb),w(end)+0.2*e/hb,10001)']));
        epsInter=epsAlGaAsSb_GonzalezCuevas(x,1,wc,T);
        epsInter(hb*wc<=Eg)=real(epsInter(hb*wc<=Eg));
        a=2*wc/c.*imag(sqrt(epsInter));
        E=hb*wc/e;
        G=zeros(size(E));
        for i=1:length(E)
            G(i)=1/4/Eu*sum((E(2:end)-E(1:end-1)).*(a(1:end-1).*exp(-abs((E(i)-E(1:end-1))/Eu))+a(2:end).*exp(-abs((E(i)-E(2:end))/Eu))));
        end
        epsInter=real(epsInter)+1i*G*c./wc.*sqrt((G*c/2./wc).^2+real(epsInter));
        epsInter=interp1(wc,epsInter,w);
        varargout{1}=max(find(imag(epsInter)>=1e-8,1)-1,1);
    elseif Eu==0
        epsInter(hb*w<Eg)=real(epsInter(hb*w<Eg));%+1i*1e-15;
        varargout{1}=find(w>=Eg/hb,1);
    else
        epsInter(imag(epsInter)<0)=real(epsInter(imag(epsInter)<0));
        varargout{1}=find(w>=Eg/hb,1);
    end

    if isempty(varargout{1})
        varargout{1}=length(w);
    end
    
    % Total
    % varargout{2} corresponds to the fraction of eps2 which corresponds to
    % interband transition. It is used to segregate between
    % electroluminescent and thermal above-bandgap radiation.
    varargout{2}=imag(epsInter*ones(1,length(N)))./imag(epsIntra+epsInter*ones(1,length(N)));
    eps=epsIntra+epsInter*ones(1,length(N));
    
end

