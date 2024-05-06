function [eps,Eg,varargout] = epsInGaP(x,w,T,Ndop,DopSign,mx,muMaj,Eu)
%EPSINGAP Computes the dielectric function of In_{x}Ga_{1-x}P
% varargout{1}: first index with interband transition (in terms of w)
% varargout{2}: fraction of absorption coming from interband transition
    arguments
        x (1,:) double {mustBeInRange(x,0,1)}               % Ga fraction (-) 
        w (:,1) double {mustBePositive}                     % frequency (rad.s^-1)
        T (1,:) double {mustBePositive}                     % temperature (K)
        Ndop (1,:) double {mustBePositive}                  % doping level per layer (m^-3)
        DopSign (1,:) double {mustBeMember(DopSign,[-1;1])} % doping type (-1 -> n, 1 -> p)
        mx (1,2) double {mustBePositive}                    % electron hole effective mass [m_electrons m_holes] (kg)
        muMaj (1,:) double {mustBePositive}                 % majority carrier mobility (m^2.V^-1.s^-1)
        Eu (1,1) double = nan                               % Urbach energy (eV)
    end

    global hb e c eps0

    % Lorentz-Drude
    weq=w/(2*pi*100*c)*ones(1,length(Ndop));
    % Lorentz
    % warning: obtained for x=0.51
    % data from Ferrini et al, The European Physical Journal B - Condensed
    % Matter (2002).
    epsInf=9.61-0.5*x;                                                             % dielectric constant (high frequency)
    sG=112.5;
    wTOG=370;
    sI=525;
    wTOI=333;
    gG=10.4;
    gI=9.2;
    % Drude
    mx=mx(1)*(1-DopSign)/2+mx(2)*(DopSign+1)/2;
    wp=ones(length(weq),1)*sqrt(Ndop*e^2./(mx*epsInf*eps0))/(2*pi*100*c);
    gammap=ones(length(weq),1)*(e./(mx.*muMaj))/(2*pi*100*c);
    eps1=epsInf*(sG^2./(wTOG^2-weq.^2-1i*weq*gG)+sI^2./(wTOI^2-weq.^2-1i*weq*gI)-wp.^2./(weq.*(weq+1i*gammap)));

    % Interband
    [eps2,Eg]=epsInGaP_Schubert(x,w,T);

    % Correction with Urbach tail
    if Eu>0
        wc=sort(unique([w;linspace(min(w(1),(Eg-0.3*e)/hb),w(end)+0.2*e/hb,10001)']));
        eps2=epsInGaP_Schubert(x,wc,T);
        eps2(hb*wc<=Eg)=real(eps2(hb*wc<=Eg));
        a=2*wc/c.*imag(sqrt(eps2));
        E=hb*wc/e;
        G=zeros(size(E));
        for i=1:length(E)
            G(i)=1/4/Eu*sum((E(2:end)-E(1:end-1)).*(a(1:end-1).*exp(-abs((E(i)-E(1:end-1))/Eu))+a(2:end).*exp(-abs((E(i)-E(2:end))/Eu))));
        end
        eps2=real(eps2)+1i*G*c./wc.*sqrt((G*c/2./wc).^2+real(eps2));
        eps2=interp1(wc,eps2,w);
        varargout{1}=max(find(imag(eps2)>=1e-8,1)-1,1);
    elseif Eu==0
        eps2(hb*w<Eg)=real(eps2(hb*w<Eg));
        varargout{1}=find(w>=Eg/hb,1);
    else
        eps2(imag(eps2)<0)=real(eps2(imag(eps2)<0));
        varargout{1}=find(w>=Eg/hb,1);
    end

    if isempty(varargout{1})
        varargout{1}=length(w);
    end
    
    % Total
    % varargout{2} corresponds to the fraction of eps2 which corresponds to
    % interband transition. It is used to segregate between
    % electroluminescent and thermal above-bandgap radiation.
    varargout{2}=imag(eps2*ones(1,length(Ndop)))./imag(eps1+eps2*ones(1,length(Ndop)));
    eps=eps1+eps2*ones(1,length(Ndop));

end

