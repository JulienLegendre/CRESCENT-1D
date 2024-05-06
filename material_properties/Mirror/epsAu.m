function eps = epsAu(w)
%EPSAU Returns the dielectric function of gold
    arguments
        w (:,1) double % angular frequency (rad.s^-1)
    end

    global hb e %#ok<GVMIS>

    % Data from Derkachova et al, Plasmonics (2016)
    eps0=9.84;
    wp=9.01; %eV
    gamma=0.072; %eV

    eps=eps0-wp^2./((hb*w/e).^2+1i*gamma*hb*w/e)+5.6i./(1+exp(-(hb*w/e-2.4)/0.17));

end