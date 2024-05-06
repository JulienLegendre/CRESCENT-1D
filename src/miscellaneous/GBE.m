function n0 = GBE(w,mu,T)
%GBE Compute the generalised Bose-Einstein distribution
    arguments
        w (:,1) double  % angular frequency (rad.s^-1)
        mu (1,:) double % electrochemical potential (eV)
        T (1,1) double  % temperature (K)
    end
    
    global hb kb e %#ok<GVMIS>

    n0 = 1./(exp((hb*w-e*mu)/(kb*T))-1);
    
end

