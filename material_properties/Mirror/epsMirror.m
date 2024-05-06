function eps = epsMirror(strMaterial,w)
%EPSMIRROR Computes the dielectric function of the mirror, depending on the
%material chosen.
    arguments
        strMaterial (1,1) string {mustBeMember(strMaterial,["Perfect","Au"])}
        w (:,1) double {mustBePositive}
    end

    global c hb e %#ok<GVMIS>
    
    switch strMaterial
        case "Perfect"
            eps=1e6*(7e8*c*hb/e)^2*ones(size(w));
        case "Au"
            eps=epsAu(w);
    end

end