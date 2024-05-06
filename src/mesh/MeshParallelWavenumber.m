function krmat = MeshParallelWavenumber(d1,d2,d,w,Nkr,options)
    %MESHPARWAVENUMBER Mesh the parallel wavenumber matrix.
    arguments
        d1      (:,1) component                 % emitter structure
        d2      (:,1) component                 % absorber structure
        d       (1,1) double {mustBePositive}   % gap distance (m)
        w       (:,1) double                    % angular frequency (rad.s^-1)
        Nkr     (1,1) uint32                    % parallel wavenumber mesh size
        options (1,1) struct                    % options of crescent1D
    end
    global e c %#ok<GVMIS>

    % Parallel wavevector matrix (rad/m)
    krmat=zeros(length(w),Nkr);
    for i=1:length(w)         
        if d<1e-5 % near field -> propagating and evanescent modes
            nmax=max(real(sqrt(d1.eps(i,1+uint32(d1.hasMirror~="no"):end)))); % max refractive index in the emitter, excluding the mirror if present
            if (w(i)<1e14 || any(d2.Eg/e<0.3)) % close to polaritonic resonances -> large impact of surface modes
                kr=[linspace(options.krrange(1),((1-1e-10)*w(i)/c),Nkr/20)';linspace((1+1e-10)*(w(i)/c),(1-1e-10)*2*(w(i)/c),4/20*Nkr)';logspace(log10(2*(w(i)/c)),log10(max(options.krrange(end),5*(w(i)/c))),1.5*Nkr/2)'];
            elseif i<min([d1.iwm d2.iwm]) % between polaritonic resonances and interband transitions -> eps2 low, mostly propagating and frustrated modes
                kr=[linspace(options.krrange(1),((1-1e-10)*w(i)/c),Nkr/2)';linspace((1+1e-10)*(w(i)/c),(1-1e-10)*1.1*nmax*(w(i)/c),Nkr/2)'];
            else % interband transition -> mostly frustrated modes
                kr=[linspace(options.krrange(1),((1-1e-10)*w(i)/c),3/20*Nkr)';linspace((1+1e-10)*(w(i)/c),(1-1e-10)*nmax*(w(i)/c),Nkr*3/4)';linspace((nmax*(w(i)/c)),(2*nmax*(w(i)/c)),2/20*Nkr)'];
            end
        else % far field -> only propagating modes
            kr=linspace(options.krrange(1),(0.9999*w(i)/c),Nkr)';
        end
        krmat(i,:)=kr;
    end
end