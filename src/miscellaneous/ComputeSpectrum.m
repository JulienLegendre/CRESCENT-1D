% Note that this function is not called in any other file, but it is kept
% as it can be useful in certain cases.
% It is aimed to be used inside crescent1D at the end of the calculation, 
% once both radiative and electrical simulation are over. 

function qw = ComputeSpectrum(d1,d2,ku,w,F2,options)
    %COMPUTESPECTRUM Compute the spectral heat flux density exchanged
    %between each layer of the emitter with each layer of the absorber
    arguments
        d1      (1,1) component % first component
        d2      (1,1) component % second component
        ku      (1,2) uint32    % position in voltage grid of point of interest
        w       (:,1) double    % angular frequency (rad.s^-1)
        F2      (:,:,:) double  % transmission function between components (m^-2)
        options (1,1) struct    % options of crescent1D
    end
    global hb %#ok<GVMIS> 

    Nw=length(w);
    wadd=sort(unique([d1.iwm d2.iwm]));
    idxeq=sort([1:Nw wadd]);
    weq=w(idxeq);
    d1.n0=zeros(length(weq),d1.Nl);
    d2.n0=zeros(length(weq),d2.Nl);
    qw=zeros(length(weq),d1.Nl,d2.Nl);
    for k1=1:d1.Nl
        if strcmp(options.Solver,"TPX")
            d1.n0inter=d1.fracInterband(idxeq,k1).*[GBE(weq(1:d1.iwm(k1)+find(wadd==d1.iwm(k1))-1),0,d1.T);GBE(weq(d1.iwm(k1)+find(wadd==d1.iwm(k1)):end),d1.mu(k1,ku(1),ku(2)),d1.T)];
            d1.n0(:,k1)=d1.n0inter+(1-d1.fracInterband(idxeq,k1)).*d1.n0intra(idxeq);
        else
            d1.n0(:,k1)=d1.n0intra(idxeq);
        end
        for k2=1:d2.Nl
            d2.n0inter=d2.fracInterband(idxeq,k2).*[GBE(weq(1:d2.iwm(k2)+find(wadd==d2.iwm(k2))-1),0,d2.T);GBE(weq(d2.iwm(k2)+find(wadd==d2.iwm(k2)):end),d2.mu(k2,ku(1),ku(2)),d2.T)];
            d2.n0(:,k2)=d2.n0inter+(1-d2.fracInterband(idxeq,k2)).*d2.n0intra(idxeq);
            qw(:,k1,k2)=hb*weq.*(d1.n0(:,k1)-d2.n0(:,k2)).*sum(F2(idxeq,d2.zrange{k2},k1),2);           
        end
    end

    %figure
    %area(hb*weq/e,permute(qw,[1 3 2])*e/hb/1e4);
    %legend("Layer "+[1:d2.Nle]+" ("+d2.mat+")");
end