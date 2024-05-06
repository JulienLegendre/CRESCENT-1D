function [w,Nw]=MeshFrequency(d1,d2,wrange,Nw,options)
    %MESHFREQUENCY Mesh the angular frequency, depending on the bandgaps of
    %materials and the mesh size Nw.
    arguments
        d1      (1,:) struct                    % emitter structure
        d2      (1,:) struct                    % absorber structure
        wrange  (1,2) double {mustBePositive}   % range of frequency considered for interband transitions, relative to Eg
        Nw      (1,1) uint32                    % size of frequency mesh related to interband transitions
        options (1,1) struct                    % options of crescent1D
    end
    global hb e %#ok<GVMIS>

    % Extract the bandgap of each material, to ensure that the related
    % frequency are part of the grid.
    if strcmp(options.Solver,"TPX")
        Egmat=unique(sort([Eg_fromFile(d1) Eg_fromFile(d2)]));
    else
        Egmat=unique(sort(Eg_fromFile(d2)));
    end
    Egmat(1:end-1)=Egmat(2:end).*((Egmat(2:end)-Egmat(1:end-1))/e<1e-4)+Egmat(1:end-1).*((Egmat(2:end)-Egmat(1:end-1))/e>=1e-4);
    Egmat=unique(Egmat);
    
    % Interband transition
    % frequency is meshed from wrange(1)*min(Eg/hb) up to wrange(2)*max(Eg/hb)
    % generally, wrange(1) can be set to 1. If Urbach tails are present,
    % interband transitions can occur below Eg and wrange(1) shall be <1.
    Nwm=zeros(length(Egmat),1); % angular frequency mesh size in each range
    Nwm(1)=round(Egmat(1)*(1-wrange(1))/(wrange(end)*max(Egmat)-wrange(1)*min(Egmat))*(Nw-1));
    w=Egmat(1)*linspace(wrange(1),1,Nwm(1)+1)'/hb;
    w=w(1:end-1);
    for i=1:length(Egmat)-1
        Nwm(i+1)=round((Egmat(i+1)-Egmat(i))/(wrange(end)*max(Egmat)-wrange(1)*min(Egmat))*(Nw-1));
        w=[w;linspace(Egmat(i),Egmat(i+1),Nwm(i+1)+1)'/hb]; %#ok<AGROW> 
        w=w(1:end-1);
    end
    w=[w;linspace(Egmat(end),wrange(end)*Egmat(end),Nw-sum(Nwm))'/hb];

    % Below-bandgap radiation is included if option is set to true
    % Minimum frequency is then 10^12 rad.s^-1
    if options.BelowBandgap
        Nwb=[1 -1]*double(round(0.8*(options.Nwb-1)))+double([1 options.Nwb]);
        w=[logspace(12,14,Nwb(1))';linspace(1e14,w(1),Nwb(2))';w];
        w=[w(1:Nwb(1)-1);w(Nwb(1)+1:sum(Nwb));w(sum(Nwb)+2:end)];
        Nw=Nw+uint32(sum(Nwb))-2;
    end
end