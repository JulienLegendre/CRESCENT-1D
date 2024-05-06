    clear

    %%%%%%%%%%%%%%%%%%%%%%%%%% System definition %%%%%%%%%%%%%%%%%%%%%%%%%%

    % Functions structure_LED and structure_PV provide the structure of the
    % two components.
    % Such functions are used to illustrate how components should be
    % defined. Generally, it is be more efficient to save the structure as
    % a .mat file, and to load it at each simulation.
    % The LED structure is flipped so that the two front layers face each
    % other.
    %run("C:/Users/jlegendre/Documents/Results/2022_12-Hétérojonctions/2023_08 InGaP - InGaAs TPX/Device/InGaP_InGaAs_Device.m");
    d1 = flip(structure_LED());
    d2 = structure_PV();

    % Vacuum gap distance between the LED and the PV cell (m)
    d = 10e-9;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Grids %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Radiation - above-bandgap angular frequency
    % The min. and max. angular frequencies are not set explicitely, but 
    % rather relatively to bandgaps: wmin = weq(1)*min(Eg)/hb and 
    % wmax = weq(2)*max(Eg)/hb. This ensures that the above-bandgap 
    % radiation is correctly defined without needing to have knowledge of 
    % the bandgap energy of each layer.
    % weq(1) can be set below one to include below-bandgap interband
    % transition caused by Urbach tails (only if Eu>0)
    Nwa = 201;                % above-bandgap frequency grid size
    weq = [1 1.1];                                      % equivalent angular frequency bounds (-)

    % Radiation - below-bandgap angular frequency
    % It has no impact on the power output, only on the efficiency.
    % Therefore, its calculation can be switched on and off with
    % computeBelowBandgap. If included, it always go from 1e12 rad.s^-1 to
    % the minimum above-bangap frequency previously defined.
    computeBelowBandgap = false;
    Nwb = 601;

    % Radiation - parallel wavenumber
    Nkr=20*15;              % grid size (should be a multiple of 20)

    % Voltage grid - explicit
    % Having a low step helps with convergence
    U{1}=(0:0.01:1.2).';    % in the LED
    U{2}=(0:0.01:1.2).';    % in the PV cell
    % It is also possible to only provide the grid size and bounds (ignored
    % if U is provided)
    % Nu = [121 121];       % voltage grid size for each component   
    % Umax = [1.2 1.2];     % maximum voltage applied to each component (V)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    opts={"Nw",Nwa,"Nwb",Nwb,"Nkr",Nkr,"U",U,"BelowBandgap",false,...
        "stopAtVoc",true,"iterEf",false,...
        "IterChemPotential",true,"TIatInterfaces",true,"TunnelDepth",50e-9,"RadRecElectrolum",true,"VariableWeights",false,"stopAtVoc",true,"iterEf",false,"kuSetPoint",[1 1;91 1;96 1;101 1;106 1;111 1;116 1;121 1]};

    [Pmax,d1,d2,P,q,eff]=crescent1D(d1,d2,d,weq,opts{:});

    %movefile('TPX_1D_Het3.mat','TPX_testaperture-top.mat');