function [Pmax,varargout] = crescent1D(d1,d2,d,wrange,options)
%CRESCENT1D Computes the performance of a 1D thermophotovoltaic (TPV) 
%/ thermophotonic (TPX) device in the near field, for homo or
%heterojunctions. Details on the solver can be found in the scientific 
%paper introducing CRESCENT-1D (ref: to be added) or in Julien Legendre's 
%thesis (https://theses.fr/2023ISAL0094)

    %% Section 1: check arguments
    
    arguments
        % structure
        d1                              (1,:) struct                                                % emitter structure
        d2                              (1,:) struct                                                % absorber structure
        d                               (1,1) double {mustBePositive}                               % gap distance (m)
        
        % frequency, parallel wavenumber, voltage meshes
        wrange                          (1,2) double {mustBePositive}                               % range of angular frequency considered for interband transition, relatively to Eg (w(1)=wrange(1)*min(Eg/hb), w(end)=wrange(end)*max(Eg/hb))
        options.Nw                      (1,1) uint32 = 101                                          % frequency mesh size relative to interband transition
        options.BelowBandgap            (1,1) logical = false                                       % if true -> computes also spectral radiative quantities for low frequency (where phonons and free carriers contribute), down to 10^12 rad.s^-1
        options.Nwb                     (1,1) uint32 = 101                                          % frequency mesh size relative to phonon and free carrier contribution to radiation
        options.Nkr                     (1,1) uint32 = 101                                          % parallel wavevector mesh size
        options.krrange                 (1,2) double {mustBePositive} = [1e3 5/d]                   % extremum parallel wavevector (rad/m)   
        options.U                       (:,1) cell = {}                                             % voltage grid (V), each column being related to a component (if only one column provided -> same mesh for both components)
        options.Umax                    (:,1) double = NaN                                          % max voltage (V) used for linear meshing, each column being related to a component (if only one column provided -> same mesh for both components). Not used in U is provided.
        options.Nu                      (:,1) uint32 = NaN                                          % voltage mesh size, each column being related to a component (if only one column provided -> same mesh for both components). Not used in U is provided.
        
        % solver options
        options.Solver                  (1,1) string {mustBeMember(options.Solver,["TPV","TPX"])} = "TPX"   % set if the emitter d1 is passive or active
        % radiative
        options.PhotonRecycling         (1,1) logical = false                                               % if false, assume the net contribution of photon recycling (i.e. net heat flux between two points of the same component) is zero; if true, compute the net contribution (warning: the implementation of photon recycling has not been thoroughly tested or studied, and must therefore be used carefully)
        % electrical - physics
        options.IterChemPotential       (1,1) logical = true                                                % if true, iterate on the chemical potential of radiation to relax the mu=eU approximation
        options.TIatInterfaces          (1,1) logical = true                                                % if true, include thermionic emission at heterointerfaces
        options.TunnelDepthMax          (1,1) double {mustBeNonnegative} = 5e-8                             % maximum barrier depth through which charge carriers can tunnel; if 0, neglect charge tunnelling
        % electrical - convergence
        options.VariableWeights         (1,1) logical = true                                                % if true, authorise reducing weighting factor value when current values does not allow reaching convergence; if false, stop computation in case of divergence to limit computation time
        options.iterEf                  (1,1) logical = true                                                % if true, authorise to iterate on quasi-Fermi level if non-physical values are obtained; if false, stop computation in case of non-physical value to limit computation time
        % electrical - resistive losses
        options.includeRseries          (1,1) logical = false                                               % if true, include the series resistance to the calculation of electrical power
        options.distanceBusbar          (1,1) double {mustBeNonnegative} = 0                                % distance between busbars for series resistance calculations (m)
        % electrical - other
        options.stopAtVoc               (1,1) logical = false                                               % if true, calculation is only performed for voltages for which the PV cell produces electrical power

        % recombination model
        options.RadRecElectrolum        (1,1) logical = true                                        % if true, radiative recombinations are included with the electroluminescent emission; if false, they are taken into account with a Brad coefficient (Rrad=Brad(np-ni^2))
        options.rmNonradRec             (1,1) logical = false                                       % if true, neglect all non-radiative recombinations
        options.SurfRecPosition         (1,1) string {mustBeMember(options.SurfRecPosition,["rear","front","both","none"])} = "both" % set where surface recombinations shall occur (at remaining interfaces, surface rec. velocity is zero)
        options.nonradTDependencyAuger  (1,1) double = 0                                            % Auger recombination coefficient varies with temperature as T^a, a being this parameter
        options.nonradTDependencySRH    (1,1) double = 0                                            % SRH recombination coefficient varies with temperature as T^a, a being this parameter
        options.nonradTDependencyMatList(:,1) string = ""                                           % list of material showing temperature dependency of non-radiative recombination coefficients
        options.nonradAlloySRH          (1,1) double = NaN                                          % value taken by SRH lifetime for alloys (supersede options.nonradTDependencySRH)
        options.nonradAlloyMatList      (:,1) string = ""                                           % list of alloys showing reduction of SRH lifetime due to alloying

        % tolerance and max iterations
        options.tol                     (1,1) double {mustBePositive} = 1e-6                        % tolerance on the relative variation of charge carrier densities between two iterations
        options.tolMu                   (1,1) double {mustBePositive} = 1e-4                        % tolerance on the absolute variation of the chemical potential between two iterations
        options.iterVmax                (1,1) double {mustBeInteger,mustBePositive} = 100           % maximum number of iterations for solving Poisson equation
        options.iterEmax                (1,1) double {mustBeInteger,mustBePositive} = 120           % maximum number of iterations for solving the set of electrical equations (drift-diffusion, continuity, Poisson)
        options.iterMumax               (1,1) double {mustBeInteger,mustBePositive} = 100           % maximum number of iterations to obtain a converged chemical potential of radiation

        % mirror definition
        options.MirrorMaterial          (1,1) string {mustBeMember(options.MirrorMaterial,["Perfect","Au"])} = "Perfect"            % material used for the back mirror (only perfect has been thoroughly tested, Au has been used for TPV only)
        options.MirrorPosition          (1,1) string {mustBeMember(options.MirrorPosition,["left","right","both","none"])} = "both" % position of the back-mirror (only the case with both mirrors being perfect has been tested for TPX, the other situations being used in the case of semi-infinite components)
        options.MirrorThickness         (1,1) double {mustBePositive} = 1e-7                                                        % mirror thickness (m)
        options.MirrorMeshsize          (1,1) uint32 = 51                                                                           % radiative spatial mesh size for mirrors
        
        % saved data
        options.kuSetPoint              (:,2) double {mustBeInteger,mustBeNonnegative} = [0 0];     % sets of voltage indices [ku1 ku2] for which the spatial quantities obtained in the drift-diffusion solver shall be saved
    end

    timerVal=tic;    
    warning('backtrace','off');
    run("FundamentalConstants.m"); % loads the fundamental constants
    
    % translates some of the input data for simplicity in usage
    Nw=options.Nw;
    Nkr=options.Nkr;
    switch options.MirrorPosition
        case "left"
            options.m=[1 0];
        case "right"
            options.m=[0 1];
        case "both"
            options.m=[1 1];
        case "none"
            options.m=[0 0];
    end

    %% Section 2: Prepare calculation                      

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Frequency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [w,Nw]=MeshFrequency(d1,d2,wrange,Nw,options); % sets the frequency mesh               
    
    %%%%%%%%%%%%%%%%%%%%%%%% Emitter properties %%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create the device corresponding to the emitter (passive or LED)
    % Computes important material properties
    clear setIndex@device
    if strcmp(options.Solver,"TPX")
        d1=activeComponent(d1,1,w,options);
    else
        d1=component(d1,1,w,options);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%% PV cell properties %%%%%%%%%%%%%%%%%%%%%%%%%

    % Create the device corresponding to the PV cell
    d2=activeComponent(d2,2,w,options);

    %%%%%%%%%%%%%%%%%%%%%%%% Parallel wavenumber %%%%%%%%%%%%%%%%%%%%%%%%%
    krmat = MeshParallelWavenumber(d1,d2,d,w,Nkr,options);

    %%%%%%%%%%%%%%%% System definition (radiative solver) %%%%%%%%%%%%%%%%

    % Dielectric function matrix
    epsMat=[ones(Nw,1)*ones(1,~d1.isSemiInf) d1.eps ones(Nw,1) d2.eps ones(Nw,1)*ones(1,~d2.isSemiInf)];
    % Thickness matrix
    t=[0 ones(~d1.isSemiInf,1)*d1.tj(1) d1.tj(2:end) d d2.tj(1:end-1) ones(~d2.isSemiInf,1)*d2.tj(end)];
    % Position matrix
    z=cumsum(t);
    
    disp(join(['Creation of the system: done (',string(duration(seconds(toc(timerVal)),'Format','hh:mm:ss')),').']));

    %% Section 3: Compute the radiative transmission coefficient/function

    % We use the fluctuational electrodynamics framework to obtain the 1D
    % radiative heat transfer in the near field. More precisely, we compute
    % here the transmission function F between two layers. It must be
    % multiplied by the difference of Bose-Einstein distributions to obtain
    % the net photon flux density. Because the local radiative generation
    % and recombination rate are required for electrical simulations, the
    % photon flux density must be spatially resolved.

    %%%%%%%%%%%%%%%% Transmission emitter - absorber %%%%%%%%%%%%%%%%%%%%

    % Determination of the transmission coefficient
    d2Phi=zeros(max(d2.Nzr),Nkr);                     	% Equivalent transmission coeff per unit kr,w (/(m^2.s.rad/m.rad/s))
    T12=zeros(Nw,Nkr);
    F=zeros(Nw,length(d1.zr)-1,length(d2.zr)-1);

    %Fback=zeros(Nw,length(d2.zr)-1);
    Fin1=zeros(Nw,length(d1.zr)-1,length(d1.zr)-1);
    Fin2=zeros(Nw,length(d2.zr)-1,length(d2.zr)-1);

    % Loop on the layers (outer loop: emitter, inner loop: absorber)
    for k1=1+(d1.hasMirror=="perfect"):d1.Nl
        for k2=1:d2.Nl-uint32(d2.hasMirror=="perfect")

            % Loop on the pulsation
            for i=1:Nw
                Fkr=zeros(d1.Nzr(k1)-1,d2.Nzr(k2)-1,Nkr);
                % Loop on parallel wavenumber
                for j=1:Nkr
                    [d2Phi(1:d2.Nzr(k2),j),varTemp]=ComputeTransmission(uint32(~d1.isSemiInf)+k1-1,uint32(~d1.isSemiInf)+d1.Nl+k2,w(i),krmat(i,j),epsMat(i,:),z,d1.zrl{k1},d2.zrl{k2});
                    % various cases are considered, depending if one of
                    % the layers is semi-infinite
                    if (d1.isSemiInf && k1==1) && (d2.isSemiInf && k2==d2.Nl)
                        Fkr(1,1,j)=d2Phi(1,j);
                    elseif (d1.isSemiInf && k1==1)
                        Fkr(1,:,j)=d2Phi(1:d2.Nzr(k2)-1,j)-d2Phi(2:d2.Nzr(k2),j);
                    elseif (d2.isSemiInf && k2==d2.Nl)
                        Fkr(:,1,j)=(d1.zrl{k1}(2:end)-d1.zrl{k1}(1:end-1))/2.*(varTemp(2:end,1)+varTemp(1:end-1,1));
                    else
                        Fkr(:,:,j)=repmat(d1.zrl{k1}(2:end)-d1.zrl{k1}(1:end-1),[1,d2.Nzr(k2)-1])/2.*(varTemp(2:end,1:end-1)+varTemp(1:end-1,1:end-1)-varTemp(2:end,2:end)-varTemp(1:end-1,2:end));
                    end
                end
                
                % Transmission coefficient between the two devices
                if k2==1
                   T12(i,:)=T12(i,:)+4*pi^2*d2Phi(1,:)./krmat(i,:); 
                end

                % Transmission function of photons between each elementary
                % layers of each device, obtained for each frequency
                F(i,d1.zrange{k1},d2.zrange{k2})=trapz(krmat(i,:),Fkr,3);

            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%% Inner transmission %%%%%%%%%%%%%%%%%%%%%%%%
    % Computes the transmission between two thin layer of a common component
    % This has not been thoroughly tested or studied, and may therefore contain errors.

    % Since components are supposed to be maintained at a constant
    % temperature, only differences in chemical potential of radiation can
    % cause a net radiative heat flux. Therefore, calculations are only
    % performed for frequency for which interband transitions occur.

    % Because the heat transfer diverges when the distance between the 
    % emitter and the absorbed approaches 0 due to evanescent modes,
    % integration over kr is only performed for modes that are propagative
    % in the component (i.e. between 0 and nw/c).

    if options.PhotonRecycling
        warning("The implementation of photon recycling has not been thoroughly tested or studied, and must therefore be used with caution.")
        if options.Solver=="TPX"
            % in the LED
            for i=min(d1.iwm):Nw
                % only the modes propagating inside at least one layer is
                % considered, as including truly evanescent modes causes
                % divergence (due to the locality approximation)
                kr=linspace(options.krrange(1),(0.9999*max(real(sqrt(d1.eps(i,:))))*w(i)/c),Nkr)';
                % outer loop: emitter, inner loop: absorber
                for k1=1:d1.Nl
                    for k11=k1+1:d1.Nl
                        Fkr1=zeros(d1.Nzr(k1)-1,d1.Nzr(k11)-1,Nkr);
                        for j=1:Nkr
                            [d2Phi,varTemp]=ComputeTransmission(~d1.isSemiInf*(options.m(1)+1)+k1-1,~d1.isSemiInf*(options.m(1)+1)+k11-1,w(i),kr(j),epsMat(i,:),z,d1.zrl{k1},d1.zrl{k11});
                            if d1.isSemiInf && k1==1
                                Fkr1(:,1,j)=d2Phi(1:d1.Nzr(k11)-1,j)-d2Phi(2:d1.Nzr(k11),j);
                            else
                                Fkr1(:,:,j)=repmat(d1.zrl{k1}(2:end)-d1.zrl{k1}(1:end-1),[1,d1.Nzr(k11)-1])/2.*(varTemp(2:end,1:end-1)+varTemp(1:end-1,1:end-1)-varTemp(2:end,2:end)-varTemp(1:end-1,2:end));
                            end
                        end
                        Fin1(i,d1.zrange{k1},d1.zrange{k11})=trapz(kr,Fkr1,3);
                    end
                end
            end
            Fin1=Fin1+permute(Fin1,[1 3 2]);
        end
        % in the PV cell
        for i=min(d2.iwm):Nw
            kr=linspace(options.krrange(1),(0.9999*max(real(sqrt(d2.eps(i,:))))*w(i)/c),Nkr)';
            for k2=1:d2.Nl
                for k22=k2+1:d2.Nl
                    Fkr2=zeros(d2.Nzr(k2)-1,d2.Nzr(k22)-1,Nkr);
                    for j=1:Nkr
                        [~,varTemp]=ComputeTransmission(~d1.isSemiInf*(options.m(1)+1)+d1.Nl+k2,~d1.isSemiInf*(options.m(1)+1)+d1.Nl+k22,w(i),kr(j),epsMat(i,:),z,d2.zrl{k2},d2.zrl{k22});
                        if d2.isSemiInf && k22==d2.Nl
                            Fkr2(:,1,j)=repmat(d2.zrl{k2}(2:end)-d2.zrl{k2}(1:end-1),[1,d2.Nzr(k22)-1])/2.*(varTemp(2:end,1)+varTemp(1:end-1,1));
                        else
                            Fkr2(:,:,j)=repmat(d2.zrl{k2}(2:end)-d2.zrl{k2}(1:end-1),[1,d2.Nzr(k22)-1])/2.*(varTemp(2:end,1:end-1)+varTemp(1:end-1,1:end-1)-varTemp(2:end,2:end)-varTemp(1:end-1,2:end));
                        end
                    end
                    Fin2(i,d2.zrange{k2},d2.zrange{k22})=trapz(kr,Fkr2,3);
                end
            end
        end
        Fin2=Fin2+permute(Fin2,[1 3 2]);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%% Simplifications %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Since we only account for one chemical potential per layer, we can do
    % some sommations independently of voltage (to decrease computational time)

    % F1 and F2 represent the simplified transmission function. For F1 for
    % instance, it is the transmission function at frequency w(i) (first coordinate) 
    % of photons between the mth elementary layer of component d1 (second
    % coordinate) and the nth complete layer of component d2 (third
    % coordinate).
    F1=zeros(Nw,length(d1.zr)-1,d2.Nl);
    F2=zeros(Nw,length(d2.zr)-1,d1.Nl);
    for k2=1:d2.Nl
        F1(:,:,k2)=sum(F(:,:,d2.zrange{k2}),3);
    end
    for k1=1:d1.Nl
        F2(:,:,k1)=permute(sum(F(:,d1.zrange{k1},:),2),[1 3 2]);
    end

    % As for F1 and F2, but for exchanges within a component.
    Fi1=zeros(Nw,length(d1.zr)-1,d1.Nl);
    for k1=1:d1.Nl
        for k11=1:d1.Nl
            Fi1(:,d1.zrange{k1},k11)=sum(Fin1(:,d1.zrange{k1},d1.zrange{k11}),3);
        end
    end
    Fi2=zeros(Nw,length(d2.zr)-1,d2.Nl);
    for k2=1:d2.Nl
        for k22=1:d2.Nl
            Fi2(:,d2.zrange{k2},k22)=sum(Fin2(:,d2.zrange{k2},d2.zrange{k22}),3);
        end
    end

    clear Fkr Fkr2 F Fin1 Fin2 PhiBack
    clear kr d2Phi varTemp
    disp(join(['Computation of the near-field transmission coefficient: done (',string(duration(seconds(toc(timerVal)),'Format','hh:mm:ss')),').']));

    %% Section 4: Compute the equilibrium electrostatic potential profile

    if strcmp(options.Solver,"TPX")
        d1=DriftDiffHetEq(d1,options.tol,options.iterVmax);    % device 1 (LED)
    end
    d2=DriftDiffHetEq(d2,options.tol,options.iterVmax);    % device 2 (PV cell)

    disp('Computation of the equilibrium profiles: done.');

    %% Section 5: Compute the electrical characteristic of the active components
    % The electrical characteristic of the active components is obtained
    % through the resolution of the drift-diffusion equations in 1D, along
    % with Poisson equation and continuity equations;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% Initialisation %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Tag changed to 'over' if convergence issues are found, which will
    % stop the calculations.
    calc_state='on';

    % Initialisation of the various voltage-dependent quantities.
    d1=initElecQuantities(d1,d1.Nu,d2.Nu);
    d2=initElecQuantities(d2,d1.Nu,d2.Nu);

    % Temporary vector for saving the chemical potential of photons in the
    % various layers of the PV cell.
    mu20=zeros(d2.Nl,1);

    % Heat flux density between the two components
    q12=zeros(d1.Nu,d2.Nu);
    % Number of iterations on the chemical potential to reach convergence
    iterMuMat=zeros(d1.Nu,d2.Nu);
    % At each iteration, Pmax can be updated if the new power output is
    % greater than the previous one. This allows ensuring that the
    % calcPoint at max. power point is saved.
    % Warning: because series resistance are not included directly within
    % the iterative process, the calcPoint at max. power is obtained
    % without series resistance, even though options are set otherwise. The
    % maximum power however includes series resistance.
    Pmax=0;

    % Identifier for the calcPoint which must be saved, as set by the
    % option kudata
    idxdata=0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% General %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for ku1=1:d1.Nu

        % The first guess to compute the solution of the drift-diffusion
        % equations at (d1.U(ku1),0) is the closest solution obtained,
        % either (d1.U(ku1-1),0) or the equilibrium conditions.
        d1=setGuessOuterLoop(d1,ku1,mu20);
        d2=setGuessOuterLoop(d2,ku1,mu20);

        % Reset the weights used to smoothen the iterative process
        % to original value
        wmu=0.3;
        d1=resetWeights(d1,1,1e-3,1e-5,2e-1);
        d2=resetWeights(d2,1,1e-3,1e-5,2e-1);

        for ku2=1:d2.Nu
            
            % The first guess for the chemical potential is eU for the LED
            % (as mainly directed by the voltage) while it is the chemical
            % potential obtained for the closest solution for the PV cell
            % (either at (ku1,ku2-1), (ku1-1,0) if ku2=1, or 0 if
            % ku1=ku2=1).
            d1=setGuessInnerLoop(d1,[ku1 ku2],options);
            d2=setGuessInnerLoop(d2,[ku1 ku2],options);
            
            % Put current case in memory, to be used in case of convergence
            % issues
            d1s=d1;
            d2s=d2;

            % Initialisation of the iterative process on the chemical
            % potential
            errMu=1;
            iterMu=0;

            % At each iteration, the chemical potential is obtained through
            % the resolution of the drift-diffusion equations. The
            % iterative process stops once the maximum variation of the
            % chemical potential becomes lower than a given tolerance,
            % which is reduced for non-zero weighting factors.
            while errMu>=options.tolMu*wmu && iterMu<options.iterMumax
                
                % Correction of Bose-Einstein distributions
                d1=d1.Computen0(w,d1.mu(:,ku1,ku2));
                d2=d2.Computen0(w,d2.mu(:,ku1,ku2));

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if strcmp(options.Solver,"TPX")
                    [d1,DDstate] = d1.SolveDriftDiff(d2,F1,Fi1,w,[ku1 ku2],options);
                    % Check if we should stop
                    if DDstate=="over - divergence"
                        calc_state='over';
                        break
                    end
                end
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PV cell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [d2,DDstate,q12(ku1,ku2)] = d2.SolveDriftDiff(d1,F2,Fi2,w,[ku1 ku2],options);
                % Check if we should stop
                if DDstate=="over - divergence"
                    calc_state='over';
                    break
                end

                %%%%%%%%%%%%%%%%%%%%%%%%% Chemical potential %%%%%%%%%%%%%%%%%%%%%%%%
                % Updates the chemical potential with the quasi-Fermi
                % levels
                [d1,errMu1]=d1.updateChemicalPotential([ku1 ku2],wmu,iterMu,options);
                [d2,errMu2]=d2.updateChemicalPotential([ku1 ku2],wmu,iterMu,options);
                errMu=max(errMu1,errMu2);
                iterMu=iterMu+1;
                % If convergence is difficult, reducing the weight and
                % do again the iterative process from the start
                if iterMu>options.iterMumax/2 && wmu==1
                    wmu=0.125;
                    d1=d1s;
                    d2=d2s;
                    iterMu=0;
                    errMu=1;
                    warning(join(['Difficult convergence (chem. potential loop, ku1=',num2str(ku1),', ku2=',num2str(ku2),'). Decreasing factor w_mu.']));
                end
                % Stop the process if required to perform the calculation
                % in one pass
                if ~options.IterChemPotential
                    errMu=options.tolMu/10;
                end
            end
            % Save the number of iterations required to reach convergence
            iterMuMat(ku1,ku2)=iterMu;
            if iterMu==options.iterMumax
                warning(join(['Not able to converge. Moving on to next voltage (chem. potential loop, ku1=',num2str(ku1),', ku2=',num2str(ku2),').']));
                calc_state='over';
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sum up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % If issues of convergence, quantities not calculated are set
            % to NaN.
            if strcmp(calc_state,'over')
                d1=d1.setRatesDivergence([ku1 ku2]);
                d2=d2.setRatesDivergence([ku1 ku2]);
                if ku2>1
                    calc_state='on';
                end
                break
            end
            
            % Set the value of total Generation-Recombination rates
            d1=d1.setRatesConvergence([ku1 ku2]);
            d2=d2.setRatesConvergence([ku1 ku2]);
            
            % save spatial information (in calcPoint) at max power point
            if options.Solver=="TPX" && d1.U(ku1)*d1.J(ku1,ku2)+d2.U(ku2)*d2.J(ku1,ku2)>Pmax
                Pmax=d1.U(ku1)*d1.J(ku1,ku2)+d2.U(ku2)*d2.J(ku1,ku2);
                d1.MaxPoint=d1.cp;
                d2.MaxPoint=d2.cp;
            elseif options.Solver=="TPV" && d2.U(ku2)*d2.J(ku1,ku2)>Pmax
                Pmax=d2.U(ku2)*d2.J(ku1,ku2);
                d2.MaxPoint=d2.cp;
            end
            % save spatial information (in calcPoint) for the points set by
            % the option kudata
            isOfInterest=([ku1 ku2]==options.kuSetPoint);
            if any(logical(prod(isOfInterest,2)))
                idxdata=idxdata+1;
                if options.Solver=="TPX"
                    d1.SetPoint(idxdata)=d1.cp;
                end
                d2.SetPoint(idxdata)=d2.cp;
            end

            % Resistive losses, modelled accordingly to Milovich et al. (2020)
            % We consider that the electrical contact covers the complete
            % back surface, so that current crowding caused by the metallic
            % grid only occurs in the front surface.
            % This model is very simple, and further work on it shall be
            % done to obtain a realistic estimation of the electrical
            % resistance.
            if strcmp(options.Solver,"TPX")
                d1.Rseries(ku1,ku2)=trapz(d1.zea,1./(e^2/kb/d1.T*(d1.cp.n.*d1.cp.D(:,1)+d1.cp.p.*d1.cp.D(:,2))))+...
                                                 1/12/(e^2/kb/d1.T*(d1.cp.n(end).*d1.cp.D(end,1)+d1.cp.p(end).*d1.cp.D(1,end)))*(options.distanceBusbar)^2/d1.tj(end);
            end
            d2.Rseries(ku1,ku2)=trapz(d2.zea,1./(e^2/kb/d2.T*(d2.cp.n.*d2.cp.D(:,1)+d2.cp.p.*d2.cp.D(:,2))))+...
                                             1/12/(e^2/kb/d2.T*(d2.cp.n(1).*d2.cp.D(1,1)+d2.cp.p(1).*d2.cp.D(1,2)))*(options.distanceBusbar)^2/d2.tj(1);

            if ku2==1
                % save this result for first guess in iteration (ku1+1,1)
                d1=d1.saveGuessOuterLoop([ku1 ku2]);
                [d2,mu20]=d2.saveGuessOuterLoop([ku1 ku2]);
            end

            % if option stopAtVoc is set to true, the iteration on the PV
            % cell voltage will stop as soon as the PV cell voltage exceeds
            % Voc. The remaining points at larger PV cell voltage are not
            % computed, and the related data therefore set to NaN.
            if options.stopAtVoc && d2.J(ku1,ku2)<0
                d1=d1.setRatesDivergence([ku1 ku2+1]);
                d2=d2.setRatesDivergence([ku1 ku2+1]);
                calc_state='on';
                break;
            end

        end
        if strcmp(calc_state,'over')
            break 
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IV and PV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Computes the electrical power input/output of each component, along
    % with the short-circuit current and open-circuit voltage.
    d1=d1.ComputeIV(d2);
    d2=d2.ComputeIV(d1);
    
    % Computes the overall power output of the device
    if strcmp(options.Solver,"TPX")
        P=(d1.P+d2.P);
    else
        P=d2.P;
    end
    % Pmax is forced to be negative so that this function can be used with
    % minimisation processes for optimisation.
    Pmax=-Pmax;
    % Modification of the power output if series resistance are included
    if options.includeRseries
        if strcmp(options.Solver,"TPX")
            d1.Us=d1.U+d1.Rseries.*d1.J;
            d2.Us=d2.U'+d2.Rseries.*d2.J;
            P=-(d1.Us.*d1.J+d2.Us.*d2.J);
        else
            d2.Us=d2.U'+d2.Rseries.*d2.J;
            P=-d2.Us.*d2.J;
        end
        Pmax=-max(P(:));
    end
    % ku1max and ku2max represent the indices for which power is maximum.
    [ku1max,ku2max]=ind2sub([d1.Nu d2.Nu],find(P==max(P(:)))); %#ok<ASGLU> 

    % Computes the overall efficiency of the device
    if strcmp(options.Solver,"TPX")
        eff=P./(q12+d1.P);
    else
        eff=P./q12;
    end
    
    %% Section 6: End
    
    % Set the function facultative outputs
    varargout{1}=d1;
    varargout{2}=d2;
    varargout{3}=P;
    varargout{4}=q12;
    varargout{5}=eff;
    
    % Show computational time
    timerVal=toc(timerVal);
    disp(join(['End of calculation, performed in',string(duration(seconds(timerVal),'Format','hh:mm:ss')),'.']));
    
    % Clear useless data and save useful data
    clear d1s d2s DDflag errMu1 errMu2 calc_state i idxdata idxeq iterMu j k1 k11 k2 k22 ku1 ku2 mu20 timerVal
    save data_crescent1D.mat
    
end
