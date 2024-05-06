classdef activeComponent < component
    %ACTIVECOMPONENT Active component of a radiative heat engine. Children
    %of class Component, used in the crescent1D function.

    properties
        Ar              (:,2) double    % effective Richardson coefficient (A.m^-2.K^-2)
        Bm              (:,1) double    % Brad on electrical spatial grid (m^3.s^-1)
        Brad            (1,:) double    % radiative recombination coefficient per layer (m^-3.s^-1)
        chi             (1,:) double    % electron affinity (J)
        Cm              (:,2) double    % [Cn Cp] on electrical spatial grid (m^6.s^-1)
        Cn              (1,:) double    % electron-electron-hole Auger recombination coefficient per layer (m^6.s^-1)
        Cp              (1,:) double    % electron-hole-hole Auger recombination coefficient per layer (m^6.s^-1)
        CTmodel         (1,:) struct    % Caughey-Thomas model for low-field mobility calculation, per layer
        CTm             (1,1) struct    % CTmodel on electrical spatial grid
        Dg              (:,1) double    % =xin+xip, on total electrical spatial grid (2 points per interface)
        Dgv             (:,1) double    % Dg, but on standard grid (1 point per interface)
        dz              (:,1) double    % spatial step in electrical spatial grid (m)
        dza             (:,1) double    % dz, but on total electical spatial grid (m)
        dzb             (:,1) double    % dz normalised by the intrinsic Debye length (m)
        Eg              (1,:) double    % bandgap energy (J)
        epsr            (1,:) double    % static dielectric constant (-)
        eta             (:,2) double    % thermionic emission - correction term due to asymmetry of properties
        etat            (:,2) double    % tunnelling - correction term due to asymmetry of properties
        fracInterband   (:,:) double    % fraction of emission/absorption at a given frequency being caused by interband transitions
        xin            (:,1) double    % heterostructure - correction term in electron concentration due to asymmetry of properties 
        xinv           (:,1) double    % xin, but on the standard electrical spatial grid
        J               (:,:) double    % current density (A.m^-2)
        Jsc             (:,:) double    % short-circuit current (A.m^-2)
        K               (:,3) double    % weights for central finite difference of second derivative
        KL              (:,3) double    % weights for backward finite difference of first derivative
        KR              (:,3) double    % weights for forward finite difference of first derivative
        xip            (:,1) double    % heterostructure - correction term in hole concentration due to asymmetry of properties 
        xipv           (:,1) double    % xip, but on the standard electrical spatial grid
        locSurfRec      (2,1) logical   % set where surface recombination occur (first value for z=0, last for z=t, true if surface recombinations occur)
        MaxPoint        (1,1) calcPoint % spatial data at maximum power point of the global device
        mx              (:,2) double    % effective masses (electron and hole) (kg)
        Nc              (1,:) double    % effective density of states in the conduction band, per layer (m^-3)
        Ncm             (:,1) double    % Nc, but on the total electrical spatial grid (m^-3)
        ni              (1,:) double    % intrinsic carrier concentration, per layer (m^-3)
        nim             (:,1) double    % ni, but on the total electrical spatial grid (m^-3)
        NleRange        (1,:) double    % indices of the layers composing the electrical system
        Nle             (1,1) double    % number of layers composing the electrical system (-)
        Nm              (:,1) double    % effective doping level (m^-3)
        Nv              (1,:) double    % effective density of states in the valence band, per layer (m^-3)
        Nvm             (:,1) double    % Nv, but on the total electrical spatial grid (m^-3)
        Nze             (:,1) double    % number of nodes in the standard electrical spatial grid
        Nzet            (1,1) double    % number of nodes in the total electrical spatial grid (2 points per interface)
        P               (:,:) double    % electrical power (P>0 : output) (W.m^-2)
        PointEq         (1,1) calcPoint % spatial data at equilibrium
        physEf          (:,:) double    % for each voltage, false if non-physical values of quasi Fermi levels are reached
        Rat             (:,:) double    % total Auger net recombination-generation current density (A.m^-2)
        Rrem            (:,:) double    % total radiative recombination current density (i.e. includes only emission) (A.m^-2)
        Rrext           (:,:) double    % total radiative net recombination-generation current density due to exchanges with the opposite component (A.m^-2) - different from Rrt only when photon recycling is included
        Rrt             (:,:) double    % total radiative net recombination-generation current density (i.e. includes both emission and absorption) (A.m^-2)
        Rseries         (:,:) double    % series resistance (Ohm.m^2)
        Rst             (:,:) double    % total Shockley-Read-Hall recombination current density (A.m^-2)
        Rsurf           (:,:) double    % total surface recombination current density (A.m^-2)
        SetPoint        (:,1) calcPoint % spatial data at the voltages set by options
        Sr              (1,2) double    % surface recombination velocity at both boundaries (m.s^-1)
        taum            (:,2) double    % [tauN tauP] on spatial grid (s)
        tauN            (1,:) double    % electron Shockley-Read-Hall lifetime (s)
        tauP            (1,:) double    % hole Shockley-Read-Hall lifetime (s)
        U               (:,1) double    % voltage (V)
        Vfb             (:,1) double    % normalised applied voltage (-)
        Voc             (:,:) double    % open-circuit voltage (V)
        wefleap         (1,1) double {mustBeInRange(wefleap,0,1)}   % weight for leaps in iterative process for quasi Fermi level convergence
        wefn            (1,1) double {mustBeInRange(wefn,0,1)}      % weight for electron quasi-Fermi level in iterative process for quasi Fermi level convergence
        wefp            (1,1) double {mustBeInRange(wefp,0,1)}      % weight for hole quasi-Fermi level in iterative process for quasi Fermi level convergence
        wnext           (1,1) double {mustBeInRange(wnext,0,1)}     % weight for quasi Fermi level and electrostatic potential in global iterative process used to prevent divergence
        ze              (:,1) double    % standard spatial grid - electrical solver (m)
        zea             (:,1) double    % total spatial grid (2 points per interface) - electrical solver (m)
        zel             (:,1) cell      % spatial grid in each layer - electrical solver (m)
        zeq             (:,1) double    % intrinsic Debye length, on spatial grid (m)
    end

    properties (Transient)
        a               (:,1) double    % ratio of consecutive spatial steps (-)
        cp              (1,1) calcPoint % current point gathering spatial data
        cp20            (1,1) calcPoint % point gathering spatial data obtained at (U1(k),0), saved to be used as first guess at (U1(k+1),0)
        n0inter         (:,:) double    % Bose-Einstein distribution related to interband transitions
    end

    methods
        function obj = activeComponent(d,index,w,opt)
            %ACTIVECOMPONENT Constructor of the class. Gather all the necessary
            %properties of the components.
            arguments
                d     (1,:) struct  % basic structure defining the component
                index (1,1) uint32 {mustBeMember(index,[1 2])} % index of the component
                w     (:,1) double  % angular frequency
                opt   (1,1) struct  % options of crescent1D
            end
            global kb e %#ok<GVMIS>

            % call constructor of parent class
            obj@component(d,index,w,opt);

            % check the spatial limit of the electrical system, obtained
            % from the layers in which Sr is provided
            for k=1:length(d)
                if isempty(obj.structure(k).Sr)
                    obj.structure(k).Sr=0;
                end
            end  
            idx=find([obj.structure.Sr]>0);
            obj.NleRange=idx(1):idx(end);
            clear idx
            dj=obj.structure(obj.NleRange);
            
            % gather the necessary data
            obj.Ar=reshape([dj.Ar],[length(dj(1).Ar) length(dj)])';
            obj.Brad=[dj.Brad];
            obj.CTmodel=[dj.CTmodel];
            obj.Cn=[dj.Cn];
            obj.Cp=[dj.Cp];
            obj.Eg=[dj.Eg];
            obj.Nc=[dj.Nc];
            obj.Nv=[dj.Nv];
            obj.Nle=length(obj.NleRange);
            obj.Nze=[dj.Nze]';
            obj.Sr=[dj(1).Sr dj(end).Sr];
            obj.chi=[dj.chi];
            obj.epsr=[dj.epsr];
            for k=1:obj.Nl
                try obj.fracInterband(:,k)=obj.structure(k).fracInterband; catch obj.fracInterband(:,k)=zeros(length(obj.structure(k).eps),1); end %#ok<SEPEX> 
            end
            obj.iwm=[obj.structure.iwm];
            obj.mx=reshape([dj.mx],[length(dj(1).mx) length(dj)])';
            obj.ni=[dj.ni];
            obj.tauN=[dj.tauN];
            obj.tauP=[dj.tauP];
        
            % Set spatial grid for the electrical solver
            obj=SetGridElec(obj);
            
            % Correction coefficient for thermionic emission and tunnelling
            % at interface with different effective masses
            Eb=[-(obj.chi(2:end)-obj.chi(1:end-1))' +(obj.chi(2:end)-obj.chi(1:end-1))'+(obj.Eg(2:end)-obj.Eg(1:end-1))']/kb/obj.T;
            theta=obj.Ar(2:end,:)./obj.Ar(1:end-1,:);
            obj.eta=and(Eb>=0,theta>1).*(theta-(theta-1).*exp(-Eb./(theta-1)))+...
                    and(Eb< 0,theta>1).*1+...
                    and(Eb>=0,theta<1).*theta+...
                    and(Eb< 0,theta<1).*(1-(1-theta).*exp(Eb./((1./theta-1))));
            obj.eta(theta==1)=Inf;
            obj.etat=(Eb>=0).*1+...
                   (Eb< 0).*theta;
            obj.etat(theta==1)=Inf;

            % set voltage grid, and check if the input is correct
            if isempty(opt.U)
                try
                    obj.Nu=opt.Nu(min(length(opt.Nu),obj.index));
                    obj.U=linspace(0,opt.Umax(min(length(opt.Umax),obj.index)),obj.Nu)';
                catch
                    error('activeComponent:voltageNotDefined',"Voltage is not properly defined. Provide either the voltage vector U in column, or the maximum voltage Umax and the voltage vector size Nu.")
                end
            else
                obj.U=opt.U{index};
                obj.Nu=length(obj.U);
            end

            if any(obj.U<0)
                warning('activeComponent:negativeVoltage',"Voltage reaches negative values. This might cause issues in electrical calculations.")
            elseif obj.U(1)~=0
                warning('activeComponent:initVoltageNonzero',"First voltage is not 0. This might cause issues in electrical calculations.")
            elseif max(obj.U)>obj.Eg/e
                error('activeComponent:voltageGreaterThanEg',"Voltage reaches values greater than Eg. This causes error in the calculation of Bose-Einstein distributions.")
            end

            obj.Vfb=e*obj.U/(kb*obj.T);

            % set surface where surface recombinations occur
            switch opt.SurfRecPosition
                case "front"
                    obj.locSurfRec=[true;false];
                case "rear"
                    obj.locSurfRec=[false;true];
                case "both"
                    obj.locSurfRec=[true;true];
                case "none"
                    obj.locSurfRec=[false;false];
            end
            if obj.index == 1
                obj.locSurfRec=flip(obj.locSurfRec);
            end

            % Initialise the current calculation point calcPoint, which
            % gathers the various spatially-varying quantities.
            obj.cp=calcPoint();
            obj.cp.T=obj.T;

            %%% Correction of recombination coefficients due to options

            % Radiative recombination (if true, they are computed using the
            % radiative terms rather than a coefficient B)
            if opt.RadRecElectrolum
                obj.Brad=zeros(1,obj.Nle);
            end

            for k=1:obj.Nle
                ke=obj.NleRange(k);
                % SRH and Auger recombination (temperature dependency)
                if ismember(obj.mat(ke),opt.nonradTDependencyMatList) && ~isnan(opt.nonradTDependencySRH) && ~isnan(opt.nonradTDependencyAuger)
                    obj.tauN(k)=obj.tauN(k)*(obj.T/300)^opt.nonradTDependencySRH;
                    obj.tauP(k)=obj.tauP(k)*(obj.T/300)^opt.nonradTDependencySRH;
                    obj.Cn(k)=obj.Cn(k)*(obj.T/300)^opt.nonradTDependencyAuger;
                    obj.Cp(k)=obj.Cp(k)*(obj.T/300)^opt.nonradTDependencyAuger;
                end
                % Shockley-Read-Hall recombination (composition dependency)
                if ismember(obj.mat(ke),opt.nonradAlloyMatList) && ~isnan(opt.nonradAlloySRH)
                    obj.tauN(k)=opt.nonradAlloySRH;
                    obj.tauP(k)=opt.nonradAlloySRH;
                end
            end
        
            % If true, force the device to operate at the radiative limit
            if opt.rmNonradRec
                obj.tauP=ones(1,obj.Nle); %Note: at the true rad. limit, tau=Inf, but 1 shall be enough
                obj.tauN=ones(1,obj.Nle);
                obj.Cp=zeros(1,obj.Nle);
                obj.Cn=zeros(1,obj.Nle);
            end
            
        end

        function obj = Computen0(obj,w,mu)
            %COMPUTEN0 Compute the effective Bose-Einstein distribution, as
            %a combination of interband and intraband contributions
            arguments
                obj (1,1) activeComponent   % current activeComponent
                w   (:,1) double            % angular frequency (rad.s -1)
                mu  (:,1) double            % chemical potential per layer (eV) 
            end
            Nw=length(w);
            obj.n0=zeros(Nw+1,obj.Nl);
            for k=1:obj.Nl
                iw=[1:obj.iwm(k) obj.iwm(k):Nw];
                obj.n0inter(:,k)=obj.fracInterband(iw,k).*[GBE(w(1:obj.iwm(k)),0,obj.T);GBE(w(obj.iwm(k):Nw),mu(k),obj.T)];
                obj.n0(:,k)=obj.n0inter(:,k)+(1-obj.fracInterband(iw,k)).*obj.n0intra(iw);
            end
        end

        function obj = initElecQuantities(obj,Nu1,Nu2)
            %INITELECQUANTITIES Init physical quantities required for
            %electrical calculations.
            arguments
                obj (1,1) activeComponent   % current activeComponent
                Nu1 (1,1) uint32            % length of voltage grid in component 1
                Nu2 (1,1) uint32            % length of voltage grid in component 2
            end

            obj=initElecQuantities@component(obj,Nu1,Nu2);
            
            obj.cp20=obj.cp;
            
            obj.J=zeros(Nu1,Nu2); 
            obj.Rrt=zeros(Nu1,Nu2);
            obj.Rrext=zeros(Nu1,Nu2);
            obj.Rrem=zeros(Nu1,Nu2);
            obj.Rst=zeros(Nu1,Nu2);
            obj.Rat=zeros(Nu1,Nu2);
            obj.Rsurf=zeros(Nu1,Nu2);
            obj.physEf=true(Nu1,Nu2);
            obj.Rseries=zeros(Nu1,Nu2);
        end

        function [obj,mu20] = saveGuessOuterLoop(obj,ku)
            %SAVEGUESSOUTERLOOP Set the first guess for various electrical
            %quantities before starting the resolution of the drift-diffusion
            %equations. Must be used at the beginning of the outer loop on
            %voltage.
            arguments
                obj (1,1) activeComponent   % current activeComponent
                ku  (1,2) uint32            % voltage index in both components
            end
            obj.cp20=obj.cp;
            if obj.index == 2
                mu20=obj.mu(:,ku(1),ku(2));
            end
        end

        function obj = setGuessOuterLoop(obj,ku1,mu20)
            %SETGUESSOUTERLOOP Set the first guess for various electrical
            %quantities before starting the resolution of the drift-diffusion
            %equations. Must be used at the beginning of the outer loop on
            %voltage.
            arguments
                obj (1,1) activeComponent   % current activeComponent
                ku1 (1,1) uint32            % voltage index in component 1
                mu20(:,1) double            % chemical potential of component 2 at previous voltage of component 1
            end
            obj.cp=obj.cp20;
            if obj.index == 1
                obj.mu(obj.NleRange,ku1,1)=repmat(obj.U(ku1),[obj.Nle,1]);
            else
                obj.mu(:,ku1,1)=mu20;
            end
        end

        function obj = setGuessInnerLoop(obj,ku,opt)
            %SETGUESSINNERLOOP Set the first guess for various electrical
            %quantities before starting the resolution of the drift-diffusion
            %equations. Must be used at the beginning of the inner loop on
            %voltage.
            arguments
                obj (1,1) activeComponent   % current activeComponent
                ku  (1,2) uint32            % voltage index in both components
                opt (1,1) struct            % options of crescent1D
            end

            if ~opt.IterChemPotential % mu=eU if this option is set to false
                obj.mu(obj.NleRange,ku(1),ku(2))=repmat(obj.U(ku(obj.index)),[obj.Nle,1]);
            elseif ku(2)>1 % if ku(2)=1, mu is set by function setGuessOuterLoop
                obj.mu(:,ku(1),ku(2))=obj.mu(:,ku(1),ku(2)-1);
            end
            if ~opt.RadRecElectrolum % if radiative recombinations are already accounted by a radiative coefficient B, must not include them with electroluminescence
                obj.mu(:,ku(1),ku(2))=0;
            end
        end

        function [obj,errMu] = updateChemicalPotential(obj,ku,wmu,iterMu,opt)
            %UPDATECHEMICALPOTENTIAL Updates the chemical potential of each
            %layer using the new data obtained from the drift-diffusion
            %resolution
            arguments
                obj     (1,1) activeComponent   % current activeComponent
                ku      (1,2) uint32            % voltage index in both components
                wmu     (1,1) double            % weight for iterative process on chemical potential of radiation
                iterMu  (1,1) uint32            % number of iterations
                opt     (1,1) struct            % options of crescent1D
            end
            global kb e %#ok<GVMIS>
            muOld=obj.mu(:,ku(1),ku(2));
            for k=1:obj.Nle
                ke=obj.NleRange(k);
                obj.mu(ke,ku(1),ku(2))=kb*obj.T/e*log(mean(obj.cp.n(sum(obj.Nze(1:k-1))+1:sum(obj.Nze(1:k))).*obj.cp.p(sum(obj.Nze(1:k-1))+1:sum(obj.Nze(1:k)))/obj.ni(k)^2));
                if abs(obj.mu(ke,ku(1),ku(2)))<1e-4
                    obj.mu(ke,ku(1),ku(2))=0;
                end
            end
             if ~opt.RadRecElectrolum % if radiative recombinations are already accounted by a radiative coefficient B, must not include them with electroluminescence
                obj.mu(:,ku(1),ku(2))=0;
             end
             errMu=max(abs(obj.mu(:,ku(1),ku(2))-muOld));
             % smoothen only after the second iteration
             if iterMu>=1
                obj.mu(:,ku(1),ku(2))=wmu*obj.mu(:,ku(1),ku(2))+(1-wmu)*muOld;
            end
        end

        function obj = resetWeights(obj,wnext,wefn,wefp,wefleap)
            %RESETWEIGHTS Reset the weighting factors used to smoothen the
            %iterative process.
            arguments
                obj     (1,1) activeComponent   % current activeComponent
                wnext   (1,1) double            
                wefn    (1,1) double
                wefp    (1,1) double
                wefleap (1,1) double
            end
            obj.wnext=wnext;
            obj.wefn=wefn;
            obj.wefp=wefp;
            obj.wefleap=wefleap;
        end

        function [obj,DDstate,q12] = SolveDriftDiff(obj,d,Fext,Fin,w,ku,opt)
            %SOLVEDRIFTDIFF Computes the different radiative generation
            %terms, and solve the drift-diffusion equations at the voltages
            %considered.
            arguments
                obj     (1,1) activeComponent   % current activeComponent
                d       (1,1) component         % opposite component
                Fext    (:,:,:) double          % transmission function with opposite component
                Fin     (:,:,:) double          % transmission function within itself
                w       (:,1) double            % angular frequency (rad.s^-1)
                ku      (1,2) uint32            % voltage index in both components
                opt     (1,1) struct            % options of crescent1D
            end

            % Computing the generation term
            % It is obtained as the sum of the radiative exchanges with the
            % opposite component (Fext) and within itself (Fin), the second
            % being non-zero only if IterChemPotential and PhotonRecycling
            % are set to true.
            [obj.cp.GaExt,q12,~,~,GemExt]=obj.ComputeGenQ(d,Fext,w);
            if opt.IterChemPotential
                [GaInt,~,~,~,GemInt]=obj.ComputeGenQ(obj,Fin,w);
            else
                GaInt=zeros(size(obj.cp.GaExt));
                GemInt=zeros(size(obj.cp.GaExt));
            end
            obj.cp.Ga=(obj.cp.GaExt+GaInt);
            obj.cp.GaEm=(GemExt+GemInt);

            if obj.index == 1
                Nu1=obj.Nu;
            else
                Nu1=d.Nu;
            end

            % Solve the DD equations, possibly several times if there are
            % convergence issues
            DDstate="on";
            while DDstate=="on"
                [obj,DDstate,obj.cp.Ga,obj.cp.GaEm]=DriftDiffHet(obj,obj.Vfb(ku(obj.index)),obj.cp.Ga,ku(1)+Nu1*(ku(2)-1),opt);
            end

            % Stopping computation if not able to converge
            if DDstate=="over - divergence"
                warning(join(['Not able to converge. Moving on to next voltage (component ',num2str(obj.index),', ku1=',num2str(ku(1)),', ku2=',num2str(ku(2)),').']));
                obj.wnext=1;
            else
                obj.cp.U = obj.U(ku(obj.index));
                obj.cp.J = obj.J(ku(1),ku(2));
            end
        end
        
        function [G,varargout] = ComputeGenQ(obj,dsec,radExMat,w)
            % COMPUTEGENQ Computes the spatially resolved generation rate
            % for the current component, coming from the exchanges with 
            % component dsec.
            arguments
                obj     (1,1)   activeComponent   % current activeComponent
                dsec    (1,1)   component         % component with which radiation is exchanged
                radExMat(:,:,:) double            % transmission function between activeComponents (m^-2)
                w       (:,1)   double            % angular frequency (rad.s^-1)
            end
        
            global hb %#ok<GVMIS>
        
            FaAbs=zeros(obj.Nzrt-1,1);
            FaEm=zeros(obj.Nzrt-1,1);
            qabs=0;
            qem=0;
        
            Nw=length(w);
        
            for ks=1:dsec.Nl
                for km=1:obj.Nl
                    % Treatment is different above and below each bandgap,
                    % as emission is then electroluminescent or thermal.
                    iwm=[1:obj.iwm(km) obj.iwm(km):Nw];
                    iws=[1:dsec.iwm(ks) dsec.iwm(ks):Nw];
                    
                    % photon flux density absorbed by the elementary layer
                    % due to emission by the thick layer
                    % Fwa represent the absorption related to interband transition
                    FwAbs=repmat(dsec.n0(:,ks),[1,obj.Nzr(km)-1]).*radExMat(iws,obj.zrange{km},ks);
                    Fwa=FwAbs.*repmat(obj.fracInterband(iws,km),1,length(obj.zrange{km}));
                    if obj.iwm(km)>=dsec.iwm(ks)
                        FaAbs(obj.zrange{km})=FaAbs(obj.zrange{km})'+trapz(w(obj.iwm(km):Nw),Fwa(obj.iwm(km)+1:Nw+1,:),1);
                    else
                        FaAbs(obj.zrange{km})=FaAbs(obj.zrange{km})'+trapz(w(obj.iwm(km):dsec.iwm(ks)),Fwa(obj.iwm(km):dsec.iwm(ks),:),1)+...
                                                                     trapz(w(dsec.iwm(ks):Nw),Fwa(dsec.iwm(ks)+1:Nw+1,:),1);
                    end
                    
                    % photon flux density emitted by interband transition
                    FwEm=repmat(obj.n0inter(obj.iwm(km)+1:Nw+1,km),[1,obj.Nzr(km)-1]).*radExMat(obj.iwm(km):Nw,obj.zrange{km},ks);
                    FaEm(obj.zrange{km})=FaEm(obj.zrange{km})'+trapz(w(obj.iwm(km):Nw),FwEm,1);
                    
                    % photon flux density emitted
                    FwEm=repmat(obj.n0(:,km),[1,obj.Nzr(km)-1]).*radExMat(iwm,obj.zrange{km},ks);
                    
                    % total heat flux density absorbed/emitted (sum of all
                    % exchanges to have total flux)
                    qwAbs=hb*repmat(w([1:dsec.iwm(ks) dsec.iwm(ks):Nw]),[1,obj.Nzr(km)-1]).*FwAbs;
                    qabs=qabs+sum(trapz(w(1:dsec.iwm(ks)),qwAbs(1:dsec.iwm(ks),:),1)+ ...
                                  trapz(w(dsec.iwm(ks):Nw),qwAbs(dsec.iwm(ks)+1:Nw+1,:),1));
                    qwEm=hb*repmat(w([1:obj.iwm(km) obj.iwm(km):Nw]),[1,obj.Nzr(km)-1]).*FwEm;
                    qem=qem+sum(trapz(w(1:obj.iwm(km)),qwEm(1:obj.iwm(km),:),1)+ ...
                                trapz(w(obj.iwm(km):Nw),qwEm(obj.iwm(km)+1:Nw+1,:),1));
                end
            end
           
            % Compute the net radiative generation-recombination rate
            G=zeros(size(obj.zea));
            Gem=zeros(size(obj.zea));
            for km=1:obj.Nle
                kme=obj.NleRange(km);
                G(sum(obj.Nze(1:km-1))+1:sum(obj.Nze(1:km)))=interp1(obj.zr(obj.zrange{kme})/2+obj.zr(obj.zrange{kme}+1)/2,(FaAbs(obj.zrange{kme})-FaEm(obj.zrange{kme}))./obj.dzr(obj.zrange{kme}),sum(obj.tj(1:kme-1))+obj.zel{km}(1:obj.Nze(km)),'linear','extrap');
                Gem(sum(obj.Nze(1:km-1))+1:sum(obj.Nze(1:km)))=interp1(obj.zr(obj.zrange{kme})/2+obj.zr(obj.zrange{kme}+1)/2,FaEm(obj.zrange{kme})./obj.dzr(obj.zrange{kme}),sum(obj.tj(1:kme-1))+obj.zel{km}(1:obj.Nze(km)),'linear','extrap');
            end
        
            % Optional outputs
            varargout{1}=qabs-qem;
            varargout{2}=FaAbs;
            varargout{3}=FaEm;
            varargout{4}=Gem;
        end

        function obj = setRatesConvergence(obj,ku)
            %SETRATESCONVERGENCE Computes the various total (i.e, spatially
            %integrated) recombination-generation rates at given voltage.
            arguments
                obj     (1,1) activeComponent   % current activeComponent
                ku      (1,2) uint32            % voltage index in both components
            end
            global e %#ok<GVMIS>
            obj.Rrt(ku(1),ku(2))=-e*trapz(obj.zea,obj.cp.Ga)+e*trapz(obj.zea,obj.cp.R(:,2));
            obj.Rrext(ku(1),ku(2))=-e*trapz(obj.zea,obj.cp.GaExt)+e*trapz(obj.zea,obj.cp.R(:,2));
            obj.Rrem(ku(1),ku(2))=e*trapz(obj.zea,obj.cp.GaEm)+e*trapz(obj.zea,obj.cp.R(:,2));
            obj.Rst(ku(1),ku(2))=e*trapz(obj.zea,obj.cp.R(:,1));
            obj.Rat(ku(1),ku(2))=e*trapz(obj.zea,obj.cp.R(:,3));
            obj.Rsurf(ku(1),ku(2))=e*sum(obj.cp.Rs);
        end

        function obj = setRatesDivergence(obj,ku)
            %SETRATESDIVERGENCE Ensures that the various rates not
            %computed by the solver are returned as NaN.
            arguments
                obj     (1,1) activeComponent   % current activeComponent
                ku      (1,2) uint32            % voltage index in both components
            end
            obj.J(ku(1),ku(2):end)=NaN;
            obj.Rrt(ku(1),ku(2):end)=NaN;
            obj.Rrext(ku(1),ku(2):end)=NaN;
            obj.Rrem(ku(1),ku(2):end)=NaN;
            obj.Rst(ku(1),ku(2):end)=NaN;
            obj.Rat(ku(1),ku(2):end)=NaN;
            obj.Rsurf(ku(1),ku(2):end)=NaN;
        end

        function obj = ComputeIV(obj,d)
            %COMPUTEIV Compute electrical power input/output of the
            %component, along with Voc and Jsc
            arguments
                obj     (1,1) activeComponent   % current activeComponent
                d       (1,1) component         % opposite component
            end
            % Must flip J to have the same treatment for any component
            % (voltage of the component being the second coordinate)
            if obj.index == 1
                obj.J=obj.J.';
            end
            obj.P=repmat(obj.U',[d.Nu,1]).*obj.J;
            obj.Jsc=obj.J(:,1);
            obj.Voc=nan(d.Nu,1);
            for ku=1:d.Nu
                kuo=find(obj.J(ku,1:obj.Nu-1).*obj.J(ku,2:obj.Nu)<0);
                if length(kuo)==1 && kuo<obj.Nu
                    obj.Voc(ku)=interp1(obj.J(ku,kuo:kuo+1),obj.U(kuo:kuo+1),0);
                end
            end
            % Second flip to obtain the expected order in coordinates
            if obj.index == 1
                obj.J=obj.J.';
                obj.P=obj.P.';
                obj.Jsc=obj.Jsc.';
                obj.Voc=obj.Voc.';
            end
        end

        function obj = setEcv(obj)
            %SETECV Compute the value of Ec and Ev.
            arguments
                obj     (1,1) activeComponent   % current activeComponent
            end
            global kb e %#ok<GVMIS>
            obj.cp.Ec=-kb*obj.T/e*log(obj.cp.n./obj.Ncm./obj.cp.Efnb);
            obj.cp.Ev=kb*obj.T/e*log(obj.cp.p./obj.Nvm./obj.cp.Efpb);
        end

        function obj = setnp(obj)
            %SETNP Compute the value of n and p.
            arguments
                obj     (1,1) activeComponent   % current activeComponent
            end
            Ve=obj.z2za(obj.cp.Vb);
            obj.cp.n=obj.ni(1)*exp(obj.xin).*exp(Ve).*obj.cp.Efnb;
            obj.cp.p=obj.ni(1)*exp(obj.xip).*exp(-Ve).*obj.cp.Efpb;
        end

        function output=z2za(obj,input)
            %Z2ZA Returns input variable on the total spatial mesh (i.e.
            %with two points at each interface).
            arguments
                obj     (1,1) activeComponent   % current activeComponent
                input   (:,1) double            % quantity of interest on initial spatial grid
            end
            for k=1:obj.Nle
                output(sum(obj.Nze(1:k-1))+1:sum(obj.Nze(1:k)))=input(sum(obj.Nze(1:k-1))+2-k:sum(obj.Nze(1:k))+1-k);
            end
            output=output.';
        end
        
        function output=za2z(obj,input)
            %ZA2Z Returns input variable on the regular spatial mesh (i.e.
            %with one point at each interface).
            arguments
                obj     (1,1) activeComponent   % current activeComponent
                input   (:,1) double            % quantity of interest on initial spatial grid
            end
            for k=1:obj.Nle
                output(sum(obj.Nze(1:k-1))+2-k:sum(obj.Nze(1:k))+1-k)=input(sum(obj.Nze(1:k-1))+1:sum(obj.Nze(1:k)));
            end
            output=output.';
        end
    end
end