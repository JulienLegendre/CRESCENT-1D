classdef component
    %COMPONENT Component used for radiative heat engine. Class used for
    %the 1D simulation of near-field radiative heat transfer and
    %drift-diffusion simulations performed in the function crescent1D. If
    %the component is active (eg LED, PV cell), use the child class
    %activeComponent instead.

    properties (GetAccess='public',SetAccess='immutable')
        index (1,1) uint32                                                                      % index of the component in the overall device
    end

    properties  
        DopSign     (1,:) double {mustBeMemberOrNaN(DopSign,[-1 1])} = NaN                      % doping sign (-1: donor, +1: acceptor)                                                                       
        dzr         (:,1) double                                                                % spatial step - radiative mesh (m)
        eps         (:,:) double                                                                % dielectric function (-)
        hasMirror   (1,1) string {mustBeMember(hasMirror,["no","imperfect","perfect"])} = "no"  
        isSemiInf   (1,1) logical = false
        iwm         (1,:) uint32                                                                % first index in frequency mesh where interband transitions occur
        mat         (1,:) string                                                                % material
        mu          (:,:,:) double                                                              % chemical potential of radiation, equal to 0 for passive component (eV)
        Ndop        (1,:) double {mustBeNonnegativeOrNaN} = NaN                                 % doping level, unsigned (m^-3)
        Nl          (1,1) uint32                                                                % number of layers
        Nu          (1,1) uint32 {mustBePositive} = 1                                           % size of the voltage grid
        Nzr         (:,1) uint32                                                                % size of the standard radiative spatial grid
        Nzrt        (1,1) uint32                                                                % size of the total radiative spatial grid (two points per interface)
        n0intra     (:,1) double {mustBePositive} = 1                                           % Bose-Einstein distribution (thermal radiation)
        structure   (1,:) struct                                                                % layer-by-layer structure
        T           (1,1) double                                                                % temperature (K)
        tj          (1,:) double {mustBeNonnegative}                                            % thickness (m)
        x           (1,:) double {mustBeInRangeOrNaN(x,0,1)} = NaN                              % alloy fraction (-)
        zr          (:,1) double                                                                % standard radiative spatial grid (m)
        zra         (:,1) double                                                                % total spatial grid (two points per interface) - radiative solver (m)          
        zrl         (:,1) cell                                                                  % spatial grid in each layer  - radiative solver (m)
        zrange      (:,1) cell                                                                  % points of the spatial grid related to each layer - radiative solver
    end

    properties (Transient)        
        n0          (:,:) double                                                                % Total Bose-Einstein distribution
    end

    methods
        function obj = component(d,index,w,opt)
            %COMPONENT Constructor of the class. Gather all the necessary
            %properties of the components.
            arguments
                d     (1,:) struct  % basic structure defining the component
                index (1,1) uint32 {mustBeMember(index,[1 2])} % index of the component
                w     (:,1) double  % angular frequency
                opt   (1,1) struct  % options of crescent1D
            end

            % set the component index
            obj.index=index;

            % check if the component is semi-infinite
            tjext=[d([1 end]).tj];
            if tjext(obj.index)==Inf
                obj.isSemiInf=true;
            end
            if any(isinf([[d(2:end-1).tj] tjext(3-obj.index)])) && length(d)>1
                error('Cannot have intermediate layer with infinite thickness.');
            end
            if opt.m(index)==1 && obj.isSemiInf
                error('Cannot have mirror on the back of semi-infinite component.');
            end 

            % load the necessary material properties
            for k=1:length(d)
                di(k)=d(k).fmat(d(k),w); %#ok<AGROW> 
            end
            % add the optional mirror to the component
            if opt.m(obj.index)==1 
                if index == 1
                    di=flip(di);
                end
                di(length(di)+1).eps=epsMirror(opt.MirrorMaterial,w);
                di(end).iwm=length(w);
                di(end).mat=opt.MirrorMaterial;
                di(end).Nzr=opt.MirrorMeshsize;
                di(end).T=di(1).T;
                di(end).tj=opt.MirrorThickness;
                if di(end).tj == Inf
                    obj.isSemiInf=true;
                end
                if index == 1
                    di=flip(di);
                end
                if opt.MirrorMaterial~="Perfect"
                    obj.hasMirror="imperfect";
                else
                    obj.hasMirror="perfect";
                end
            else
                obj.hasMirror="no";
            end

            % gather all properties in component
            obj.Nl=length(di);
            obj.Nzr=[di.Nzr]';
            obj.T=di(1).T;
            obj.eps=reshape([di.eps],[length(di(1).eps) length(di)]);
            obj.iwm=length(di(1).eps)*ones(1,length(obj.Nzr));
            obj.mat=[di.mat];
            obj.n0intra=GBE(w,0,obj.T);
            obj.tj=[di.tj];

            % set spatial grid for the radiative solver
            obj=SetGridRad(obj);

            % load additional data, if applicable
            for k=1:obj.Nl
                try obj.x(k)=di(k).x;               catch obj.x(k)=NaN;       end %#ok<SEPEX> 
                try obj.DopSign(k)=di(k).DopSign;   catch obj.DopSign(k)=NaN; end %#ok<SEPEX> 
                try obj.Ndop(k)=di(k).Ndop;         catch obj.Ndop(k)=NaN;    end %#ok<SEPEX> 
            end

            % save layer-by-layer structure
            obj.structure=di;

        end

        function obj = Computen0(obj,w,~)
            %COMPUTEN0 Compute Bose-Einstein distribution
            arguments
                obj (1,1) component         % current activeComponent
                w   (:,1) double            % angular frequency (rad.s -1)
                ~
            end
            Nw=length(w);
            obj.n0=zeros(Nw+1,obj.Nl);
            for k=1:obj.Nl
                iw=[1:obj.iwm(k) obj.iwm(k):Nw];
                obj.n0(:,k)=obj.n0intra(iw);
            end
        end

        function obj = initElecQuantities(obj,Nu1,Nu2)
            %INITELECQUANTITIES Init physical quantities required for
            %electrical calculations.
            arguments
                obj (1,1) component         % current activeComponent
                Nu1 (1,1) uint32            % length of voltage grid in component 1
                Nu2 (1,1) uint32            % length of voltage grid in component 2
            end
            obj.mu = zeros(obj.Nl,Nu1,Nu2);
        end

        function obj = saveGuessOuterLoop(obj,~)
            %SAVEGUESSOUTERLOOP Set the first guess for various electrical
            %quantities before starting the resolution of the drift-diffusion
            %equations. Must be used at the beginning of the outer loop on
            %voltage.
        end

        function obj = setGuessOuterLoop(obj,~,~)
            %SETGUESSOUTERLOOP Set the first guess for various electrical
            %quantities before starting the resolution of the drift-diffusion
            %equations. Must be used at the beginning of the outer loop on
            %voltage.
        end

        function obj = setGuessInnerLoop(obj,~,~)
            %SETGUESSINNERLOOP Set the first guess for various electrical
            %quantities before starting the resolution of the drift-diffusion
            %equations. Must be used at the beginning of the inner loop on
            %voltage.
        end

        function [obj,errMu] = updateChemicalPotential(obj,~,~,~,~)
            %UPDATECHEMICALPOTENTIAL Updates the chemical potential of each
            %layer using the new data obtained from the drift-diffusion
            %resolution
            errMu=0;
        end

        function obj = resetWeights(obj,~,~,~,~)
            %RESETWEIGHTS Reset the weighting factors used to smoothen the
            %iterative process.
        end

        function obj = setRatesConvergence(obj,~)
            %SETRATESCONVERGENCE Computes the various total (i.e, spatially
            %integrated) recombination-generation rates at given voltage.
        end

        function obj = setRatesDivergence(obj,~)
            %SETRATESDIVERGENCE Ensures that the quantities not
            %computed by the solver are returned as NaN.
        end

        function obj = ComputeIV(obj,~)
            %COMPUTEIV Compute electrical power input/output of the
            %component, along with Voc and Jsc
        end
    end
end

function mustBeNonnegativeOrNaN(value)
    %MUSTBEMEMBERORNAN value must be between the values given in S,
    %or NaN
    if ~all(value>=0 | isnan(value))
        eidType = 'mustBeNonnegativeOrNaN:notNonnegativeOrNaN';
        msgType = "Value must be greater than or equal to 0, or NaN.";
        error(eidType,msgType)
    end
end

function mustBeMemberOrNaN(value,S)
    %MUSTBEMEMBERORNAN value must be between the values given in S,
    %or NaN
    if ~all((value==S(1) | value==S(2)) | isnan(value))
        eidType = 'mustBeMemberOrNaN:notMemberOrNaN';
        msgType = "Value must be a member of the provided set, or NaN.";
        error(eidType,msgType)
    end
end

function mustBeInRangeOrNaN(value,lower,upper)
    %MUSTBEMEMBERORNAN value must be between the values given in S,
    %or NaN
    if ~ all((value>=lower & value<=upper) | isnan(value))
        eidType = 'mustBeInRangeOrNaN:notInRangeOrNaN';
        msgType = "Value must be greater than or equal to "+S(1)+" and less than or equal to "+S(2)+", or NaN.";
        error(eidType,msgType)
    end
end