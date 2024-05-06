classdef calcPoint
    %CALCPOINT Electrical data at a given voltage, used in the function
    %crescent1D and related to a object of class activeComponent.

    properties
        D       (:,2) double % diffusion coefficient (m^2.s^-1)
        Ec      (:,1) double % lowest energy in the conduction band (eV)
        Ev      (:,1) double % highest energy in the valence band (eV)
        Efnb    (:,1) double % Slotboom variable for electrons (-)
        Efpb    (:,1) double % Slotboom variable for holes (-)
        J       (1,1) double % total current density (A.m^-2)
        Ga      (:,1) double % spatially-resolved net radiative generation-recombination rate (m^-3.s^-1)
        GaEm    (:,1) double % spatially-resolved radiative recombination rate (related to emission) (m^-3.s^-1)
        GaExt   (:,1) double % spatially-resolved net radiative generation-recombination rate exchanged with the opposite component (m^-3.s^-1)
        Jn      (:,1) double % electron current (A.m^-2)
        Jp      (:,1) double % hole current (A.m^-2)
        n       (:,1) double % electron density (m^-3)
        p       (:,1) double % hole density (m^-3)
        R       (:,3) double % spatially-resolved recombination rate (column 1: SRH, 2: radiative (if computed with B), 3: Auger) (m^-3.s^-1) 
        Rs      (:,1) double % surface recombination rate (m^-2.s^-1)
        T       (1,1) double % temperature (K)
        U       (1,1) double % voltage (V)
        Vb      (:,1) double % normalised electrostatic potential (-)
    end

    properties (Dependent)
        Efn     (:,1) double % electron quasi-Fermi level (eV)
        Efp     (:,1) double % hole quasi-Fermi level (eV)
        V       (:,1) double % electrostatic potential (V)
    end

    methods
        function cp = calcPoint
            %CALCPOINT Constructor
        end

        function value = get.Efn(cp)
            global kb e %#ok<GVMIS>
            value=kb*cp.T/e*log(cp.Efnb);
        end
        function value = get.Efp(cp)
            global kb e %#ok<GVMIS>
            value=-kb*cp.T/e*log(cp.Efpb);
        end
        function value = get.V(cp)
            global kb e %#ok<GVMIS>
            value=kb*cp.T*cp.Vb/e;
        end

    end
end