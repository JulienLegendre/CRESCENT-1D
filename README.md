# CRESCENT-1D
CRESCENT-1D is a Coupled Radiative and Electrical Solver for Efficient Near-field Thermophotonics in 1D, developed in MATLAB.

## Main information about the solver

This solver can estimate the performance of near-field thermophotonic devices used in a heat engine configuration. Such devices are made of two facing optoelectronic components: a heated light-emitting diode, and a photovoltaic cell kept at ambient temperature. By carefully controlling the voltage applied to each component, the overall system is capable of converting part of the heat delivered to the light-emitting diode into electrical power.

CRESCENT-1D, as its name suggests, couples a radiative and an electrical solver. The radiative solver is based on the fluctuational electrodynamics framework, and provides the spatially- and spectrally-resolved photon transmission coefficient between the two components. Then, the electrical solver simulates charge transport in each component based on the drift-diffusion equations, which allows obtaining their current-voltage characteristic. These two segments of the overall solver are coupled through the radiative heat flux density flowing in the system, which depends both on the photon transmission coefficient and on the local charge concentration (through the photon chemical potential), and is used as an input of the electrical solver.

More information of CRESCENT-1D can be found in the scientific paper introducing it (ref. not yet available) and in Julien Legendre's thesis (https://theses.fr/2023ISAL0094).

## How is the repository organised?

The repository is divided into three main folders:
* **src**, which gathers the main function **crescent1D** along with all the functions necessary to run it (apart from the material properties)
* **material_properties**, in which are gathered the MATLAB functions used to set the **material properties** of the component layers
* **example**, in which we show an **example** of how to call crescent1D

## Licence

The code is made available under licence GNU-GPL 3.0.

## Questions and remarks

In case you have any question of remark, feel free to contact me at: julien.legendre@icfo.eu.
