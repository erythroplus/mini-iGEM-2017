# mini-iGEM-2017

The repository contains the computer modelling performed by Group 3
of the Imperial College London mini iGEM project in 2016-17.

The project devises a synthetic gene system to sense progesterone-to-estradiol ratios, and produce erythroprotein to
minimimze anaemia during cyclic blood loss, part of the female reproductive cycle.
The system is based on a ratiometric hormone detector circuit with auto-sensitisation, and a recombinase switch.

Released under the CRAPL open-source license, provided in the LICENSE.md file.

Contents:
- Invertase modelling: stochastic modelling of PhiC31 invertase kinetics using a Gillespie algorithm
                       hill_fit.m can be used to fit a sigmoid function to the resulting invertase trajectories
                       
                       
- Circuit modelling:   contains dynamic simulations of the entire circuit presented on our poster
                       run using the hormones_menstruation.m wrapper script

All simulation code was written for MATLAB 2016b.
