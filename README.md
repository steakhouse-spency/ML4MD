# Prediction of Relaxation Time of Water in Carbon Nanotubes
  
## Authors
Jose Cobena (jcobena@usc.edu) - PhD CHE  
Spencer Ortega (sportega@usc.edu) - MS CS 

## Abstract
Molecular dynamics (MD) simulations are typically used to compute the Intermediate scattering function (ISF) which is directly correlated to the structural relaxation time. The advantage of MD data over scattering experiments is the access to very low wave numbers that are usually inaccessible in experiments. However, to achieve accurate results, long simulations (40 ns) sometimes are required. Here, we attempt to predict the relaxation time using supervised learning algorithms on the trajectories obtained from MD. This methodology will help to reduce computational time in the study of glassy systems.

![Image of workflow](https://github.com/spencer-ortega/cs653-final/blob/master/images/workflow.jpg)
  
## Problem
Features (input): nano tube diamater, density (water), temperature, displacement at time t
Label (output): relaxation time (water)

## Methods
LAMMPS - Parallel MD simulator
Supervised Learning - regression, decesion tree, knn, NN, etc. 

## Expected Results
Model that can be used to predict relaxation time of water.  
