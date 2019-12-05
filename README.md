# ML for MD simulations
  
## Authors
Jose Cobena (jcobena@usc.edu) - PhD CHE  
Spencer Ortega (sportega@usc.edu) - MS CS 

## Abstract
Molecular dynamics (MD) simulations are typically used to compute the Intermediate scattering function (ISF) which is directly correlated to the structural relaxation time. The advantage of MD data over scattering experiments is the access to very low wave numbers that are usually inaccessible in experiments. However, to achieve accurate results, long simulations (40 ns) sometimes are required. Here, we attempt to predict the relaxation time using supervised learning algorithms on the trajectories obtained from MD. This methodology will help to reduce computational time in the study of glassy systems.
  
## Problem
Features (input): nano tube diamater, water density, temperature  
Label (output): relaxation time.  

## Methods
LAMMPS - Parallel MD simulator
Supervised Learning - regression, decesion tree, knn, NN, etc. 

## Expected Results
Algorithm/Model that can be used for to predict this property.  
