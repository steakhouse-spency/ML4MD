# Prediction of the Relaxation Time of Water in Carbon Nanotubes
  
## Authors
Jose Cobena (cobenare@usc.edu) - PhD CHE  
Spencer Ortega (sportega@usc.edu) - MS CS 

## Abstract
Molecular dynamics (MD) simulations are typically used to study glassy systems such as confined water in carbon nanotubes (CNT). The Intermediate scattering function (ISF) which is directly correlated to the structural relaxation time can be computed directly from the trajectory obtained from MD. The advantage of MD data over scattering experiments is the access to very low wave numbers that are usually inaccessible in experiments. However, to cover the whole range of time required to plot the ISF, long simulations (40 ns) sometimes are required. Here, we attempt to predict the relaxation time using supervised learning algorithms on the trajectories obtained from MD. The parameters used in the model are the diameter of the CNT, density of water, temperature and the average displacement of the water molecules up to the divergence time, which is the time where the ISF start to differentiate from each others. This methodology will help to reduce computational time in the study of glassy systems and in particular in the estimation of the relaxation time.

![Image of workflow](https://github.com/spencer-ortega/cs653-final/blob/master/images/workflow.jpg)
  
## Problem
Features (input): nano tube diamater, density (water), temperature, displacement at time t
Label (output): relaxation time (water)

## Methods
LAMMPS - Parallel MD simulator
Supervised Learning - regression, decesion tree, knn, NN, etc. 

## Expected Results
Model that can be used to predict relaxation time of water.  
