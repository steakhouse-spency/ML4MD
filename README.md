# Supervised Learning Techniques to Predict the Relaxation Time of Water confined in Nanotubes
  
## Authors
Jose Cobena (cobenare@usc.edu) - PhD CHE  
Spencer Ortega (sportega@usc.edu) - MS CS 

## Abstract
Molecular dynamics (MD) simulations are typically used to study systems that are difficult to analyze in experiments. One MD system that has gained the attention of researchers in recent years is the study of confined water in Nanotubes (NT). Among the properties that can be computed to characterize such a system is the relaxation time. The Intermediate scattering function (ISF), which is directly correlated to the structural relaxation time, can be computed directly from the trajectory obtained from MD by fitting the ISF curves to the Kohlrausch–Williams–Watts (KWW) function. The advantage of MD data over scattering experiments is the access to very low wavenumbers that are usually inaccessible in experiments. However, to cover the whole range of time required to plot the ISF, long simulations (40 ns) are required. Here, we attempt to predict the relaxation time of confined water in NTs by computing the ISF for multiple NTs, all varying by features that globally describe the system (tube material, tube diameter, temperature), and feeding this labeled data to some supervised learning algorithm. 

![Image of workflow](https://github.com/spencer-ortega/cs653-final/blob/master/images/workflow.jpg)

## Methodology

1. Create and fill NT Structures
    What kind?
      C, SiC, ...
      270 - 320 by 20
      diameter
    Wanted to figure out an open-source way to create structures
      used VMD
        creates conects for tube
    
    Optimally fill tube with water
      Jose had a python script he would use for his research
        But hardcoded for 1 type of tube
      optimized and generalized his script
        tons of redundant I/O
        takes in arguments for NT properties
        
    At this point: Able to create nanotube with water
      
      
      

2. Choose and validate force fields

  We initially wanted to use norebo
    Ice4p?
    All NT atoms are fixed
      Water interaction was terrible using open-source tubes
  
  While trying to tweak and test our no-rebo approach, we simultaneously started testing rebo
    other model name EW
    initially thought it was norebos fault
    rebo allows the tube to move as well
      water still was not acting right
    Added anchors to stabilize the tube
      still not right
      
  Since both types of simulations were bad, we then focused on the tube structure
    decided to use an empty Carbon NT from Joses previous research to see if out structure was the issue
      ran perfectly for both norebo and rebo
      the way we created tubes via VMD was the issue
        still dont know exactly why
      
      
  

3. Simulate confined water in NT

4. Compute relaxation time

5. Train and test ML

## Challenges

## Future Work




## Problem
Features (input): nano tube diamater, density (water), temperature, average displacement at time t (divergence time), intial cage-rearrangement time.
Label (output): relaxation time (water)

## Methods
LAMMPS - Parallel MD simulator
Supervised Learning - regression, decesion tree, knn, NN, etc. 

## Expected Results
Model that can be used to predict relaxation time of water.  

## References
[1] Cobena-Reyes, J., Kalia, R. K., & Sahimi, M. (2018). Complex behavior of ordered and icelike water in carbon nanotubes near its bulk boiling point. The journal of physical chemistry letters, 9(16), 4746-4752.
