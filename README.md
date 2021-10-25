# Input data and source code for the paper "European Gas Infrastructure Expansion Planning: An Adaptive Robust Optimization Approach"
### Iegor Riepin, Matthew Schmidt, Luis Baringo, and Felix Müsgens (2021). in review

> The European natural gas market is undergoing fundamental changes, fostering uncertainty regarding both supply and demand. This uncertainty is concentrated in the value of strategic infrastructure investments, e.g., Projects of Common Interest supported by European Union public funds, to safeguard security of supply. This paper addresses this matter by suggesting an adaptive robust optimization framework for the problem of gas infrastructure expansion planning that considers long-term uncertainties. This framework confronts the drawbacks of main-stream methods of incorporating uncertainty in gas market models (i.e., stochastic scenario trees), in which the modeler predefines the probabilities and realization paths of unknown parameters. Our mathematical model endogenously identifies the unfortunate realizations of unknown parameters, and suggests the optimal investments strategies to address them. We use this feature to assess which infrastructure projects are valuable in maintaining system resilience amid cold-winter demand spikes, supply shortages, and budget constraints. The robust solutions point to consistent preferences for specific projects. We find that real-world construction efforts have been focused on the most promising projects from a business perspective. However, we also find that most Projects of Common Interest are unlikely to be realized without financial support, even if they would serve as a hedge against stresses in the European gas system.

### Keywords:
Adaptive robust optimization; Capacity planning; European gas market; Uncertainty

### Links: 
Working paper: http://www.optimization-online.org/DB_FILE/2021/10/8654.pdf

### The code reproduces the benchmarks from the paper 
Note that model output files are uploaded into 'results' folder. The folder contains the results used in the paper in the [GAMS Data eXchange (GDX)](https://www.gams.com/latest/docs/UG_GDX.html) format 

### Citing us

The ARO model used in our paper is free: you can access, modify and share it under the terms of the <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>. This model is shared in the hope that it will be useful for further research on topics of robust optimization and gas transmission network expantion, but without any warranty of merchantability or fitness for a particular purpose. 

If you use the model or its components for your research, we would appreciate it if you
would cite us as follows:
```
@techreport{Riepin2021,
author = {Riepin, Iegor and Schmidt, Matthew and Baringo, Luis and M{\"{u}}sgens, Felix},
institution = {Optimization Online: an eprint site for the optimization community},
keywords = {Adaptive robust optimization,Capacity planning,European gas market,Uncertainty},
number = {8654},
series = {OO Working Paper},
title = {{European Gas Infrastructure Expansion Planning: An Adaptive Robust Optimization Approach}},
url = {http://www.optimization-online.org/DB_HTML/2021/10/8654.html},
year = {2021}
}
```
