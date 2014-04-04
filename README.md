emo_2013_viz
============

Codebase for EMO 2013 visualisation paper "Visualising high-dimensional Pareto relationships in two-dimensional scatterplots, Evolutionary Multi-criterion Optimization" by Jonathan E. Fieldsend and Richard M. Everson

Readme for release 1.0 
2/4/2013
Jonathan Fieldsend

This is research quality code, it undoubtedly could benefit from some refactoring, and there is some wasted computation, however I hope the comments (and the published work), make it relatively easy to follow. If you find a bug (and especially if you have a fix for it) please email the author at <J.E.Fieldsend @ exeter.ac.uk>.

There are a number of supporting functions that have been sparsely commented (if commented at all), as they provide some processing which is shared across approaches; the three functions you will want to invoke directly are:

deterministic_compression_and_visualisation_closeness
deterministic_compression_and_visualisation_dominance
deterministic_compression_and_visualisation_koppenyoshida

Please use the help option for each of these methods. 

The first two undertake the visualisations described in:

Fieldsend JE, Everson RM. 

The last undertakes the visualisation described in 

M. Koppen and K. Yoshida. 

*given* an already optimised permutation (Koppen and Yoshida used NSGA-II). If you have not provided an permutation for the non-dominated subset, then spectral seriation can be used instead (and is a built in option).