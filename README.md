# STNMF-SNN

This code is to demonstrate dissecting cascade computational components of spiking neural components with spike-triggered non-negative matrix factorization.

Accompanying publication:
“Dissecting cascade computational components in spiking neural networks”
by Shanshan Jia, Dajun Xing, Zhaofei Yu, Jian K. Liu, PLoS Comput Biol (2021).


### Usage:
run demo.m in Matlab

### SNN model:
The program simulates a 3-layer SNN model to generate Data.mat and STE.mat. The program can just be run directly as is, using the supplied code. <br>
One should manually adjust iglist to order the inferred subunits as modeled subunits. <br>
(6 subunits in layer 1, 2 cells in layer 2, 1 cell in layer 3, as in Fig. 3 ).

### STNMF Analysis:
Once the simulated data are generated, STNMF will be conducted to generate the following outputs:

The program generates 4 figures:<br>
1. The STA of layer 3, displayed as a 2D color map.<br>
2. The current best set of modules (of layer 1) and the subunits of the model cell, displayed as 2D color maps.<br>
3. The spike train of inferred and modeled subunits.<br>
4. The spike correlation between inferred and modeled spike train.



### License issues:
The code in this package is distributed under the GNU General Public License.<br>

There are some utility codes taken from various sources and used for modeling and analysis.

This code is a further extension based on STNMFanalysis.  STNMFanalysis, website: https://github.com/jiankliu/STNMFanalysis.<br>

The semi-NMF algorithm is based on modified code from the NMF MATLAB Toolbox by Yifeng Li, Alioune Ngom (https://sites.google.com/site/nmftool/, citation: Y. Li and A. Ngom, The non-negative matrix factorization toolbox for biological data mining, BMC Source Code for Biology and Medicine, vol 8, pp. 10), which is also distributed under the GNU General Public License. See the accompanying license file and license notes in the source code files.<br>
