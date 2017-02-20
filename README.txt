Cartan-Sync is a certifiably correct algorithm for synchronization over the special Euclidean group,
a common problem arising in the context of 2D and 3D geometric estimation (e.g. pose-graph SLAM and camera motion estimation).
This method presents a different approach from that presented in
"Rosen et al. SE-Sync: A Certifiably Correct Algorithm for Synchronization over the Special Euclidean Group",
arXiv Preprint arXiv:1612.07386.

The present code extends that of github.com/david-m-rosen/SE-Sync
corresponding to the reference work providing a simple interface
to call and compare both methods.

We are making this software freely available in the hope that it will be useful to others. 
If you use Cartan-Sync in your own work, please cite our paper: 

@article{briales17CartanSync,
title = {{Cartan-Sync}: Fast and Global {SE(d)-Synchronization}},
author = {Briales, J. and Gonzalez-Jimenez, J.},
booktitle = {Submitted to Robotics and Automation Letters (RA-L)},
}

and the following paper of Absil et al., which describes the Riemannian trust-region (RTR) method that Cartan-Sync employs:

@article{Absil2007Trust,
title = {Trust-Region Methods on {Riemannian} Manifolds},
author = {Absil, P.-A. and Baker, C.G. and Gallivan, K.A.},
journal = {Found.\ Comput.\ Math.},
volume = {7},
number = {3},
pages = {303--330},
year = {2007},
month = jul,
}

If you use the MATLAB implementation of Cartan-Sync, please also cite the following reference for the Manopt toolbox, which provides the MATLAB implementation we use for RTR:

@article{Boumal2014Manopt,
  title={{Manopt}, a {MATLAB} Toolbox for Optimization on Manifolds.},
  author={Boumal, N. and Mishra, B. and Absil, P.-A. and Sepulchre, R.},
  journal={Journal of Machine Learning Research},
  volume={15},
  number={1},
  pages={1455--1459},
  year={2014}
}

==== Copyright and License ====

The MATLAB implementations of Cartan-Sync contained herein are copyright (C) 2017 by Jesus Briales, and are distributed under the terms of the GNU General Public License (GPL) version 3 (or later).  Please see the files LICENSE.txt and COPYING.txt for more information.

Contact: jbriales@uma.es


==== Getting Started ====

To use the MATLAB implementation of Cartan-Sync, simply place the 'MATLAB' folder in any convenient (permanent) location, and then run the script MATLAB/setup.m.  For a minimal working example, see MATLAB/main.m
