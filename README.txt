Synchronization over the Special Euclidean group SE(d),
also commonly known as Pose Graph Optimization (PGO),
is a common problem arising in the context of 2D and 3D geometric estimation
(e.g. pose-graph SLAM and camera motion estimation).

Certifiably correct algorithms are those that, with high probability,
are capable of attaining global optimality and certifying having done so, a posteriori.

This repository features two different certifiably correct algorithms
for synchronization over the special Euclidean group:

- Cartan-Sync, our contribution,
  a novel algorithm based upon the application of the Riemannian staircase
  directly on the domain of Pose Graph Optimization (list of poses).
  
- SE-Sync, proposed by Rosen et al., which we also dub Stiefel-Sync,
  that marginalizes PGO and then applies the Riemannian staircase
  in the resulting Rotation Synchronization problem.

The present code extends that of github.com/david-m-rosen/SE-Sync
corresponding to the reference work providing a simple interface
to call and compare both methods.

We are making this software freely available in the hope that it will be useful to others.
Note however Cartan-Sync is ongoing work, so important changes may still happen in the code.
Also if you find any issues, let us know!

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

If you use SE-Sync (Stiefel-Sync) instead, please refer to the relevant works
pointed in https://github.com/david-m-rosen/SE-Sync.

If you use the MATLAB implementations, please also cite the following reference for the Manopt toolbox, which provides the MATLAB implementation used for RTR:

@article{Boumal2014Manopt,
  title={{Manopt}, a {MATLAB} Toolbox for Optimization on Manifolds.},
  author={Boumal, N. and Mishra, B. and Absil, P.-A. and Sepulchre, R.},
  journal={Journal of Machine Learning Research},
  volume={15},
  number={1},
  pages={1455--1459},
  year={2014}
}

If you use the C++ implementations, please also cite the following reference for the ROPTLIB library, which provides the C++ implementation of RTR employed within the C++ library:

@techreport{Huang16ROPTLIB,
title = {{ROPTLIB}: An Object-Oriented {C++} Library for Optimization on {Riemannian} Manifolds},
author = {Huang, W. and Absil, P.-A. and Gallivan, K.A. and Hand, P.},
institution = {Florida State University},
number = {FSU16-14},
year = {2016},
}


==== Copyright and License ====

The MATLAB and C++ implementations of Cartan-Sync contained herein are copyright (C) 2017 by Jesus Briales, and are distributed under the terms of the GNU General Public License (GPL) version 3 (or later).  Please see the files LICENSE.txt and COPYING.txt for more information.

Contact: jbriales@uma.es


==== Getting Started ====

MATLAB:

To use the MATLAB implementation of Cartan-Sync (and SE-Sync), simply place the 'MATLAB' folder in any convenient (permanent) location, and then run the script MATLAB/setup.m.  For a minimal working example, see MATLAB/main.m

C++:

The C++ implementation of Cartan-Sync (and SE-Sync) can be built and exported as a CMake project.
For a minimal working example, see C++/examples/main, which provides a simple command-line utility for processing .g2o files.
For a fast deployment just go to the root folder and run:
```
cd C++
cmake ..
make
cd ..
C++/bin/se-sync data/kitti_00.g2o CartanSync
```
The general syntax for the main application is
`se-sync <g2oFile> <Algorithm>`,
where `Algorithm` can be CartanSync or StiefelSync.

