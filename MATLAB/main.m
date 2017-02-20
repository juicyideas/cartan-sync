% This file contains a minimum working example demonstrating the use of the
% MATLAB distribution of Cartan-Sync, a certifiably correct algorithm for 
% synchronization over the special Euclidean group
% 
% The code provides a simple interface to allow calling also the similar
% approach SE-Sync [1] for comparison.
%
% Copyright (C) 2017 by Jesus Briales
% Modified from homonym file by David Rosen in github.com/david-m-rosen/SE-Sync

% [1] Rosen, D. M., Carlone, L., Bandeira, A. S., & Leonard, J. J. (2016).
%     SE-Sync: A Certifiably Correct Algorithm for Synchronization over the Special Euclidean Group.
%     arXiv Preprint arXiv:1612.07386, 1–49. Retrieved from https://arxiv.org/abs/1612.07386

function main( approach )

% Choose approach to solve SE(d)-synchronization
if nargin < 1
  approach = 'Cartan'; % Our approach (by default)
%   approach = 'Margin'; % Reference approach [1]
end

%% Import SE-Sync
run setup.m;


%% Select dataset to run
data_dir = '../data/';  % Relative path to directory containing example datasets

% 3D datasets
sphere2500 = 'sphere2500';
spherenoise = 'sphere_bignoise_vertex3';
torus = 'torus3D';
grid = 'grid3D';
garage = 'parking-garage';
cubicle = 'cubicle';
rim = 'rim';

% 2D datasets
CSAIL = 'CSAIL';
manhattan = 'manhattan';
city10000 = 'city10000';
intel = 'intel';
ais = 'ais2klinik';

% KITTI datasets:
% 'kitti_00'
% 'kitti_02'
% 'kitti_05'
% 'kitti_06'
% 'kitti_07'
% 'kitti_08'
% 'kitti_09'

% Pick the dataset to run here
file = 'kitti_00';
g2o_file = fullfile( data_dir, strcat(file,'.g2o') );

%% Read in .g2o file
tic();
fprintf('Loading file: %s ...\n', g2o_file);
measurements = load_g2o_data(g2o_file);
t = toc();
num_poses = max(max(measurements.edges));
num_measurements = length(measurements.kappa);
d = length(measurements.t{1});
fprintf('Processed input file %s in %g seconds\n', g2o_file, t);
fprintf('Number of poses: %d\n', num_poses);
fprintf('Number of measurements: %d\n', num_measurements);

%% Set Manopt options (if desired)
Manopt_opts.tolgradnorm = 1e-2;  % Stopping tolerance for norm of Riemannian gradient
Manopt_opts.rel_func_tol = 1e-5;  % Additional stopping criterion for Manopt: stop if the relative function decrease between two successive accepted iterates is less than this value
Manopt_opts.miniter = 1;  % Minimum number of outer iterations (i.e. accepted update steps) to perform
Manopt_opts.maxiter = 300;  % Maximum number of outer iterations (i.e. accepted update steps) to perform
Manopt_opts.maxinner = 150;  % Maximum number of iterations for the conjugate-gradient method used to compute approximate Newton steps

%% Set SE-Sync options (if desired)
SE_Sync_opts.r0 = 5;  % Initial maximum-rank parameter at which to start the Riemannian Staircase
SE_Sync_opts.rmax = 10;  % Maximum maximum-rank parameter at which to terminate the Riemannian Staircase
SE_Sync_opts.eig_comp_rel_tol = 1e-4;  % Relative tolerance for the minimum-eigenvalue computation used to test for second-order optimality with MATLAB's eigs() function
SE_Sync_opts.min_eig_lower_bound = -1e-3;  % Minimum eigenvalue threshold for accepting a maxtrix as numerically positive-semidefinite
SE_Sync_opts.Cholesky = true;  % Select whether to use Cholesky or QR decomposition to compute orthogonal projections (in SE-Sync [1])
SE_Sync_opts.Precon = true; % Use preconditioning, if available

% SE_Sync_opts.init = 'origin'; % Initialize at conventional origin point in the manifold
% SE_Sync_opts.init = 'random'; % Random initialization in the manifold
SE_Sync_opts.init = 'chordal';  % Use chordal initialization for rotations

% Setup SE-Sync problem (object)
switch approach
  case 'Cartan'
    sync_problem = SE_Sync_Cartan(measurements, Manopt_opts, SE_Sync_opts);
  case 'Margin'
    sync_problem = SE_Sync_Margin(measurements, Manopt_opts, SE_Sync_opts);
  otherwise
    error('Unknown approach')
end
% Run whole solver pipeline for SE-Sync
[SDPval, Xopt, xhat, Fxhat, SE_Sync_info] = run( sync_problem );

% Plot resulting solution
plot_loop_closures = true;
if plot_loop_closures
    plot_poses(xhat.t, xhat.R, measurements.edges, '-b', .25);
else
    plot_poses(xhat.t, xhat.R);
end
axis tight;
% view(90, -90);  % For plotting 3D but nearly-planar datasets

