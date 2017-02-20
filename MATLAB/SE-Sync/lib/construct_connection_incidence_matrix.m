function [A_,Omega_] = construct_connection_incidence_matrix(measurements)
% [A_,Omega_] = construct_connection_incidence_matrix(measurements)
%
% This function computes and returns the oriented *connection* incidence
% matrix of the underlying directed graph of measurements.
% The k-th block column of this matrix is defined as:
%
%   A[k,b_ij,vi,vj] = -vi, if edge k *leaves* node i
%                     +vj, if edge k *enters* node j,
%                     zeros(...), otherwise.
%
% (see eq. (?) in the paper).
% 
% This also provides the *connection* (diagonal) weight matrix
% that has the *matrix* weights for all the measurements:
% 
%   Omega = blkdiag( Om_k )
% 
% (see eq. (?) in the paper).
% 
% Copyright (C) 2017 by Jesus Briales

% get main parameters defining the connection graph
d = length(measurements.t{1}); % d = dimension of SE(d)
n = max(max(measurements.edges));  % N = number of nodes in the pose graph
m = size(measurements.edges,1); % M = number of edges in the pose graph

% define parameters for blocks in the matrix
d_  = d+1;  % True dimension of the SE(d) block (+1 for translation)
Os  = [zeros(1,d),1]; % convenient vector for defined blocks
Id_ = eye(d_); % convenient block

% preallocate indices for the sparse blocks:
% cell array with lists of scalar indeces corresponding to each block-index
NNZ = d_^2*m;
c_blkidx = mat2cell( 1:m*d_^2, 1, d_^2*ones(1,m) );
% preallocate arrays for storage of sparse description
[rows_in,cols_in,rows_out,cols_out,vals_in,vals_out] = deal(zeros(1,NNZ));
c_Om = cell(1,m);

%Iterate over the measurements in the pose graph
for k = 1:m
   
    %EXTRACT MEASUREMENT DATA
    i = measurements.edges(k, 1);  %The node that this edge leaves
    j = measurements.edges(k, 2);  %The node that this edge enters
    
    Rij = measurements.R{k};      %Rotation matrix for this observation
    kap = measurements.kappa{k};  %Precision for this rotational observation    
    tij = measurements.t{k};      %Translation corresponding to this observation
    tau = measurements.tau{k};    %Precision for this translational observation
    % build SE(d) block: this is where we set specific values for PGO case
    Tij = [Rij,tij ; Os];
    vi = Tij;
    vj = Id_; % Note we treat it as a complete block for generality
    om = spdiags([kap*ones(d,1);tau],0,d_,d_);
    
    
    %PROCESS MEASUREMENT DATA TO STORE IT IN SPARSE MATRICES
    % Note the blkidx for the column is k in both cases,
    % and the same list of column indeces could be duplicated
    
    % Set the k-th d_xd_ SE(d) block leaving node i (NOTE: NEGATIVE)
    [r, c, v_Tij] = rcvize_matrix(Tij, i,k);
    rows_in(c_blkidx{k})  = r;
    cols_in(c_blkidx{k})  = c;
    vals_in(c_blkidx{k})  = -v_Tij;
    
    % Set the k-th d_xd_ SE(d) block entering node j (NOTE: POSITIVE)
    [r, c, v_Id_] = rcvize_matrix(Id_, j,k);
    rows_out(c_blkidx{k}) = r;
    cols_out(c_blkidx{k}) = c;
    vals_out(c_blkidx{k}) = v_Id_;
    
    % Set the k-th d_xd_ weight block
    c_Om{k} = om; % om is already sparse, otherwise sparse(om)
    
end


%% Create the sparse matrices

% concatenate all data for connection incidence: input and output nodes
rows = [rows_in, rows_out];
cols = [cols_in, cols_out];
vals = [vals_in, vals_out];
A_ = sparse(rows, cols, vals, d_*n, d_*m);

% save all weight blocks into a block diagonal matrix
Omega_ = blkdiag(c_Om{:});

end