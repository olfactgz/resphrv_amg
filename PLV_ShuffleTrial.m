function [r, perm_plv, permz] = PLV_ShuffleTrial( x, y, n_perm, b_ind)
% Phase locking value between x and y across trials.
%
% Usage
%   [r, ~, permz] = PLV_ShuffleTrial( x, y, nperm);
%   r = PLV_ShuffleTrial( x, y, nperm, b_ind);
%
% Input
%   x, freq x time x trial analytic time series
%   y, freq x time x trial analytic time series
%   nperm, number of permutations, a scalar
%   boot_ind, n_trl x n_boot bootstrap index, see the example below
%
% Output
%   r, a data structure with the following fields
%       z, z score, useful for group-level random effect analysis
%       plv, real plv (complex number) if no bootstapping requested
%       absplv, average plv across bootstraps
%       n, number of trials
%
%   perm_plv, permuted plv (with the amplitude of the real plv being the first).
%       permutation (nperm + 1) x freq x time
%
%
% Example of boot_ind
%   % total number of trials
%   n_total = 200;
%   % Bootstrapping: 100 bootstraps, sample size of 60
%   n_boot = 100;
%   b_sz = 60;
%   b_ind = arrayfun( @(x) randsample(n_total, b_sz, 0), 1:n_boot, 'un', 0);
%   b_ind = cell2mat( b_ind);
%
%
% GZ

if nargin < 3
    error( 'Not enough input arguments.');
end

if nargin < 4
    b_ind = [];
end

if ~any( ~isreal( x(:))) || ~any( ~isreal( y(:)))
    error( 'x and y must be analytic signals');
end

if ndims(x) ~= 3 || ndims(y) ~= 3
    error( 'x and y must be freq x time x trial.');
end

sz_x = size( x);
sz_y = size( y);
if ~isequal( sz_x( 2:end), sz_y( 2:end))
    error( 'x and y must be matched in time and trial length.');
end

n_pnt = sz_x(2);
n_trl = sz_x(3);
n_freq = max( [sz_x(1), sz_y(1)]);

% Complex conjugate of y
y = conj( y);

[boot_sz, n_boots] = size( b_ind);
if isempty( b_ind)
    boot_sz = n_trl;
    n_boots = 1;
    b_ind = (1:n_trl)';
end

perm_plv = nan( n_perm+1, n_freq, n_pnt);
permz = zeros( n_boots, n_freq, n_pnt);
real_plv = complex( zeros( n_freq, n_pnt));
if n_boots > 1
    real_absplv = zeros( n_freq, n_pnt);
else
    real_absplv = [];
end

perm_ind = nan( boot_sz, n_perm);
for pidx = 2 : n_perm+1
    perm_ind( :, pidx) = randperm( boot_sz, boot_sz);
end

for b_idx = 1 : n_boots
    ind = b_ind( :, b_idx);
    boot_plv = mean( exp( 1i*angle( x( :, :, ind) .* y(:, :, ind))), 3);
    real_plv = real_plv + boot_plv;
    if n_boots > 1
        real_absplv = real_absplv + abs( boot_plv);
    end

    if n_perm < 1
        continue;
    end

    perm_plv(1, :, :) = boot_plv;

    for pidx = 2 : n_perm+1
        if n_boots < 2
            fprintf( 'Permutation %d/%d    \r', pidx-1, n_perm)
        else
            fprintf( 'Bootstrapping: %d/%d; Permutation %d/%d     \r', b_idx, n_boots, pidx-1, n_perm)
        end
        % x_ind = ind( randperm( boot_sz, boot_sz));
        perm_plv( pidx, :, :) = mean( exp( 1i*angle( x(:,:, ind( perm_ind( :, pidx))) .* y(:,:, ind))), 3);
    end

    % squared root transformation for normal distribution
    perm_plv_sq = sqrt( abs( perm_plv));

    % z score
    [avg, sd] = normfit( perm_plv_sq( 2:end, :, :));
    permz( b_idx, :, :) = (perm_plv_sq( 1, :, :) - avg)./sd;
end
fprintf( '                                                   \n');

% average across bootstraps
real_plv = real_plv/n_boots;
real_absplv = real_absplv/n_boots;
z = permute( mean( permz, 1), [2, 3, 1]);

r = [];
r.absplv = [];
r.plv = real_plv;
r.n = boot_sz;
r.z = z;
if n_boots > 1
    r.absplv = real_absplv;
    r.plv = [];
end

end % function

