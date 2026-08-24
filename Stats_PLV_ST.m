function [r, perm_plv] = Stats_PLV_ST( x, y, nperm, boot_ind)
% Phase locking value between x and y.
% Test the significance of the PLV using trial shuffling permutation method.
%
% Usage
%   [r, perm_plv] = Stats_PLV_ST( x, y, nperm, boot_ind);
% 
% Input
%   x, freq x time x trial analytic time series
%   y, freq x time x trial analytic time series
%   nperm, number of permutations
%   boot_ind, n_trl x n_boot bootstrap index
%
% Output
%   r, data structure with the following fields
%       plv, real plv (complex number)
%       z, Rayleigh's z of the real_plv
%       p, p value of the Rayleigh's z
% 
%   perm_plv, permuted plv (with the real plv being the first).
%       permutation (nperm + 1) x freq x time
% 
% GZ


% Complex conjugate
y = conj( y);

[boot_sz, n_boots] = size( boot_ind);
[~, n_pnt, n_trl] = size( x);
n_freq = max( [size( x, 1), size( y, 1)]);

perm_plv = nan( nperm+1, n_freq, n_pnt);
if n_boots < 1
    real_plv = mean( exp( 1i*angle( x .* y)), 3);
    for pidx = 2 : nperm+1
        fprintf( 'Permutation %d/%d         \r', pidx-1, nperm);
        ind = randperm( n_trl, n_trl);
        perm_plv( pidx, :, :) = mean( exp( 1i*angle( x .* y(:,:,ind))), 3);
    end

else
    real_plv = complex( zeros( n_freq, n_pnt));
    for b_idx = 1 : n_boots
        fprintf( 'Boostrapping %d/%d         \r', b_idx, n_boots);
        ind = boot_ind(:, b_idx);
        real_plv = real_plv + mean( exp( 1i*angle( x(:,:,ind) .* y(:,:,ind))), 3);
    end
    real_plv = real_plv/n_boots;

    for pidx = 2 : nperm+1
        fprintf( 'Permutation %d/%d         \r', pidx-1, nperm);
        sub_samp = sort( randsample( 1:n_trl, boot_sz, false));
        x_ind = randsample( sub_samp, boot_sz, false);
        y_ind = randsample( sub_samp, boot_sz, false);
        perm_plv( pidx, :, :) = mean( exp( 1i*angle( x(:,:,x_ind) .* y(:,:,y_ind))), 3);
    end
end

perm_plv(1, :, :) = real_plv;
perm_plv = abs( perm_plv);

% squared root transformation for normal distribution
perm_plv_sq = sqrt( perm_plv);

% z score
[avg, sd] = normfit( perm_plv_sq( 2:end, :, :));
perm_plvz = (perm_plv_sq - avg)./sd;

r = [];
r.plv = real_plv;
[r.p, r.z] = plv_rayleiz( real_plv, n_trl);
r.permz = permute( perm_plvz( 1, :, :), [2, 3, 1]);

end % function


function [p, z] = plv_rayleiz( plv, n)
% Rayleigh's z-statistic of phase locking value
%
% [p, z] = RayleighZ( plv, n)
%
% Input
%   plv, plv value
%   n, number of samples
%
% Output
%   p, p value
%   z, z score
%
% Script was adapted from CircStat2012 toolbox circ_rtest.
% https://github.com/poe-lab/CircStat2012a

r = abs( plv);
% compute Rayleigh's R (equ. 27.1)
R = n*r;
% compute Rayleigh's z (equ. 27.2)
z = R.^2 / n;
% compute p value using approxation in Zar, p. 617
p = exp(sqrt( 1 + 4*n + 4*(n^2 - R.^2)) - (1 + 2*n));
end
