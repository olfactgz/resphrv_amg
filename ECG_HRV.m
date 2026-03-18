function [r, ax] = ECG_HRV( rpk, srate, varargin)
% Retrieve continuous R-R duration time series from R peaks.
%
%
% Usage
%   r = ECG_HRV( rpk, srate);
%   r = ECG_HRV( rpk, srate, ecg);
%   r = ECG_HRV( rpk, srate, ecg, KEY, VALUE,...);
%   r = ECG_HRV( rpk, srate, KEY, VALUE,...);
%
%
% Input
%    rpk,  R peak location in samples
%   srate, sampling rate in Hz
%     ecg, ECG time series or scalar of the total number of time points
%          If ecg is not set, the resulting HRV time series will be shorter
%          than the original ECG signal, which equals to the last R peak location
%
%   Key-Value pairs
%       'interpolation_method', interpolation method (default to spline)
%       'location', time point which the R-R interval will be assigned at
%                   'left', 1st R peak (default)
%                   'right', 2nd R peak
%                   'middle', center of between two R peaks
%       'rr_fix_method', not used yet 'none' (default), 'interpolation' (preferred),
%                        or 'remove' (not yet implemented)
%       'rr_interpolation_method', interpolation method for RR intervals
%       'rr_threshold', scalar (z score threshold) or vector;
%                   Replace outliers of RR interval using interpolation if set
%       'show', 'yes' | 'no'
%       'toi', time of interesting for plotting ([lo, hi] in seconds)
%       'keepshape', 'no' (default) | 'yes'
%
%
% Output
%   r, data structure with the following fields
%      'hrvec', interpolated rrd time series
%       'hrvt', time vector of interpolated rrd time series
%        'rrt', time point where the rrd value was assigned at
%      'rrloc', sample location where the rrd value was assigned at
%     'rrorig', raw RR interval
%   'rrinterp', same to rrdorig if orig_rr_dur had no outliers
%
%
%
% This program is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation (version 3).
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with
% this program. If not, see <https://www.gnu.org/licenses/>.
%
%
%
% G. Zhou
%
% 05/13/2025 Wrote this script.
%


opt = struct( 'interpolation_method', 'spline', ... see interp1
    'location', 'left', ... 'left', 'right', 'middle'
    'rr_fix_method', 'none',...'none', 'remove', 'interpolation'
    'rr_threshold', [],... % z score threshold or vector
    'rr_interpolation_method', 'spline',...
    'show', 'no',...
    'toi', [],...
    'keepshape', 'no',... 'yes', or 'no';
    'param', []);

if nargin < 2
    error( 'R peak location (in samples) and sampling rate must be set.');
end

if ~isvector( rpk) || numel( rpk) < 2
    error( 'Not enough R peaks.');
end

if ~isempty( varargin) && ~ischar( varargin{1})
    ecg = varargin{1};
    varargin = varargin( 2:end);
else
    ecg = [];
end

opt = G_SparseArgs( opt, varargin);

% Total number of data time points
if isempty( ecg)
    N = rpk(end);
else
    if isscalar( ecg)
        N = ecg;
    else
        N = length( ecg);
    end
end

if N < rpk( end)
    error( 'Something went wrong.');
end

% Midlle point between R Peaks
mid_rpk_samp = (rpk( 1 : end-1) + rpk( 2 : end)) / 2;
mid_rpk_samp = round( mid_rpk_samp);

t = linspace( 0, (N-1)/srate, N);
rpk_t = t( rpk);
mid_rk_t = t( mid_rpk_samp);

% r-r interval
rr_dur = diff( rpk_t);
orig_rrd = rr_dur;

% rrd_z = (rr_dur - mean(rr_dur)) / std( rr_dur);
[orig_avg, orig_sd] = normfit( rr_dur);
rrd_z = (rr_dur - orig_avg) / orig_sd;

% Replace outliers using interpolation method
% This seems to be an important step before HRV-Respiration coherence analysis
rr_thr = opt.rr_threshold;
if isempty( rr_thr)
    good_ind = true( size( rr_dur));

else
    if isscalar( rr_thr)
        good_ind = abs( rrd_z) < rr_thr;

    else
        if numel( rr_thr) ~= numel( rr_dur)
            error( ['RR duration threshold must either be a scalar of zscore ' ...
                'or a logical vector with the same length to RR duration.']);
        end
        good_ind = logical( rr_thr);
    end

    if all( good_ind)
        fprintf( 'RR duration has no outliers exceeding the threshold.\n');
    end
end

if any( ~good_ind)
    warning( 'Original RR duration was interpolated to fix outliers.');
    tmp = 1 : numel( rr_dur);
    rr_dur = interp1( tmp( good_ind), rr_dur( good_ind), tmp, opt.rr_interpolation_method);
end

switch lower( opt.location)
    case 'middle'
        % assign RR duration to the middle point between two R peaks
        rpk_t4interp = mid_rk_t;
        rr_loc = mid_rpk_samp;
    case 'left'
        % assign RR duration to the 1st R peak
        rpk_t4interp = rpk_t( 1 : end-1);
        rr_loc = rpk( 1 : end-1);
    case 'right'
        % assign RR duration to the 2nd R peak
        rpk_t4interp = rpk_t( 2:end);
        rr_loc = rpk( 2:end);
    otherwise
        error( 'location might be one of the folowing: left, right or middle');
end

switch lower( opt.keepshape)
    case 'yes'
        % assign the median to the first and last data points
        med_dur = median(rr_dur);

        if rr_loc(1) > 1
            rpk_t4interp = [t(1); rpk_t4interp(:)];
            rr_loc = [1; rr_loc(:)];
            rr_dur = [med_dur; rr_dur(:)];
        end

        if rr_loc(end) < N
            rpk_t4interp = [rpk_t4interp(:); t(end)];
            rr_loc = [rr_loc(:); N];
            rr_dur = [rr_dur(:); med_dur];
        end

    case 'no'
        % do nothing
    otherwise
        error( 'unknown value for keepshape');
end

rrd_t = rpk_t4interp;
rrd_interp = interp1( rpk_t4interp, rr_dur, t, opt.interpolation_method);
rrd_interp( isnan( rrd_interp)) = 0;
rrd_interp( t < rpk_t4interp(1)) = 0;
rrd_interp( t > rpk_t4interp(end)) = 0;

r = struct( 'rrt', rrd_t, ... time point where the rrd value was assigned at
    'rrloc', rr_loc, ...sample location where the rrd value was assigned at
    'rrorig', orig_rrd, ...raw RR duration
    'rrinterp', rr_dur, ... same to rrdraw if orig_rr_dur had no outliers
    'hrvec', rrd_interp,...interpolated rrd time series
    'hrvt', t, ...% time vector of interpolated rrd time series
    'srate', srate);

fprintf( ['Original RR duration range: %.3f ([%.3f, %.3f]) seconds\n' ...
    'Interpolated RR duration range: %.3f ([%.3f, %.3f]) seconds\n'], ...
    range( r.rrorig), min( r.rrorig), max( r.rrorig), ...
    range( r.rrinterp), min( r.rrinterp), max( r.rrinterp));


%% To implement: Calculate HRV measurements



%% illustration
ax = [];
if ~strcmpi( opt.show, 'yes')
    return;
end

if any( ~good_ind)
    figure( 'name', 'RR duration outlier interpolation');
    ax = subplot( 1, 6, 1);
    h = histfit( orig_rrd);
    set( h(1), 'facecolor', 0.1*ones( 1, 3), 'EdgeColor', 0.6*ones( 1, 3), 'linewidth', 0.25);
    set( h(2), 'color', 'b');
    % % Matlab syntax
    %    h(1).FaceColor = 0.1*ones( 1, 3);
    %    h(1).EdgeColor = 0.6*ones( 1, 3);
    %    h(1).LineWidth = 0.25;
    %    h(2).Color = 'b';
    ylm = max( get(h(1), 'YData'))*1.05;
    if isscalar( rr_thr)
        line( orig_avg + [-1, 1]*rr_thr*orig_sd, ylm( [1, 1]), ...
            'LineWidth', 2, 'Color', 'r');
    end
    lh = legend( {'', 'fit line', sprintf('|z| < %g', rr_thr)}, 'box', 'off', 'Location', 'northwest');
    %    lh.IconColumnWidth = lh.IconColumnWidth/4;

    ylabel( 'Count');
    xlabel( 'Time (s)');
    set( gca, 'box', 'off', 'tickdir', 'out');

    subplot( 1, 6, 2:6);  hold on;
    h1 = plot( rpk_t( 1:end-1), orig_rrd, 'k-o');
    b = G_Boundary1d( ~good_ind);
    for k = 1 : size( b, 1)
        tmp = b( k, 1) : b( k, 2);
        h2 = plot( rpk_t( tmp), orig_rrd( tmp), 'r-o', 'linewidth', 1, 'markerfacecolor', 'r');
        h3 = plot( rpk_t( tmp), rr_dur( tmp), '-*', 'linewidth', 1, 'color', [0.2 0.8 0.1]);
    end

    legend( [h1, h2, h3], {'original', 'before fixiation', 'fixed'});
    if ~isempty( ecg) && ~isscalar( ecg)
        tmp_ecg = G_Rescale( ecg, [-1, 0]);
        h4 = plot( t, tmp_ecg, 'color',[122 4 156]/255, 'linewidth', 1);
        legend( [h1, h2, h3, h4], {'original', 'before fixiation', 'fixed', 'ecg'}, ...
            'box', 'off', 'Location', 'west');
    end
    set( gca, 'xlim', t( [1, end]));
end

% data range with valid RR intervals
if isempty( opt.toi)
    % len = 0.05*(t(end) - t(1));
    % t_win = [-1, 1] * len + mean(t);
    t_win = rrd_t( [1, end]);
else
    t_win = opt.toi;
end


figure( 'Name', 'Example ECG_HRV'); hold on;
h1 = plot( t, rrd_interp, 'k-', 'LineWidth', 1);
h2 = plot( rpk_t4interp, rr_dur, 'ro', 'markerfacecolor', 'r', 'LineWidth', 1);
ylabel( 'RR interval (s)');
xlabel( 'Time (s)');

if isempty( ecg) || isscalar( ecg)
    h3 = [];
    warning( 'set ecg to plot');
else
    toi_loc = G_RangeLoc( t, t_win);
    bd = [min( ecg( toi_loc)), max( ecg( toi_loc))];
    ref_range = rrd_interp( toi_loc);
    rescale_ecg = G_Rescale( ecg, ref_range, 'bounds', bd);
    re_pk = G_Rescale( ecg( rpk), ref_range, 'bounds', bd);
    plot( rpk_t4interp, re_pk(1:end-1), 'ro', 'markerfacecolor', 'r', 'LineWidth', 1);
    sb = G_Rescale( [0, 1], ref_range, 'bounds', bd);
    h3 = plot( t, rescale_ecg, 'linewidth', 1, 'color', [146, 8, 165]/255);
    line( t_win( [end, end]), sb, 'LineWidth', 6, 'color', [146, 8, 165]/255);

    for k = 1 : numel( rr_dur)
        if prod( t( rpk( k)) - t_win) < 0
            continue;
        end
        text( t( rpk( k)), rescale_ecg( rpk( k)), num2str( rr_dur(k)));
    end

    % yyaxis right
    % plot( t, ecg, 'linewidth', 1, 'color', [146, 8, 165]/255);
    % for k = 1 : numel( rr_dur)
    %     if prod( t( rpk( k)) - t_win) > 0
    %         continue;
    %     end
    %     text( t( rpk( k)), ecg( rpk( k)), num2str( rr_dur(k)));
    % end
    % ylabel( 'ECG amplitude (mv)');
end

set( gca, 'xlim', t_win, 'box', 'off', 'TickDir', 'out');

if isempty( h3)
    legend( [h1, h2], {'Interpolated RR', 'Data for interpolation'}, 'box', 'off');
else
    legend( [h1, h2, h3], {'Interpolated RR', 'Data for interpolation', 'ECG'}, 'box', 'off');
end

end % function

