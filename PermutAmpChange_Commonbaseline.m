function [stats, epoch_time, ev_epochs, comp_data, sd] = PermutAmpChange_Commonbaseline( ts, para, varargin)
% Within-conidtion event-induced amplitude change
%
% Usage
%   [stats, epoch_time, ev_epochs, amp] = PermutAmpChange_Commonbaseline( ts, para)
%
% Input
%   ts, time series vector
%       para = [];
%       para.freq = 2:2:200; % center frequencies for filtering
%       para.bandwidth = linspace( 2, 50, length( para.freq)); % bandwidth of each center frequency
%       para.srate = 500; % sampling rate, in Hz
%
%       para.events = your_events; % in samples
%       para.epoch = [-1, 2]; % in seconds, relative to events
%       para.baseline = [-0.5, 0]; % in seconds, relative to events
%       para.smooth_win = 0.01; % temporal smoothing, in seconds
%       para.nbpermuts = 10000; % number of permutations,
%       para.measure = 'amplitude'; % or 'power'
%       para.baseline_condition = 1; % which condition to use as baseline
%       para.trl = 'mean'; % 'mean' o r'median' value across trials
%
%
% Output
%   stats, a structure with the following fields
%       Z, z-score
%       %resp, baseline substracted spectrogram (amplitude)
%       %resp_percentchange, signal percent change
%   epoch_time, time vector for epochs
%   ev_epochs, raw amplitude epochs (large file size) freq x time x trial
%   comp_data, complex time series from Hilbert transform
%   sd, standard deviation of permutated distribution
%
% Statistical method reference
%   Ryan T. Canolty et al., Spatiotemporal dynamics of word preocessing in
%   the human brain. Front Neurosci. 2007 1(1) 185-196.
%

stats = [];
epoch_time = [];
ev_epochs = [];
comp_data = [];
sd = [];
amp = [];
nbpermuts = 10000;


if ~isempty( varargin)
    amp = abs( varargin{1});
    nbsamples = size( amp, 2);
    comp_data = varargin{1};
end

events = [];
if isfield( para, 'events')
    events = para.events;
end

epoch = [];
if isfield( para, 'epoch')
    epoch = para.epoch;
end

smooth_win = [];
if isfield( para, 'smooth_win')
    smooth_win = para.smooth_win;
end

if isfield( para, 'nbpermuts')
    nbpermuts = para.nbpermuts;
end

freq = para.freq;
nbfreqs = length( freq);
bandwidth = para.bandwidth;
srate = para.srate;

if ~isfield( para, 'measure')
    measure = 'amplitude';

else
    measure = para.measure;
    switch measure
        case {'amplitude', 'power'}
            % do nothing
        otherwise
            error( 'Unknown measurement.');
    end
end

if ~isfield( para, 'trl')
    para.trl = 'mean';
end

if ~ismember( para.trl, {'mean', 'median'})
    error( 'trl must be meand or median');
end

if isempty( amp)
    if ~isvector( ts)
        error( 'Time series must be a vector.');
    end

    nbsamples = length( ts);
    ts = reshape( ts, [1, nbsamples]);
    comp_data = zeros( nbfreqs, nbsamples);
    for freq_idx = 1 : nbfreqs
        bpfreq = freq( freq_idx) + [-1, 1]*bandwidth( freq_idx)/2;
        fprintf( 'Bandpass filtering %d/%d: %g - %g Hz               \r', freq_idx, nbfreqs, bpfreq(1), bpfreq(2));
        if bpfreq(1) < eps
            filt_ts = ft_preproc_lowpassfilter( ts, srate, bpfreq(2), [], 'fir');
        else
            filt_ts = ft_preproc_bandpassfilter( ts, srate, bpfreq, [], 'fir');
        end

        comp_data( freq_idx, :) = hilbert( filt_ts);
    end

    % amplitude
    amp = abs( comp_data);
    fprintf( 'Bandpass filtering                                       Done.\n');
end

% amplitude -> power
if strcmpi( measure, 'power')
    amp = amp .^2;
end

% temporal smoothing
if ~isempty( smooth_win)
    smooth_kernel = round( smooth_win * srate);
    if mod( smooth_kernel, 2) == 0
        smooth_kernel = smooth_kernel + 1;
    end
%    amp = G_Smooth( amp, 2, smooth_kernel);

    for k = 1 : size( amp, 1)
        amp( k, :) = conv( amp( k, :), ones( 1, smooth_kernel)/smooth_kernel, 'same');
    end

end

if isempty( events)
    return;
end

if isempty( epoch) || length( epoch) ~= 2
    error( 'Epoch must be given as [lo, hi] in seconds.');
end

if ~isfield( para, 'baseline_condition')
    baseline_condition = [];
else
    baseline_condition = para.baseline_condition;
    if ~isempty( baseline_condition)
        if ~ismember( baseline_condition, 1:length(events))
            error( 'Baseline condition is out of range.');
        end
    end
end

seg = floor( epoch * srate);
nbevs = length( events);
resp = cell( nbevs, 1);
resp_percentchange = cell( nbevs, 1);
Z = cell( nbevs, 1);
ev_epochs = cell( nbevs, 1);
for ev_idx = 1 : nbevs
    current_event = events{ ev_idx};
    if isempty( current_event)
        ev_epochs{ ev_idx} = [];
        resp{ ev_idx} = [];

    else
        % freq, samples x epochs
        ev_epochs{ ev_idx} = G_EpochTS( amp, current_event, seg);
        % average across trials
        switch para.trl
            case 'mean'
                resp{ ev_idx} = mean( ev_epochs{ev_idx}, 3);
            case 'median'
                resp{ ev_idx} = median( ev_epochs{ev_idx}, 3);
            otherwise
                error( 'error');
        end
    end
end % event

% baseline correction
epoch_time = epoch(1) : 1/srate : epoch(2);
baselineloc = [];
if isfield( para, 'baseline')
    baseline = para.baseline;
    baselineloc = G_RangeLoc( epoch_time, baseline);
end

if ~isempty( baseline_condition)
    fprintf( 'Using condition %d as a common baseline\n\n', baseline_condition);
    if ~isempty( events{ baseline_condition})
        if strcmpi(  para.trl, 'mean')
            avg = mean( resp{ baseline_condition}( :, baselineloc, :), 2);
        else
            avg = median( resp{ baseline_condition}( :, baselineloc, :), 2);
        end

        for ev_idx = 1 : nbevs
            resp{ ev_idx} = bsxfun( @minus, resp{ ev_idx}, avg);
            resp_percentchange{ ev_idx} = 100* bsxfun( @rdivide, resp{ ev_idx}, avg);
        end % event

    else
        error( 'Baseline condition has no events.');
    end

else
    fprintf( 'Perform baseline correction for each condition separately.\n\n');
    for ev_idx = 1 : nbevs
        if ~isempty( events{ ev_idx})
            if strcmpi( para.trl, 'mean')
                avg = mean( resp{ ev_idx}( :, baselineloc, :), 2);
            else
                avg = median( resp{ ev_idx}( :, baselineloc, :), 2);
            end

            resp{ ev_idx} = bsxfun( @minus, resp{ ev_idx}, avg);
            resp_percentchange{ ev_idx} = 100* bsxfun( @rdivide, resp{ ev_idx}, avg);
        else
            resp_percentchange{ ev_idx} = [];
        end

    end % event
end

sd = cell( nbevs, 1);
for ev_idx = 1 : nbevs
    if ~isempty( events{ ev_idx})

        current_event = events{ ev_idx};
        permut_avg = zeros( nbfreqs, nbpermuts);
        for pidx = 1 : nbpermuts
            fprintf( 'Permutation %d/%d               \r', pidx, nbpermuts);
            shift_amount = randi( [int32( 0.1*srate), nbsamples], 1);
            permut_onsets = mod( current_event + shift_amount, nbsamples);
            % make sure there' enough sample at the beginning and end
            while any( permut_onsets(:) < -1*epoch(1)*srate + 1) || any( permut_onsets(:) > nbsamples - epoch(2)*srate - 1)
                shift_amount = randi( [int32( 0.1*srate), nbsamples], 1);
                permut_onsets = mod( current_event + shift_amount, nbsamples);
            end
            switch para.trl
                case 'mean'
                    permut_avg( :, pidx) = mean( amp( :, permut_onsets), 2);
                otherwise
                    permut_avg( :, pidx) = median( amp( :, permut_onsets), 2);
            end
        end

        sd{ ev_idx} = std( permut_avg, [], 2);
        Z{ ev_idx} = bsxfun( @rdivide, resp{ ev_idx}, sd{ ev_idx});

    else
        Z{ ev_idx} = [];
    end
end % event
fprintf( 'Permutation                           Done.\n');

stats.resp_percentchange = resp_percentchange;
stats.Z = Z;
stats.resp = resp;

end % function
