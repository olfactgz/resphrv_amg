function [pxx, hrv_resp_coh, perm_coh] = ECG_HRV_Resp_Coh( srate, resp, hrv, varargin)
% Coherence between heart rate variability and respiratory signal.
%
%
% Input
%   srate, sampling rate in Hz
%   resp, respiratory time series
%   hrv, heart rate variability time series
%
%   Key-Value pairs
%       'hrv_plomb_maxfreq', 10,... Calculate PSD using plomb method
%       'rr', Raw RR duration (without interpolation, not equally distributed,
%                           outliers could have been interpolated or fixed)
%       'rrt', time point of rr
%       'hrv_psd_nfft', [],... Calculate PSD using FFT method
%       'hrv_psd_overlap', 0.95,... % [0-1] hrv_psd_nfft * hrv_psd_overlap
%       'hrv_psd_win', [],...
%       'coh_nfft', [],... Coherence analysis
%       'coh_overlap', 0.95, ...
%       'coh_win', [], ...
%       'coh_threshold', [], ...
%       'npermutations', 0,... number of permutations of coherence significance test
%       'prctile', [95, 99] significance threshold
%
%
% Steps for coherence analysis
% #1. Identify R peaks on ECG signal.
% #2. Calculate heart-rate variability based on R peaks (See ECG_HRV.m)
% #3. Epoching ECG and Respiratory signals.
% #4. Reject bad epochs.
% #5. Calculate coherence.
%
%
%rrd = ECG_HRV( data.Rpeak, srate, data.ecg, 'param', param.hrv);
%
%% data window with valid RR duration
%t = rrd.hrvt;
%samp_loc = G_RangeLoc( t, rrd.rrt( [1, end]));
%hrv = rrd.hrvec( samp_loc);
%resp = data.resp( samp_loc);
%
%%% high-pass filtering
%%hrv = ft_preproc_highpassfilter( hrv, srate, 0.01, [], 'fir');
%%resp = ft_preproc_highpassfilter( resp, srate, 0.01, [], 'fir');
%% for speed consideration
%hrv = hrv - mean( hrv);
%resp = resp - mean( resp);
%
% GZ


opt = struct( ...
    'hrv_plomb_maxfreq', 10,...
    'rr', [],...
    'rrt', [],... % time vector for rr_duration
    'hrv_psd_nfft', [],...
    'hrv_psd_overlap', 0.9,... % [0-1] hrv_psd_nfft * hrv_psd_overlap
    'hrv_psd_win', [],...
    'hrv_zmethod', 'overallnormfit',... 'overallnormfit' | 'overallz' | 'epochz' | 'epochnormfit'
    'resp_zmethod', 'overallnormfit',... 'overallnormfit' | 'overallz' | 'epochz' | 'epochnormfit'
    'coh_nfft', [],...
    'coh_overlap', 0.9, ...
    'coh_win', [], ...
    'coh_threshold', [], ...
    'perm_method', '',... to implement
    'npermutations', 0, ...% number of permutations
    'prctile', [95, 99],...
    'param', []);

opt = G_SparseArgs( opt, varargin);


% power spectral density of heart rate variability
% f_max = 0.5 / ( (locs(end) - locs(1)) / (length( locs)-1))
pxx1 = [];
f1 = [];
pxx2 = [];
f2 = [];
if exist( 'plomb', 'file')
    rr_dur = opt.rr;
    plm_freq = opt.hrv_plomb_maxfreq;

    if ~isempty( rr_dur)
        rr_dur_t = opt.rrt;
        if isempty( rr_dur_t)
            error( 'rrt must be set to calculate the psd of rr using plomb method.');
        end

        [pxx1, f1] = plomb( rr_dur, rr_dur_t, plm_freq);
    end

    N = length( hrv);
    hrv_t = linspace( 0, (N-1)/srate, N);
    [pxx2, f2] = plomb( hrv, hrv_t, plm_freq);
end


psd_nfft = opt.hrv_psd_nfft;
if isempty( psd_nfft)
    pxx3 = [];
    f3 = [];
    pxx4 = [];
    f4 = [];
else
    psd_overlap = round( opt.hrv_psd_overlap * psd_nfft);
    if exist('OCTAVE_VERSION', 'builtin') ~= 0
        [pxx3, f3] = pwelch( hrv, opt.hrv_psd_win, opt.hrv_psd_overlap, psd_nfft, srate);
        [pxx4, f4] = pwelch( resp, opt.hrv_psd_win, opt.hrv_psd_overlap, psd_nfft, srate);
    else
        % matlab
        [pxx3, f3] = pwelch( hrv, opt.hrv_psd_win, psd_overlap, psd_nfft, srate);
        [pxx4, f4] = pwelch( resp, opt.hrv_psd_win, psd_overlap, psd_nfft, srate);
    end
end


pxx = [];
pxx.hrv.raw_rrd.psd = pxx1;
pxx.hrv.raw_rrd.freq = f1;
pxx.hrv.rrd_interp.plomb.psd = pxx2;
pxx.hrv.rrd_interp.plomb.freq = f2;
pxx.hrv.rrd_interp.fft.psd = pxx3;
pxx.hrv.rrd_interp.fft.freq = f3;

pxx.resp.psd = pxx4;
pxx.resp.freq = f4;

% % Calculate coherence
coh_nfft = opt.coh_nfft;
if isempty( coh_nfft)
    coh = [];
    freq = [];
    coh_thr = [];
else
    overlap_len = round( coh_nfft * opt.coh_overlap);
    hrv_mat = squeeze( G_EpochTS( hrv, [], coh_nfft, overlap_len));
    resp_mat = squeeze( G_EpochTS( resp, [], coh_nfft, overlap_len));

    % remove bad epochs
    coh_zthr = opt.coh_threshold;
    if ~isempty( coh_zthr)
        hrv_mat_z = MyZscore( hrv, opt.hrv_zmethod, coh_nfft, overlap_len);
        resp_mat_z = MyZscore( resp, opt.resp_zmethod, coh_nfft, overlap_len);
        hrv_mat_z = max( abs( hrv_mat_z), [], 1);
        resp_mat_z = max( abs( resp_mat_z), [], 1);
        epoch_z = max( [hrv_mat_z; resp_mat_z], [], 1);
        seg_ind = epoch_z < coh_zthr;
    else
        seg_ind = true( size( hrv_mat, 2), 1);
    end

    if sum( seg_ind) < 10
        error( 'Not enough epochs to calculate coherence.');
    else
        fprintf( '%d (total %d) epochs for coherence calculation\n', sum( seg_ind), numel( seg_ind));
    end

    [coh, freq] = G_Coh( hrv_mat( :, seg_ind), resp_mat( :, seg_ind), ...
        'srate',srate, 'nfft', coh_nfft, 'window', opt.coh_win);

    % significance of coherence
    perm_coh = [];
    nb_perm = opt.npermutations;
    if isempty( nb_perm) || nb_perm < 1
        coh_thr = [];
    else
        % circular shifting, bad epochs?
        x = hrv_mat( :, seg_ind);
        y = resp_mat( :, seg_ind);
        nb_trl = size( x, 2);

        perm_coh = nan( nb_perm, numel( coh));
        for k = 1 : nb_perm
            fprintf( 'Permuation %d/%d\n', k, nb_perm);
            perm_coh( k, :) = G_Coh( x, y( :, randperm( nb_trl, nb_trl)), ...
                'srate',srate, 'window', opt.coh_win);
            % y_epoch = squeeze( G_EpochTS( circshift( s_resp, randperm( npnt, 1)), [], nfft, overlap_len));
            % perm_coh( k, :) = G_Coh( s_hrv_epoch, y_epoch, ...
            % 'fs',srate, 'hwinlen', hanning( nfft));
        end

        coh_thr = prctile( max( perm_coh, [], 2), opt.prctile);
    end
end

hrv_resp_coh = [];
hrv_resp_coh.coherence = coh;
hrv_resp_coh.freq = freq;
hrv_resp_coh.prctile = opt.prctile;
hrv_resp_coh.threshold = coh_thr;
end % function


function x_mat_z = MyZscore( x, method, seg_len, overlap)
switch lower( method)
    case 'overallnormfit'
        [avg, sd] = normfit( x);
        xz = (x - avg) / sd;
        x_mat_z = G_EpochTS( xz, [], seg_len, overlap);
        x_mat_z = permute( x_mat_z, [2, 3, 1]);

    case 'overallzscore'
        xz = zscore( x);
        x_mat_z = G_EpochTS( xz, [], seg_len, overlap);

    case 'epochnormfit'
        x_mat = G_EpochTS( x, [], seg_len, overlap);
        x_mat = permute( x_mat, [2, 3, 1]);
        x_mat_z = nan( size( x_mat));
        for k = 1 : size( x_mat, 2)
            [avg, sd] = normfit( x_mat( :, k));
            x_mat_z( :, k) = (x_mat(:,k) - avg)/sd;
        end

    case 'epochzscore'
        x_mat = G_EpochTS( x, [], seg_len, overlap);
        x_mat_z = zscore( x_mat, [], 2);
        x_mat_z = permute( x_mat_z, [2, 3, 1]);

    otherwise
        error( 'Unknown z score method.');
end
end
