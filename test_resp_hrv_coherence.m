% This scripts contains examples for:
%#1. Calculate heart rate variability (HRV)
%#2. Calculate coherence between respiratory and HRV signals
%
%
% The example data is: example_ecg_resp.mat
% The mat file contains the following variables:
%    resp, respiratory signal (airflow, > 0 inhale)
%     ecg, ECG signal
%   srate, sampling rate (in Hz) of the ECG and respiratory signals
%     rpk, locations (in sample) of the R-peaks
% resp_ev, locations (in sample) of inhale and exhale onsets


%%## Retrieve heart rate variability (HRV) time series from R-peaks
clear; clc
load example_ecg_resp.mat
% execution time < 1 s
rrd = ECG_HRV( rpk, srate, ecg, 'show', 'yes');

% The resulting HRV time series can be accessed as rrd.hrvec


%%##2 Calculate coherence between respiratory and HRV signals
% execution time < 1 s
rrd = ECG_HRV( rpk, srate, ecg, 'rr_threshold', 2);

% data window with valid RR duration
samp_loc = rrd.rrloc(1) : rrd.rrloc(end);
hrv = rrd.hrvec( samp_loc);
resp_sub = resp( samp_loc);

% Remove DC component of the signals: 1) high-pass filtering; 2) demean
% high-pass filtering
%hrv = ft_preproc_highpassfilter( hrv, srate, 0.01, [], 'fir');
%resp_sub = ft_preproc_highpassfilter( resp_sub, srate, 0.01, [], 'fir');
% or use demean for speed consideration
hrv = hrv - mean( hrv);
resp_sub = resp_sub - mean( resp_sub);


% FFT length and window settings for power spectral density (PSD) and coherence analysis
param = [];
param.coh = struct( 'coh_threshold', 5, 'npermutations', 200);
param.coh.hrv_psd_nfft = 2^nextpow2( 3*60*srate);
param.coh.hrv_psd_win = hamming( param.coh.hrv_psd_nfft);
param.coh.coh_nfft = 2^nextpow2( 20*srate);
param.coh.coh_win = hamming( param.coh.coh_nfft);
param.npermutations = 200;

% Calculate PSD and coherence
% Requires signal toolbox for pwelch function
% The number of permutations determines the execution time.
% execution time is typically < 30 s for 200 perms
[hrv_psd, hrv_resp_coh] = ECG_HRV_Resp_Coh( srate, resp_sub, hrv, 'param', param.coh);


% To plot the coherence
figure; hold on;
plot( hrv_resp_coh.freq, hrv_resp_coh.coherence);
xlabel( 'Frequency (Hz)'); ylabel( 'Coherence');
set( gca, 'xlim', [0, 2]);
% significance thresold
plot( [0, 2], [1, 1] * hrv_resp_coh.threshold(1), 'color', 'r');


