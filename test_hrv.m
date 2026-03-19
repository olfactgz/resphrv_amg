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


