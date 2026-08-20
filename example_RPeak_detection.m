% Raw ECG time sereis, 1xN vector
ecg = [];

% original sampling rate, in Hz
raw_srate = 2000;

% new sampling rate, in Hz
srate = 500;

% organize to fieldtrip data structure
N = length( ecg);

cfg = [];
cfg.trial{1} = ecg;
cfg.label{1} = 'ECG';
cfg.time{1} = linpace( 0, (N-1)/raw_srate, N);
cfg.fsample = raw_srate;
ecg_data = ft_datatype_raw( cfg);

% band-pass filtering
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfilttype = 'fir';
cfg.bpfreq = [0.5, 48];
ecg_data = ft_preprocessing( cfg, ecg_data);

% down-sampling
cfg = [];
cfg.resamplefs = srate;
cfg.detrend = 'no';
dn_ecg_data = ft_resampledata( cfg, ecg_data);

% preprocessed ECG
preproc_ecg = dn_ecg_data.trial{1};

% z score normalization
z = zscore( preproc_ecg);

% find R peaks
[~, pkloc] = findpeaks( z, 'MinPeakHeight', 2, 'MinPeakDistance', 0.4*srate);

% visualize results and add/remove detected events if necessary
EMarker( 'mat', preproc_ecg, 'srate', srate, 'events', pkloc, 'type', 'R Peak');
