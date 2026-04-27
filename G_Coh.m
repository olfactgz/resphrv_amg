function [coh, freq, img_coh, raw_psd, fft_result] = G_Coh( x, y, varargin)
% calculate magnitude squared coherence between x and y
%
% Output
%   coh - magnetidue squared coherence
%   freq - frequencies
%   Pxx/Pyy - power spectral density of x/y
%   nbsegs - number of segments
%
%
% Example Usage:
%   [coh, freq] = G_Coh( x, y, 'srate', 1000, 'nfft', 2048, 'window', hamming(2048));
%
% x - 1st signal
% y - 2nd signal
%     samples x segments
%
% options with default values
% srate - sampling frequency in Hz, default 2*pi
% nfft - nfft, default seg_len
% window - windowing, default ones( 1, nfft)
% plot - 'yes'|'no'
%
%
% Haishen

warning( 'See G_Coherence for newer version.');

opt = struct( 'srate', [],...
    'nfft', [],...
    'window', [],...
    'plot', 'no',...
    'foi', []);

opt = G_SparseArgs( opt, varargin);


if ~ismatrix( x) || ~ismatrix( y) || ~isequal( size(x), size(y))
    error( 'x and y must be matrix with the same dimension.');
end

[seg_len, nb_segs] = size( x);

srate = opt.srate; % Hz
if isempty( srate)
    srate = 2 * pi;
end

nfft = opt.nfft; % in samples
if isempty( nfft)
    nfft = seg_len;
end

if round(nfft/2) ~= nfft/2
    error( 'nfft has to be an even number.');
end

if seg_len > nfft
    % nfft = 2 .^ nextpow2( seg_len);
    error( 'nfft is shorter than segment length.');
end

% window/taper, coulde be a scalar or a vector
win = opt.window;
if isempty( win)
    win = ones( seg_len, 1);
end

if isscalar( win)
    hwin = ones( win, 1);
elseif isvector( win)
    hwin = win;
else
    error( 'Unknown window type.');
end

win_len = length( hwin);

% nfft is always made to be even
freq = linspace( 0, srate/2, nfft/2 + 1);

% normalization factor of power spectrum
NG = mean( hwin .^ 2);

% window segmentations
hwin = repmat( hwin(:), [ 1, nb_segs]);

% the window length should be the same as the segmentation length, if not padding with zeros
if win_len > seg_len
    windiff = win_len - seg_len;
    pre_pad = floor( windiff/2);
    post_pad = windiff - pre_pad;
    x = padarray( x, [pre_pad, 0], 'pre');
    x = padarray( x, [post_pad, 0], 'post');
    y = padarray( y, [pre_pad, 0], 'pre');
    y = padarray( y, [post_pad, 0], 'post');
end

x = hwin .* x;
y = hwin .* y;

% spectrum from fft
Fx = fft( x, nfft, 1);
Fy = fft( y, nfft, 1);

nfft = double( nfft);
Fx = Fx( 1 : ( nfft/2 + 1), :);
Fy = Fy( 1 : ( nfft/2 + 1), :);

% power spectral density: (1/(FsxN))*F(f)*conj(F(f))
norm_factor = (1 /( srate*nfft*NG));
Pxx = norm_factor * ( Fx .* conj( Fx));
Pxx( 2:end-1, :) = Pxx( 2:end-1, :) .* 2;
Pyy = norm_factor * ( Fy .* conj( Fy));
Pyy( 2:end-1, :) = Pyy( 2:end-1, :) .* 2;

% cross-spectral density
Pxy = norm_factor * ( Fx .* conj( Fy));
Pxy( 2:end-1, :) = Pxy( 2:end-1, :) .* 2;

% magnitude squared coherence estimation
coh = ((abs( mean( Pxy, 2))) .^2 ) ./ (mean( Pxx, 2) .* mean( Pyy, 2));
img_coh = (mean( Pxy, 2)) ./ sqrt( (mean( Pxx, 2) .* mean( Pyy, 2)));

raw_psd = [];
raw_psd.x = Pxx;
raw_psd.y = Pyy;
raw_psd.xy = Pxy;
raw_psd.norm_factor = norm_factor;

fft_result = [];
fft_result.x = Fx;
fft_result.y = Fy;

% plot result
if ~strcmpi( opt.plot, 'yes')
    return;
end

if ~isempty( opt.foi)
    f_oi = opt.foi;
else
    f_oi = freq( [1, end]);
end

if abs( srate - 2*pi) < 1e-10
    xlbl = 'Frequency (radian)';
else
    xlbl = 'Frequency (Hz)';
end

ax_param = {'xlim', f_oi, 'box', 'off', 'tickdir', 'out'};
ln_prop = {'linewidth', 1, 'color', 'k'};

figure;
subplot( 311);
plot( freq, coh, ln_prop{:});
title( 'Magnitude squared coherence');
ylabel( 'Coherence'); xlabel( xlbl);
set( gca, ax_param{:});

% psd
subplot( 312);
plot( freq, 10*log10( mean( Pxx, 2)), ln_prop{:});
title( 'Power Spectral Density of signal x');
ylabel( 'PSD (db)'); xlabel( xlbl);
set( gca, ax_param{:});

subplot( 313);
plot( freq, 10*log10( mean( Pyy, 2)), ln_prop{:});
title( 'Power Spectral Density of signal y');
ylabel( 'PSD (db)'); xlabel( xlbl);
set( gca, ax_param{:});


% randn('state',0);
% h = fir1(30,0.2,rectwin(31));
% h1 = ones(1,10)/sqrt(10);
% r = randn( 16384,1); % 16384
% x = filter(h1,1,r);
% y = filter(h,1,x);
% % mscohere(x,y,ones( 1, 1024),512,1024)
% xmat = squeeze( G_EpochTS( x, [], 1024, 512));
% ymat = squeeze( G_EpochTS( y, [], 1024, 512));
%
% [c1, f1] = G_Coh( xmat, ymat, 'nfft', 1024);


% end
%
% % significance of coherence (M. Bourguignon et al. / NeuroImage 55 (2011) 1475�C1479)
% Nf = 185;
% Nc = 306;
% P = 0.01;
% L = 300; % number of dijoint epochs
% N = Nf * Nc; % number of requencies x number of channels
% p = 1 - (1-P) ^ (1/N); % P-global false positive rate for all frequencies and all channels .05;
%  % p - significance value for individual channels p, L-the number of disjoint epochs used
% Ct = 1 - p ^ (1/(L-1))


