function [pval, fval, df, es, misc] = G_RAnova1( mat, varargin)
% One-way repeated anova
%
% Input
%   x, sample x obervations (conditions) x ...
%
% Output
%   pval, between-subject p value
%   fval, between-subject f value
%   df, degree of freedom, data structure with the following fields
%       total
%       between
%       error
%       subject
%   es, effect size of a data structure with the following fields
%       etas, eta squared
%       partial_etas, partial eta squared
%       omegs, omeg squared
%       partial_omegs, partial omeg squared
%   misc, a data structure with the following fields
%       pval_within, within-subejct p value
%       fval_within, within-subject f value
%
% Steps to calculate one-way repeated ANOVA
% https://www.statology.org/repeated-measures-anova-by-hand/
% 1. Calculate total sum of squares (SST).
%   SST = s2_total( n_total - 1);
%   s2_total, the variance for the entire dataset
%   n_total, total number of obervations in the entire dataset
%
% 2. Calculate between sum of squares (SSB)
%    SSB = sum( n_j .* ( x_j - x_total)^2)
%    n_j, total number of obervations in the jth group
%    x_j, the mean of the jth group
%    x_total, the mean of the entire dataset
%
% 3. Calculate subject sum of squares (SSS)
%   SSS = sum( r2k/c) - (N2/rc);
%   r2k: squared sum of the kth subject
%   N: the grand total of the entire dataset
%   r: total number of subjects
%   c, total number of groups
%
% 4. Calculate SSE
%   SSE = SST - SSB - SSS
%
% df between = number of group - 1
% df subject = number of subjects - 1
% df error: df_between * df_subject
% MS between: SSB / df_between
% MS subject: SSS / df_subject
% MS error: SSE / df_error
% F: MS_between / MS_error
%
% Example
% mat = [30 28 16
%     14 18 10
%     24 20 18
%     38 34 20
%     26 28 14];
%
% mat = repmat( mat, [1, 1, 2]);
%
%
% % Testing (rm_anova1 from)
%  Trujillo-Ortiz, A., R. Hernandez-Walls and R.A. Trujillo-Perez. (2004).
%  RMAOV1:One-way repeated measures ANOVA. A MATLAB file. [WWW document].
%  URL http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=5576
% nb_samp = size( mat, 1);
% nb_cond = size( mat, 2);
% fact_cond = arrayfun( @(x) x*ones( nb_samp, 1), 1:nb_cond, 'un', 0);
% fact_cond = cell2mat( fact_cond(:));
% fact_subj = repmat( 1:nb_samp, [1, nb_cond])';
%
% X = [mat(:), fact_cond, fact_subj];
% [pp, ff] = rm_anova1( X)
%

if nargin < 2
    es_meas = 'all';
else
    es_meas = varargin{1};
end

all_es = {'etas', 'partial_etas', 'omegs', 'partial_omegs'};
if ~iscell( es_meas)
    if strcmpi( es_meas, 'all')
        es_meas = all_es;
    else
        es_meas = {es_meas};
    end
end

if ~all( ismember( es_meas, all_es))
    error( 'Unsupported effect size measurement.');
end

sz = size( mat);
nb_samp = sz( 1);
nb_cond = sz( 2);
perm_dim = [3 : numel( sz), 1, 2];

total_avg = mean( mean( mat, 2), 1);

% total sum of squares
sst = sum( sum( (mat - total_avg) .^ 2, 2), 1);
sst = permute( sst, perm_dim);

n_total = nb_samp*nb_cond;
df_total = n_total - 1;
df_subj = nb_samp - 1;
df_between = nb_cond - 1;
df_error = df_subj * df_between;

% mean of each condition/group
cond_avg = mean( mat, 1);

% between sum of squares
ssb = sum( sz(1) .* ( cond_avg - total_avg) .^2, 2);
ssb = permute( ssb, perm_dim);

% subject sum of squares (SSS)
sss = sum( sum( mat, 2).^2, 1)/nb_cond - (sum( sum( mat, 2), 1).^ 2)/(df_total + 1);
sss = permute( sss, perm_dim);
sse = sst - ssb - sss;

msa = ssb / df_between;
mss = sss / df_subj;
mse = sse / df_error;

% between subject f value
F_between = msa ./ mse;

% within-subject f
F_subj = mss ./ mse;

% Probability associated to the F-statistics.
p_between = 1 - fcdf( F_between, df_between, df_error);
p_subj = 1 - fcdf( F_subj, df_subj, df_error);

% partial eta squared
eta2 = ssb ./ (ssb + sse);
es = [];
if ismember( 'partial_etas', es_meas)
    es.partial_etas = eta2;
end

if ismember( 'etas', es_meas)
    es.etas = ssb ./ sst;
end

% pval.between = p_between;
% pval.subject = p_subj;
pval = p_between;

misc = [];
misc.pval_within = p_subj;

% fval.between = F_between;
% fval.subject = F_subj;
fval = F_between;
misc.fval_within = F_subj;

df.total = df_total;
df.between = df_between;
df.error = df_error;
df.subject = df_subj;

% % standard omeg squared, cannot be used for repeated meausres
% std_omegs = (ssb - df_between*mse) / (sst + mse);

if ~any( ismember( es_meas, {'omegs', 'partial_omegs'}))
    return;
end

% see mes1way.m. Harald Hentschke 2026
% github.com/hhentschke/measures-of-effect-size-toolbox
x = (df_between/n_total) * (msa - mse);
% omeg squared
if ismember( 'omegs', es_meas)
    es.omegs = x ./ (x + (1/nb_cond)*(mss - mse) + mse);
end

% partial omeg squared
if ismember( 'partial_omegs', es_meas)
    es.partial_omegs = x ./ (x + mse);
end

end % function
