function [opt, key] = G_SparseArgs( opt, varargin)
% Function sparse key-value pairs input.
% 
% Input
%   opt: structure with keys and default values
% 	args: user input key-value pairs. All keys are case-insensitive.
%       'param', reseved fieldname (See Example Usage below)
% 
% Output
% 	opt: a structure of which each field has a user-specified value
%       The fieldnames are in lower case.
%   key, all fields that were associated with a new values
%       
% Example Usage
%     function test_G_SparseArgs( varargin)
%     opt = struct( 'Width', [],... % case insensitive
%         'height', [],...
%         'param', []);
%     
%     [opt, key] = G_SparseArgs( opt, varargin);
%     
%     % Note. opt.width instead of opt.Width. All fieldnames are in lower case.
%     fprintf( 'Width: %g; Height: %g\n', opt.width, opt.height);
%     fprintf( 'keys: %s\n', strjoin( key, ', '));
%     
%     end % function
% 
%   The input arguments can be passed to the function in the following ways
%   1). test_G_SparseArgs( 'width', 10, 'hegith', 5);
% 
%   2). param = [];
%       param.width = 10;
%       param.height = 5;
%       test_G_SparseArgs( 'param', param);
% 
%   3). Overwrite those parameters in the 'param' structure
%       param = [];
%       param.width = 10;
%       param.height = 5;
%       test_G_SparseArgs( mat, 'param', param, 'height', 20); 
%       % height will be 20
% 
% GZ

if nargin ~= 2
    error( 'Syntax: G_SparseArgs( opt, varargin);');
end

if ~isstruct( opt)
    error( 'OPT must be a structure.');
end

% sparse input arguments
opt = LowerFname( opt);
opt_nam = fieldnames( opt);

args = varargin{1};
if mod( length( args), 2) ~= 0
   error('KEY-VALUE pairs needed.')
end

% 2 x N cell array
args = reshape( args, 2, []);

% keys
key = args( 1, :);

% key must be char
key_ok = cellfun( @ischar, key);
if ~all( key_ok)
    error( 'KEY must be of type char.');
end

% to lower case
key = cellfun( @lower, key, 'UniformOutput', false);

% Are keys all found int OPT
ia = ismember( key, opt_nam);
if any( ~ia)
    error( 'Unrecognized parameter (s): %s', strjoin( key( ~ia), ', '));
end

for k = 1 : length( key)
    opt.( key{ k}) = args{ 2, k};
end

if ~isfield( opt, 'param') || ~isstruct( opt.param)
    return;
end

opt.param = LowerFname( opt.param);
fn = fieldnames( opt.param);
ia = ismember( fn, opt_nam);
if any( ~ia)
    error( 'Unrecognized parameter (s) in PARAM: %s', strjoin( fn( ~ia), ', '));
end

% key exists in key-value already
ia = ismember( fn, key);
fn = fn( ~ia);
for k = 1 : length( fn)
    n = fn{ k}; 
    opt.(n) = opt.param.( n);
end

key = [key(:); fn(:)];
key = unique( key, 'stable');

end % function


%% sub-function
function new_s = LowerFname( s)
new_s = [];
if isempty( s)
    return;
end

fn = fieldnames( s);
for k = 1 : length( fn)
    new_s.( lower( fn{ k})) = s.( fn{ k});
end
end % sub-func
