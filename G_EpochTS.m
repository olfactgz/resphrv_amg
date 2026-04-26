function mat_epochs = G_EpochTS( mat, ev, seg, varargin)
% Epoch time series into segments
%
% Input
%   mat, channel x samples
%   ev, vector of event locations (in samples)
%   seg, [lo_boundary, hi_boundary] relative to each event (in
%         samples) where a negative number means before the event
%
%   To epoch data into continuous segments, set ev to [] and seg to
%   segment length, the overlap between segments can be specified by a
%   fourth input argument.
%
% Output
%   mat_epochs, channel x samples x epoch
%
% Example
%   % epoch continous data into segments centered at ev
%   mat_epochs = G_EpochTS( rand( 1, 1000), [ 400, 600], [-399, 300]);
%   % epoch continous data into segments with no overlaps
%   mat_epochs = G_EpochTS( rand( 1, 1000), [], 100);
%   % epoch continous data into segments with overlap
%   mat_epochs = G_EpochTS( rand( 1, 1000), [], 100, 10);
%
% 01-24-2023, use for loop to replace arrayfun to make it faster
%
% 20-Feb-2018, fix bugs in continuous segmentation
%
% GZ


if nargin < 3
    error( 'Not enought input arguments');
end

% TODO validate input arguments
% seg can be cell array
if ~ismatrix( ev) || ~ismatrix( mat)
    error( 'error');
end

% overlap between continuous segments
if nargin > 3
    overlap = varargin{ 1};
    if ~ismatrix( overlap)
        error( 'The overlap must be a number.');
    end
else
    overlap = 0;
end

% make sure that input is channel x samples
if isvector( mat)
    mat = reshape( mat, 1, []);

else
    if ~ismatrix( mat)
        error( 'Input data must be channel x sample');
    end
end

% if ~issorted( ev)
%     warning( 'Events are not sorted.');
% end

[nb_chnls, nb_pnts] = size( mat);
% it's not likely that number of channels is greater than that of samples
if nb_chnls > nb_pnts
    warning( 'The number of channels is greater than signal length.');
end

if isempty( ev)
    % continuous segmentation
    if ~ismatrix( seg)
        error( 'For empty ev, segment length must be a matrix.');
    end

    nbsegs = length( seg);
    if isscalar( overlap)
        overlap = repmat( overlap, 1);
    else
        if length( overlap) ~= nbsegs
            error( 'overlap and seg must be vectors of the same length.');
        end
    end

    pad_seg = 'no';
    center_shift = 'no';
    mat_epochs = cell( nbsegs, 1);
    for sidx = 1 : nbsegs
        buf = buffer( (1 : nb_pnts), seg( sidx), overlap( sidx), 'nodelay');
        if strcmpi( pad_seg, 'yes')
            % pad the last not-full segment with zeros
            if any( buf(:) == 0)
                buf( buf == 0) = nb_pnts + 1;
                if strcmpi( center_shift, 'yes')
                    % shift real data samples to center
                    shift_amount = floor( seg - sum( buf( :, end) < nb_pnts+1));
                    buf( :, end) = circshift( buf( :, end), shift_amount);
                end
            end

        else
            buf( :, ~all( buf > 0)) = [];
        end

        if isempty( buf)
            mat_epochs{ sidx} = [];
            continue;
        end

        for chnl_idx = 1 : nb_chnls
            chnl_data = mat( chnl_idx, :);
            chnl_data( nb_pnts+1) = 0;
            mat_epochs{ sidx}( chnl_idx, :, :) = chnl_data( buf);
        end
    end

    if length( seg) == 1
        mat_epochs = mat_epochs{1};
    end

else
    % epoching according to event markers
    nb_evs = length( ev);

    if iscell( seg)
        nbsegs = length( seg);
        mat_epochs = cell( nbsegs, 1);
        for sidx = 1 : nbsegs
            cur_seg = seg{ sidx};
            valid_segment( nb_pnts, ev, cur_seg);
            % cur_seg = cur_seg(1) : cur_seg(2);
            % tmp = arrayfun( @(x) mat( :, x + cur_seg), ev(:), 'un', 0);
            % mat_epochs{ sidx} = cat( 3, tmp{:});
            seg_len = diff( cur_seg) + 1;
            mat_epochs{ sidx} = nan( nb_chnls, seg_len, nb_evs);
            for k = 1 : nb_evs
                mat_epochs{ sidx}( :, :, k) = mat( :, ev( k) + (cur_seg(1) : cur_seg(2)));
            end
        end

    else
        valid_segment( nb_pnts, ev, seg);
        % seg = seg(1) : seg(2);
        % mat_epochs = arrayfun( @(x) mat( :, x + seg), ev(:), 'un', 0);
        % mat_epochs = cat( 3, mat_epochs{:});
        seg_len = diff( seg) + 1;
        mat_epochs = nan( nb_chnls, seg_len, nb_evs);
        for k = 1 : nb_evs
            mat_epochs( :, :, k) = mat( :, ev( k) + (seg(1) : seg(2)) );
        end
    end
end % event

end % function


function valid_segment( nb_pnts, ev, seg)
if length( seg) ~= 2
    error( 'Segment length must be a 1x2 or 2x1 vector.');
end

if any( ev + seg(1) < 1)
    error( 'Not enough samples at beginning of time series.');
end

if any( ev + seg(2) > nb_pnts)
    error( 'Not enough samples at the end of time series.');
end
end
