function [dtype, stats] = dtype_check(vals)

% automatically choose the appropriate to-be-written file data type based
% on images values (vals), and return the appropriate SPM dtype value
%
% see also: spm_type.m
%
% UINT8 - integers from 0 to +255  
% INT16 - integers from -32768 to +32767 
% INT32 - integers from -2.1M to +2.1M 
% FLOAT32 - single precision float 
% FLOAT64 - double precision float 
%
% RJE | 2017.11.09

vals = vals(:);

% get unique to save time
vals_u = unique(vals);

% how many unique?
val_u = numel(vals_u);

% NaNs present?
nan_check = sum(isnan(vals_u));

if nan_check == 0
    is_nan = 0;
else
    is_nan = 1;
end

% min value?
val_min = min(vals_u);

% max value?
val_max = max(vals_u);

% integers?
rem = mod(vals_u,1);
if max(rem) == 0 && min(rem) == 0
    is_int = 1;
else
    is_int = 0;
end

%% decision time

if is_int == 1
    if val_min >= 0         && val_max <= 2^8-1  && is_nan == 0
        % UINT8
        dtype = [2 0];
    elseif val_min >= -2^15 && val_max >= 2^15-1 && is_nan == 0
        % INT16
        dtype = [4 0];
    elseif val_min >= -2^31 && val_max <= 2^31   && is_nan == 0
        % INT32
        dtype = [8 0];
    else
        % we must have NaNs, which requires single-precision float (FLOAT32)
        dtype = [16 0];
    end
elseif is_int == 0
    % we have real numbers
    dtype = [16 0]; % use single-precision float (FLOAT32)
end


%% output
dtype = dtype;
stats.min = val_min;
stats.max = val_max;
stats.val_u = val_u;
