function [output] = imcalc(varargin)

%
% A batch script to do SPM imcalc transformations
%
% 1. Make sure the directory containing the program is added to your MATLAB path
% 2. At the prompt, type: imcalc
% 3. Select the list of files to transform
% 4. Select the operation to perform
%    (see details at: http://robjellis.net/tools/imcalc_documentation.pdf)
% 5. Each file will be written as a .nii to the same directory as the
%    original (except summed images, which will be written to the directory specified by the user), 
%    appended by a suffix which is the actual transformation applied
% 6. The transformed images can then be used in SPM design files, as would
%    regular .con files. 
% 
% For questions or feedback: email robjellis@gmail.com
%
% Copyright (C) 2011 by Robert J Ellis | http://robjellis.net

%
%    This program is a set of image calculator for MR images, and is 
%    made available the neuroimaging community as copyright freeware. 
%
%    You can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.  
%    See the GNU General Public License for more details.
%
%    A copy of the GNU General Public License may be found at:
%    http://www.gnu.org/licenses/gpl.html
%          
%    NOTES:
%
%    * Portions of the code call SPM functions (http://www.fil.ion.ucl.ac.uk/spm), 
%      which is also released under the GNU General Public License.  
%      SPM must be installed for imcalc to work


% must start here to get imcalc at the top of the path!
spm_path_fix('imcalc.m'); % spm_path_fix is by RJE

% get save date of this file name
loc = which('imcalc.m'); % note: must have a unique name and NOT an internal variable name

file_info = dir(loc); % note: there cannot be a variable "dir" active in Matlab or this won't work
save_date = file_info.date;
save_date = save_date(1:11); % trim off the timestamp for clarity



%% launch the program
if nargin == 1
    % skip right to it
    method = varargin{1};
else

    % get spm version and revision
    [vv, rr] = spm('ver');
    
fprintf('  _                     _      ')
fprintf('\n (_)                   | |     ')
fprintf('\n  _ _ __ ___   ___ __ _| | ___ ')
fprintf(['\n | | ''_ ` _ \\ / __/ _` | |/ __|   Image calculations and transformations (using ' vv ' rev. ' rr ')']')
fprintf(['\n | | | | | | | (_| (_| | | (__    This software version: ' save_date ]')
fprintf('\n |_|_| |_| |_|\\___\\__,_|_|\\___|   (C) Robert J Ellis (http://tools.robjellis.net)')  
    
% Nov 2017

%% primary method choice

M1 = input(['\n\n Choose a category of options to perform, or hit <ENTER> to quit: ',...
'\n     [1] Binarize/Threshold operations',...
'\n     [2] Mathematical/Statistical operations',...
'\n     [3] Other image transformations',...
'\n     [4] Masks and Masking',...
'\n     [5] Combine multiple images into one',...
'\n     [6] Split cluster image(s)',...
'\n     [7] User-entered expression (group-wise or image-wise)',...
'\n     --> ']);

if isempty(M1)
   % remove copies of SPM functions from the path to avoid errors
   remove_subfolder('imcalc.m','spm_fxns')
      
   fprintf('    Goodbye! \n\n');
   return 
   
elseif M1 < 0 || M1 > 7
    fprintf(' Error: Input not recognized. Terminating. \n')
    return
end

%% secondary method choice
M2 = []; % leave this here

fprintf('\n Choose a specific option to perform: ');

if M1 == 1 % binarize/threshold operations
    M2 = input(['',...
    '\n     [1] Binarize all non-zero voxels to 1.0 <ENTER>',...
    '\n     [2] Threshold (+ optional binarize)',...
    '\n     --> ']);  

    if isempty(M2)
        M2 = 1;
    end

elseif M1 == 2 % math and statistics
    M2 = input(['',...
    '\n     [1] Flip sign of all non-zero voxels',...
    '\n     [2] Transform T-values to Z-values',...
    '\n     [3] Winsorize (cap) using raw thresholds',...
    '\n     [4] Winsorize (cap) using percentile thresholds',...
    '\n     [5] Z-score transform voxels using image mean and SD',...
    '\n     [6] Homotopic voxel comparisons',...
    '\n     --> ']);    

elseif M1 == 3 % other image transformations
    M2 = input(['',...
    '\n     [1] Replace zeros or NaNs with another value',...
    '\n     [2] Flip image horizonally',...
    '\n     [3] Write single voxels to image(s)',...
    '\n     [4] Pad image with extra voxels',...
    '\n     [5] Linear XYZ translation of non-zero voxels',...
    '\n     --> ']);  

elseif M1 == 4 % masking
    M2 = input(['',...
    '\n     [1] Apply an inclusive or exclusive mask',...
    '\n     [2] Write hemisphere masks from a source image',...
    '\n     --> ']);  

elseif M1 == 5 % combine images (e.g., create cluster image)
    M2 = input(['',...
    '\n     [1] Operations on image pairs (add, sub, mult, div, etc.)',...
    '\n     [2] Create a cluster image',...
    '\n     [3] Threshold -> Binarize -> Sum (-> Rescale)',...
    '\n     --> ']);  
    
elseif M1 == 6 % split cluster image
    M2 = 1; % will yield "6.1"
elseif M1 == 7 % user-entered expression
    M2 = 1; % will yield "7.1"
end

if isempty(M2)
    return
end

% get the final choice
meth = str2num([num2str(M1) '.' num2str(M2)]); % will yield M1.M2

% *******************************************
% remap matrix of [old new] options
remap = [
    0   7.1;
    1   1.1;
    2   1.2;
    3   5.3;
    4   NaN;
    5   NaN;
    6   5.1;
    7   2.1;
    8   3.2;
    9   4.1;
    10  2.2;
    11  2.3;
    12  2.5;
    13  2.4;
    14  NaN;
    15  3.3;
    16  5.2;
    17  6.1;
    18  3.5;
    19  NaN;
    20  4.2;
    21  2.6;
    22  3.1;
    23  NaN;
    24  3.4
    
    ];
% *******************************************

% remmap it to older options
new = remap(:,2);
old = remap(:,1);
    
method = old(new == meth);

 
end

% leave this here to make code work
grpmat = 2; 

% get current directory
curdir = pwd;

%% method 0
if method == 0
   groups = input(' How many groups of images? \n   [1] One <ENTER> \n   [2] Two (pairwise operations) \n   --> ');
   if isempty(groups)
       groups = 1;
   end


   if groups == 1
       % select the images
       [files1] = spm_select([1 Inf],'image','Select image(s) to transform:',[],pwd,'.*'); 
       numf1 = size(files1,1);
       
       if numf1 == 0
           fprintf('    Error: No files selected. Terminating. \n\n')
           return
       end
       
       % do we have 4D data?
       % details on 4D volumes:
       % https://en.wikibooks.org/wiki/SPM/Working_with_4D_data
       % update 2017.11.03 - don't worry about 4D calculations for now
   
%       grpmat = input(' Calculation on \n   [1] 4D volume or \n   [2] each 3D volume <ENTER> \n   --> ');
%        if isempty(grpmat)
%            grpmat = 2;
%        end

        grpmat = 2;

   elseif groups == 2  
        [files1] = spm_select([1 Inf],'image','Select Group 1 images (e.g., "cond A") to transform:');
        numf1 = size(files1,1);
        
        if numf1 == 0
           fprintf('    Error: No files selected. Terminating. \n\n')
           return
        end
        
        files2 = spm_select([1 numf1],'image','Select Group 2 images (e.g., "cond B") to transform:');
        numf2 = size(files2,1);
        
        if numf1 ~= numf2
           fprintf(' Error: Number of files in Group 1 and Group 2 are not identical. Please respecify. \n\n')
           return
        end        
         
   end
   
%   % make sure to replace zeros ...
%   rep_zeros = input(' Replace zero-valued voxels with NaNs (otherwise inaccurate mean, std, etc.): \n   [1] Yes <ENTER> \n   [2] No \n   --> ');
   
%   if isempty(rep_zeros)
       rep_zeros = 1; % do this by default ...
%   end
   
   tag = ''; % leave this here as default
                  
   if groups == 1
       if grpmat == 1 
           expr = input([' Operation to perform at each voxel (ignoring NaNs): ',...
                        '\n   [0] User-entered expression on the volume (v) ',...
                        '\n   [1] Mean of all images ', ...
                        '\n   [2] Standard deviation of all images', ...
                        '\n   [3] Coefficient of variation (std / mean) ',...
                        '\n   [4] One-sample T-test ', ...
                        '\n   --> ']);
           if expr == 1
               expr = 'nanmean(V,4)';
               tag = '_mn';
           elseif expr == 2
               expr = 'nanstd(V,[],4)'; % sample standard deviation using N - 1
               tag = '_sd';
           elseif expr == 3
               expr = 'nanstd(V,[],4) ./ nanmean(V,4)';
               tag = '_cv';
           elseif expr == 4
               expr = 'nanmean(V,4) ./ (nanstd(V,[],4) / sqrt(numf1))';
               tag = '_t';
           elseif expr == 0
               expr = input(' Enter the operation to perform on volume V - e.g., max(V,[],4) ','s'); 
           end

        newdir = spm_select(1,'dir','Select the target directory for the output image:');
        newdir = char(newdir);
           % make sure it has a filesep
           if newdir(end) ~= filesep
               newdir = [newdir,filesep];
           end
        
        cd(newdir) 
   
       elseif grpmat == 2
           expr = input(' Enter the expression (no quotes) to perform \n on each image volume (v) (e.g., v == 10): ','s'); 
           
       end       
   elseif groups == 2
       expr = input(' Enter the expression (no quotes) to perform \n on each *pair* of image volumes (e.g., v1 + v2): ','s');
   end
   
   exprname = input(' Enter a suffix to write to the image: ','s');
   exprname = [exprname tag];
end % method == 0

%% next set of methods
this_set = [1 4 5 7 8 9 10 11 12 13 14 18 21 22 24];
if ismember(method,this_set)
   groups = 1;
   files1 = spm_select([1 Inf],'image','Select images to transform:'); 
   numf1 = size(files1,1);
   
   if numf1 == 0
       fprintf('    Error: No files selected. Please respecify. \n\n')
       return
   end

   if method == 4
      binit = input('\n Binarize images before summing?: [1] yes; [2] no: ');
   else
      binit = 1;
   end
   
   if method == 10
      df = input(' Enter the d.f. common to *all* to-be transformed images: ');  
   end
  
   if method == 22
      
      rpwith = input('\n Replacement to perform: \n   [1] 0 -> NaN \n   [2] NaN -> 0 \n   [3] 0 -> Max value \n   [4] NaN -> Max value \n   [5] 0 -> [Other value] \n   [6] NaN -> [Other value] \n   --> ');
      if rpwith == 1
          rwval = NaN;
      elseif rpwith == 2
          rwval = 0;
      elseif rpwith == 3 || rpwith == 4
          % defined during each file read
      elseif rpwith == 5 || rpwith == 6
          rwval = input(' Enter the desired value: ');
      end
   end
   
elseif method == 2 || method == 3
   groups = 1;
   files1 = spm_select([1 Inf],'image','Select file(s) to transform:'); 
   numf1 = size(files1,1);
   
   if numf1 == 0
       fprintf('    Error: No files selected. Please respecify. \n\n')
       return
   end   
   
   nthr = input('\n Thresholds at which to binarize: \n   [1] Identical across images <ENTER> \n   [2] Unique for each image \n   -->  ');
   
   if isempty(nthr)
       nthr = 1;
   end

    % percentile or raw?
    use_prc = input(' Threshold type: \n   [1] Manual <ENTER> \n   [2] Percentile-based \n   --> ');

    if isempty(use_prc)
       use_prc = 1;
    end
    
    if use_prc == 1
        if nthr == 1
           vec2 = input(['\n Enter threshold to use for all ' num2str(numf1) ' images: ']);
           vec = zeros(numf1,1) + vec2; 
        else
           vec = input(['\n Enter a vector of ', num2str(numf1), ' thresholds, in [ ]: ']); 
        end
    elseif use_prc == 2
        if nthr == 1
           vec2 = input(['\n Enter percentile threshold (0 to 100) to use for all ' num2str(numf1) ' images: ']);
           vec = zeros(numf1,1) + vec2; 
        else
           vec = input(['\n Enter a vector of ', num2str(numf1), ' percentile thresholds (0 to 100), in [ ]: ']);  
        end
    end
    
    nvec = numel(vec);
    
    if nvec ~= numf1
       % vector is not the same size as number of files
       fprintf('\n Error: the vector of thresholds and the number of files does not match. Respecify.\n\n');
       return
    end  
    
    UorD = input(['\n Retain values that are: ',...
             '\n   [1]  > threshold <ENTER>',...
             '\n   [2] >= threshold ',...
             '\n   [3]  < threshold ',...
             '\n   [4] <= threshold ',...
             '\n   [5] == threshold ',...
             '\n   -->  ']);
   
   if isempty(UorD)
       UorD = 1;
   end
   
   % switch between methods
   if method == 2
       thrit = input('\n Method: \n   [1] Threshold (preserve values) \n   [2] Binarize all non-zero voxels to 1.0 values <ENTER> \n   -->  ');
   elseif method == 3
       thrit = 2; % for simplicity; we don't care about raw values, only binarized
       
       % do we want to do a final rescaling step?
       rescale_it = input(' Rescale final image from 0.0 to 1.0? \n   [1] Yes \n   [2] No <ENTER> \n   --> ');
       
       if isempty(rescale_it)
           rescale_it = 2;
       end
   end
   
   if isempty(thrit)
       thrit = 2;
   end
   
  
      
elseif method == 6
   groups = 2;
   input_ch = input('\n Input style: \n   [1] a single pair of images \n   [2] separate Group 1 and Group 2 images \n   -->  ');
   
   if input_ch == 1
        files = spm_select([1 2],'image','Select a pair of files:');
        files1 = files(1,:);
        files2 = files(2,:);
        numf1 = 1;
        numf2 = 1;
   elseif input_ch == 2
        files1 = spm_select([1 Inf],'image','Select Group 1 files (e.g., "cond A") to transform:'); 
        numf1 = size(files1,1);
        
        if numf1 == 0
           fprintf('    Error: No files selected. Please respecify. \n\n')
           return
           
        end
        
        files2 = spm_select([1 numf1],'image','Select Group 2 files (e.g., "cond B") to transform:');
        numf2 = size(files2,1);
       
   end

       % check to see if same number of files
       if numf1 == numf2
           % OK
       else
           fprintf('\n Error: The number of Group 1 and Group 2 files does not match. Please respecify. \n\n');
           return
       end
   
   op = input('\n Select the desired operand: \n   [1] G1 + G2; \n   [2] G1 - G2; \n   [3] G1 .* G2; \n   [4] G1 ./ G2: \n   [5] min(G1,G2) \n   [6] max(G1,G2) \n   [7] null(G1,G2) \n   [8] mean(G1,G2) \n   -->  ');
   
   if op == 3 || op == 4
        pow = input('\n Select the *powers* to scale G1 and G2 to, respectively, in [ ] (default = [1 1]): '); % e.g., [1 1/2] it we are dividing a beta map by variance map
   else
      pow = [1 1];
   end
   
   newname = input('\n Enter a short descriptive name for this condition (no spaces), or hit <ENTER> for none: ','s');
   
   if isempty(newname)
       newname = '';
   else 
       newname = ['_' newname];
   end
   
   if op == 1
      op1 = strcat(newname,'_add'); 
   elseif op == 2
      op1 = strcat(newname,'_sub'); 
   elseif op == 3
      op1 = strcat(newname,'_mlt'); 
   elseif op == 4
      op1 = strcat(newname,'_div');
   elseif op == 5
      op1 = strcat(newname,'_min');
   elseif op == 6
      op1 = strcat(newname,'_max');
   elseif op == 7
      op1 = strcat(newname,'_null');
   elseif op == 8
      op1 = strcat(newname,'_mean');
   end
   
elseif method == 15
    files1 = spm_select(1,'image','Select a template image to define coordinate space:');
    numf1 = size(files1,1);
    
    if numf1 == 0
       fprintf('    Error: No files selected. Please respecify. \n\n')
       return
    end
    
    XYZmm = input('\n\n Enter a [N x 3] matrix of voxels (in mm coordinates), in [ ], or "0" to enter a range of values: '); % rje edit 
    if numel(XYZmm) == 1 
        % make a grid
        xin = input('  X coordinates; "0" yields [ -72:4:72]): ');
        if xin == 0
           xin = [-72:4:72];
        end
        nx = numel(xin);
        yin = input('  Y coordinates; "0" yields [-104:4:72]): ');
        if yin == 0
           yin = [-104:4:72]; 
        end
        ny = numel(yin);
        zin = input('  Z coordinates; "0" yields [ -60:4:76]): ');
        if zin == 0
           zin = [-60:4:76];  
        end
        nz = numel(zin);
        % now we have to create the list
        XYZmm = nan(nx*ny*nz,3);
        ct = 1;
        for xx = 1:nx
            for yy = 1:ny
                for zz = 1:nz
                    crd = [xin(xx) yin(yy) zin(zz)];
                    XYZmm(ct,:) = crd;
                    ct = ct+1;
                end
            end
        end
        
    end
    numcoord = size(XYZmm,1);
    
    imval = input('\n Enter the intensity value to write at each voxel (<ENTER> for 1.0): ');
    
    if isempty(imval)
        imval = 1;
    end
    
    if numcoord == 1
        numim = 1;
    elseif numcoord > 1
        
        numim = input(['\n Write \n   [1] A single image or \n   [2] ', num2str(numcoord), ' individual images \n   --> ']);
    end
    
    if numim == 2
        numf1 = numcoord; % to make the loop work
    elseif numim == 1
        numf1 = 1;
        imsuf = input('\n Suffix to append to image (<ENTER> to skip): ','s');
        
        if isempty(imsuf)
            imsuf = '';
        else
            imsuf = ['_' imsuf];
        end
    end
    
    % read in the template volume
    tvol = spm_vol(files1);
    vt = spm_read_vols(tvol);
    dimm = size(vt);
    xdim = dimm(1); % size of x-dimension
    
    v1mat=tvol.('mat');
    xflip = sign(v1mat(1,1));    % -1 means values read in from L to R

    origin2 = inv(v1mat);
    origin2 = origin2(1:3,4);
    origin2 = origin2';
    origin = origin2; % rje confirmed this is correct
    
    % we need to get from mm to vox
    vox = [abs(v1mat(1,1)), abs(v1mat(2,2)), abs(v1mat(3,3))];
    
    XYZvox(:,1) = XYZmm(:,1) / vox(1);
    XYZvox(:,2) = XYZmm(:,2) / vox(2);
    XYZvox(:,3) = XYZmm(:,3) / vox(3);
    
    % ok, we better make sure we have integers
    XYZvox = round(XYZvox);
    
    % now we need to get coordinates into matrix space
      XYZmat = zeros(size(XYZmm));
      XYZmat(:,2) = XYZvox(:,2) + origin(2); % y coord are fine
      XYZmat(:,3) = XYZvox(:,3) + origin(3); % z coord are fine

      % note: this may not agree with how MRIcron displays points, but it
      % *will* write correctly into Matlab space!
    if xflip == 1
        XYZmat(:,1) = XYZvox(:,1) + origin(1); % no problem
    elseif xflip == -1
        XYZmat(:,1) = (xdim - XYZvox(:,1)) - (xdim - origin(1)); % have to count in from the upper edge of x
    end
   
elseif method == 16
    files1 = spm_select([1 Inf],'image','Select images to sum into cluster image:');
    numf1 = size(files1,1);

    if numf1 == 0
       fprintf('    Error: No files selected. Please respecify. \n\n')
       return
    end
    
    indx = input(['\n Paste in a [' num2str(numf1) ' x 1] vector of index values: ']);
    indx = indx(:); % just to be sure
    if numel(indx) == numf1
        % OK
    else
        fprintf(' The number of input values does not match the number of files. Respecify. \n\n');
        return
    end

    newname = input('\n Enter a short descriptive name for this image (no spaces), or hit <ENTER> for none: ','s');
   
       if isempty(newname)
           newname = '';
       end
       
elseif method == 17
    % Updated Nov 2017: now can handle cluster images in batch mode
    
    files1 = spm_select([1 Inf],'image','Select cluster image(s) to split:');
    tot_files = size(files1,1);

    if tot_files == 0
       fprintf('    Error: No files selected. Please respecify. \n\n')
       return
    end    
    
    % if we have multiple cluster images selected, assume we run everything
    % on full auto mode (i.e., no subsequent masking, etc.)
    
    if tot_files == 1
        split_type = input('\n Split method: \n   [1] auto <ENTER> \n   [2] manual \n   --> ');
    else
        split_type = 1;
    end

    if isempty(split_type)
        split_type = 1;
    end
    
    if split_type == 2
        regs = input('\n Enter a [N x 2] matrix of all L and R values to split. \n If L and R have the same value, then enter [1 1; 2 2] and so on: ');
            
        % to make later code work   
        numf1 = size(regs,1);
        num_type = 0; 
        cont = 1;
        
    elseif split_type == 1
        num_type = input('\n Round voxel values in each image? \n   [0] No, leave unaltered <ENTER> \n   [1] Yes, to integers  \n   [2] Yes, to one decimal \n   [3] Yes, to two decimals \n   --> ');
        
        if isempty(num_type)
            num_type = 0;
        end
    end
    
    keep_reg = input('\n For the split images: \n   [1] Preserve original voxel values <ENTER> \n   [2] Binarize all non-zero voxels to 1.0 \n   --> ');

    if isempty(keep_reg)
        keep_reg = 1;
    end
          
    % now we add the new loop
    for t = 1:tot_files

        if t > 1
            progressbar(0,0)
        elseif t == 1
            progressbar(0)
        end

        % read in the first volume
        file1 = files1(t,:);

        v1n = spm_vol(file1);
        v1  = spm_read_vols(v1n);  % the actual volume
        v1(isnan(v1)) = 0;   % replace NaN if they are present

        if num_type == 0
           % don't need to do anything to v1
        elseif num_type == 1
           v1 = round(v1);
        elseif num_type == 2
           v1 = (round(v1)*10)/10;
        elseif num_type == 3
           v1 = (round(v1)*100)/100; 
        end
        v1all = v1(:);
        regs = unique(v1all); % get the unique values in this image that we will later split apart
        regs = regs(regs>0); % don't care about zero
        numf1 = numel(regs); % to make the next loop work 
          
        if split_type == 1
            % careful with lots of unique voxels values (i.e., maybe we should have rounded the image ...)
            if numf1 > 300
                if tot_files == 1
                    fprintf([' *Warning*: the number of to-be-written images is ' num2str(numf1) '.']);
                elseif tot_files > 1
                    fprintf([' *Warning*: the number of to-be-written images for file ' num2str(t) ' is ' num2str(numf1) '.']);
                end

                cont = input(' Do you wish to continue? \n   [1] Yes \n   [2] No <ENTER> \n   --> ');
                if isempty(cont)
                    cont = 2; % just so we don't accidentally start writing tons of files
                end
            else
                % everything is fine, assume write all images
                cont = 1;
            end
        end
              
        % the splitting loop
        if cont == 1
            if t > 1
                progressbar(t/tot_files,[])
            end
            
            if tot_files > 1
                fprintf(['\n Working on file ' num2str(t) '...'])
            end

            for f = 1:numf1
                if t > 1
                    progressbar([],f/numf1)
                elseif t == 1
                    progressbar(f/numf1)
                end

               reg = regs(f);  
               vol = (v1 == reg);       
               if keep_reg == 1
                    vol = vol * reg; % now we actually encode the region to be clear
               elseif keep_reg == 2
                    % keep the binarized version
               end

                % prepare to write the file
                if isempty(findstr(file1,',1'))
                % OK 
                else
                % get rid of the suffix
                file1 = file1(1:findstr(file1,',1')-1);
                end    

                % check for appropriate dtype
                dtype = dtype_check(vol(:));
                v1n.('dt') = dtype;

                % write this volume
                v1n.('fname') = strcat(file1, '_reg_', num2str(reg), '_bin.nii'); 
                v1n = spm_write_vol(v1n,vol);

            end % f loop
        else
            % we skip this file and move on to the next
            fprintf(['\n Skipping file ' num2str(t) '. \n'])
        end
    end % t loop 
    
    progressbar(1)

elseif method == 19 % deprecated (code never written) - 2017.11.01
    files1 = spm_select([1 Inf],'image','Select all images:');
    numf1 = size(files1,1);
    n = numf1;
    % now set up the venn regions. Not including the "universe", there are
    % (2^n - 1) intersection regions defined for n elements
    % (a positive value means add binarized images; a minus value means exclude
    % binarized images)
    
    venn = zeros(2^n - 1,n);
    
    if n == 2
       venn = [1 -2; -1 2]; 
    elseif n == 3
        
    elseif n == 4
        
    end 
    
elseif method == 20
    files1 = spm_select(1,'image','Select a single source image:');
    
    if size(files1,1) == 0
       fprintf('    Error: No files selected. Please respecify. \n\n')
       return
    end  
    
    numf1 = 1;
                
elseif method == 23
    files1 = spm_select([1 Inf],'image','Select the N cluster maps:');
    regs = input('\n Enter the vector of values to isolate, in [ ]: ');
    regs = regs(:);
    numf1 = size(regs,1); % to make the write loop work correctly
    
    nsubs = size(files1,1);
    
    if nsubs == 0
       fprintf('    Error: No files selected. Please respecify. \n\n')
       return
    end  
    
end % method


% =============================================
%% additional calculations
if method == 4 || method == 3 || method == 5 || method == 16 || method == 23
   
   [same_dir, direc] = dir_check(files1);
   
   if same_dir == 1
       % if all images are from the same directory, don't need to ask
       newdir = direc;
   else
       newdir = spm_select(1,'dir','Select the target directory for the output image(s):');
   end
   
    if isempty(newdir)
       fprintf('    Error: No directory selected. Please respecify. \n\n')
       return
    end   
    
   newdir = char(newdir);
   % make sure it has a filesep
   if newdir(end) ~= filesep
       newdir = [newdir,filesep];
   end
   cd(newdir) 
   
elseif method == 10
   if exist('norminv.m','file') == 2
       z_meth = input(' Use [1] Matlab norminv and tcdf or [2] spm_t2z: '); 
   else
       fprintf(' (Note: MATLAB stats toolbox not available. Using spm_t2z function.)\n\n ')
       z_meth = 2;
   end
     
elseif method == 11
   thrL = input(' Enter the lower threshold, L, or <ENTER> for no lower threshold. \n (Values < L will become L): ');
   thrU = input(' Enter the upper threshold, U, or <ENTER> for no upper threshold. \n (Values > U will become U): ');
   
   if isempty(thrL)
       thrL = -Inf;
   end
   
   if isempty(thrU)
       thrU = Inf;
   end
   
elseif method == 13
   prcL = input(' Enter the lower percentile threshold, L (0 to 100), or <ENTER> for no lower threshold. \n Values < prc(L) will become prc(L): ');
   prcU = input(' Enter the lower percentile threshold, U (0 to 100), or <ENTER> for no upper threshold. \n Values > prc(U) will become prc(U): ');
   
   if isempty(prcL)
       prcL = 0;
   end
   
   if isempty(prcU)
       prcU = 100;
   end
   

elseif method == 14
   nmult = input(' How many SDs to divide by? (e.g., 3): ');  
   
elseif method == 18
    shift = input(' Enter the desired translation, in *voxels*, for [x y z]: ');
    if numel(shift) ~= 3
      error('Error: requires 3 values in [ ].')
    end
    
elseif method == 21
    op = input('\n Select the desired operand to perform at all homotpic voxel locations (O = original; X = x-flipped): \n   [1] O - X (signed) \n   [2] X - O (signed) \n   [3] min(O,X) \n   [4] max(O,X); \n   [5] mean(O,X); \n   [6] compare values in O vs. X \n   --> ');
   
   if op == 1
      op1 = strcat('_hdiff_O-X'); 
   elseif op == 2
      op1 = strcat('_hdiff_X-O'); 
   elseif op == 3
      op1 = strcat('_hmin'); 
   elseif op == 4
      op1 = strcat('_hmax');
   elseif op == 5
      op1 = strcat('_hmean');
   elseif op == 6
      op1 = strcat('_comp');
   end

elseif method == 24
    padvox = input('\n Number of voxels to add on all sides of image (single number, or [x y z]: ');
    if numel(padvox) == 1
        padvox = [padvox padvox padvox];
    end
end

% =======================================
% general masks

if method == 9
    use_mask = 1;
    
elseif method == 17 % for splittting out a cluster image, assume we don't need a mask (makes code simplier)
    use_mask = 2;
    maskch = 1; % just to make code work
else
    use_mask = input('\n Apply an additional inclusive/exclusive mask? \n   [1] Yes \n   [2] No <ENTER> \n   --> ');
    if isempty(use_mask)
        use_mask = 2;
    end
    maskch = 1; % just to make code work
end

if use_mask == 1
    mask = spm_select(1,'image','Select inclusive/exclusive mask:',[],pwd,'.*');
    
    if isempty(mask)
       fprintf('    Note: No mask selected. Continuing ... \n\n')
       return
    end     
    
    % get mask directory
    maskdir = fileparts(mask);

    
    % first, let's make sure that the images and masks have the same
    % dimensions. We just check this once, and hope for the best.

    chkmask = strtrim(mask(1,:));
    cmv = spm_read_vols(spm_vol(chkmask));

    chkfile = strtrim(files1(1,:));
    cfv = spm_read_vols(spm_vol(chkfile));

    if sum(size(cmv) - size(cfv)) == 0
      % we are ok; use the original mask
     
    else
      % need to reslice

      fprintf(' \n Note: Mask has a different dimension than the target image(s). Reslicing ...');

      % masks will have the prefix "rs_"
      % assume nearest neighbor for simplicity!
      % some code borrowed from Tor Wager's reslice_imgs.m

        flags = struct('interp', 0, ... % rje: neearest neighbor
            'mask', 0, ...              % do not mask
            'mean', 0, ...              % do not write mean image
            'hold', -1, ...             % TW: i don't think this is used anymore
            'which', 1, ...             % reslice 2nd-nth only
            'wrap', [1 1 0]', ...       % the default (rje: SPM says [1 1 0] for fMRI)
            'prefix','rs_' ...           
            );

        P = str2mat(chkfile,mask);
        spm12_reslice(P, flags) % copy of spm_reslice from SPM12

        % use the resliced mask
        mask = strtrim(mask(1,:)); % cut out text whitespace in filename

        xx = find(mask == filesep);
        dd = mask(1:max(xx));
        ff = mask(max(xx)+1:end);
        mask = [dd 'rs_' ff];

    end
    
    fprintf(' Done.\n') % done with reslice
    
    % decide whether inclusive or exclusive
    maskch = input('\n Treat mask as \n   [1] Inclusive <ENTER> \n   [2] Exclusive \n   --> ');
    if isempty(maskch) 
        maskch = 1;
    end
    
       
elseif use_mask == 2
    % just use the full bounding box
    mask = files1(1,:);
    maskdir = fileparts(mask);
end

%% ready to go
if method == 17 && split_type == 1
    % code is handled above
else
    fprintf(' Working ...')  % for all options
    
    % Nov 2017: use progressbar.m (freeware)
    if numf1 > 1
    	progressbar(0)
    end
end

%% mask non-brain voxels   
vm = spm_read_vols(spm_vol(mask));
vm(isnan(vm)) = 0;  % replaces NaNs with zero
vm = vm ~= 0;  % binarize; all active voxels will now equal 1.0
vm = double(vm);
dimm = size(vm);

% now we have to "flip" the mask if exclusive
if maskch == 2
   vm = abs(vm - 1); % 0s become 1s and 1s become 0s
end

if use_mask == 2
    vm = ones(dimm(1), dimm(2), dimm(3));
end

% ==================================
%% THE LOOP

% ----------------
if method <= 1 || method >= 10

   for f=1:numf1
       
        if numf1 > 1
           progressbar(f/numf1)
        end
   
       if method ~= 0 && method ~= 15 && method ~= 17 && method ~= 23
           file1 = files1(f,:);
           v1n = spm_vol(file1);

            v1mat=v1n.('mat');
            xflip = sign(v1mat(1,1));    % -1 means values read in from R to L

            mat_inv = inv(v1mat);
            origin = mat_inv(1:3,4);
            origin = origin'; % rje confirmed this is correct

            vox = [abs(v1mat(1,1)), abs(v1mat(2,2)), abs(v1mat(3,3))];


           v1  = spm_read_vols(spm_vol(file1));  % the actual volume
           v1(isnan(v1)) = 0;   % replace NaN if they are present
           voldim = size(v1);

           if method == 8 || method == 20 || method == 21
               % confirm that the x-value of the origin is in the center of the
               % x-dimension

               x1 = 1: 1 :(origin(1)-1);
               x2 = voldim(1): -1 :origin(1)+1;
               if numel(x1) == numel(x2)
                  % OK
               else
                  fprintf('\n Warning: the x-value of the origin is not in the center of the image. Terminating. \n\n');
                  return
               end
           end

       elseif method == 15
           v1n = tvol;
           v1 = zeros(size(vm));

       elseif method == 0

           if groups == 1
                file1 = files1(f,:);
                v1n = spm_vol(file1); % always use the first group
                v  = spm_read_vols(spm_vol(file1));  % the actual volume

                if rep_zeros == 1
                    v(v == 0) = NaN; % replace 0s with NaN
                else
                    % do nothing; keep the zeros
                end

                if grpmat == 1

                    if f == 1
                       sizev = size(v);
                       V = nan(sizev(1),sizev(2),sizev(3),numf1);
                    end
                       V(:,:,:,f) = v;
                end

           elseif groups == 2
                file1 = files1(f,:);
                v1n = spm_vol(file1); % always use the first group
                v1  = spm_read_vols(spm_vol(file1));  % the actual volume

                if rep_zeros == 1
                    v1(v1 == 0) = NaN;
                else
                    % do nothing; keep the zeros
                end

                file2 = files2(f,:);
                v2  = spm_read_vols(spm_vol(file2));  % the actual volume

                if rep_zeros == 1
                    v2(v2 == 0) = NaN; % rje: error saying "v2 may be unused" is fine, because the "expr" takes care of it (and v1 is used elsewhere by other options)
                else
                    % do nothing; keep the zeros
                end            
           end

            % leave this here
            if f == 1
              zero_vox = NaN;
              k = 1;
            end      

           if grpmat == 1 && f == numf1 % only after all volumes are read in ...

              volx = eval(expr); % finally we evaluate the expression on the 4D volume
              vol = volx .* vm; % reapply the mask


              % check - are there any voxels?
              if sum(vol(:) > 0) == 0
                  zero_vox(k) = f;
                  k = k + 1;
              end
           elseif grpmat == 2
              volx = eval(expr); % evaluate the expression on v1 and v2, defined above
              vol = volx .* vm; % reapply the mask

                        % check - are there any voxels?
              if sum(vol(:) > 0) == 0
                  zero_vox(k) = f;
                  k = k + 1;
              end
           end


       end % method choice
   
   % ===== for all methods, we only do calculations based on voxels in the
   % implicit mask 
   if method ~= 0 && method ~= 17 && method ~= 23
        v1 = v1 .* vm;
   end    
        
   if method == 1
        vol = v1~=0;  % binarize all non-zero voxels to 0
        
   elseif method == 10
       
       if z_meth == 1
           % will only work if "stats" toolbox is installed
           zmeth = 'spm-t2z';
            vol = norminv(tcdf(v1,df)); % transforms t to z using two MATLAB functions
       elseif z_meth == 2
           zmeth = 'mat-t2zt';
            vol = rje_spm_t2z(v1,df);
       end
       
   elseif method == 11
       
       % apply lower threshold
       v1(v1 < thrL) = thrL;
       
       % apply upper threshold
       v1(v1 > thrU) = thrU;
       
       vol = v1; % reassign for later use; will mask non-brain voxels later

   elseif method == 12
       v = v1(v1~=0); % all non-zero values (so we don't have to use nanmean or nanstd); will be a vector
       m  = mean(v);
       sd = std(v);
       
       vol = (v1 - m) / sd;       
   
   elseif method == 13
       
       % note: let's assume we don't have a 4D volume here ...
       
       v = v1(v1~=0); % all non-zero values (so percentiles are not biased); will be a vector
       
       thrL = prctile_nist(v,prcL); % doesn't require stats toolbox
       thrU = prctile_nist(v,prcU);
       
       % apply lower threshold
       v1(v1 < thrL) = thrL;
       
       % apply upper threshold
       v1(v1 > thrU) = thrU;
       
       vol = v1; % reassign for later use; will mask non-brain voxels later
       
   elseif method == 14
       nonz = v1(v1~=0);
       vol = v1 / (nmult * std(nonz));
   
   elseif method == 15
       file1 = files1; % just copy the name
       
       if numim == 2
           coor = XYZmat(f,:); % a single voxel
           ocoor = XYZmm(f,:); % the original coordinate (mm)
           v1(coor(1),coor(2),coor(3)) = 1.0; % write this value to that voxel
           vol = v1;
       elseif numim == 1
           ct = 1;
         
           for g = 1:numcoord;
               coor = XYZmat(g,:); % a single voxel
               ocoor = XYZmm(g,:); % the original coordinate (mm)
               if vm(coor(1),coor(2),coor(3)) == 1 % for voxels in the mask ...
                  v1(coor(1),coor(2),coor(3)) = imval; % write this value to that voxel
                  coor_in(ct,:) = ocoor; % print this out for reference; variable will grow on each increase but it can't be helped
                  ct = ct + 1;
               end
                 
           end
           
           vol = v1; % this will contain all the coords
           %output.coor_in = coor_in;
    
       end
       
   elseif method == 16
       if f == 1     
            % create the output volume
            vol = zeros(size(vm));
       end
           
       tvol = v1~=0;  % binarize all non-zero voxels to 0
       tvol = tvol .* indx(f); % so that we have this value at all voxels
       
       vol = max(vol,tvol); % easiest way to push masks together without introducing new voxel intensities
       
   elseif method == 17
       if split_type == 2 % this stays here
              reg_pair = regs(f,:);
              reg = reg_pair(1); % for the write label  
              volL = (v1 == reg_pair(1)); % binarize matches to this value
              volR = (v1 == reg_pair(2));
              vol = volL + volR;
              
       elseif split_type == 1
             % code was moved 
       end
       
   elseif method == 18    
       vol = zeros(voldim(1),voldim(2),voldim(3)); % empty
       for i = 1:voldim(1)
           for j = 1:voldim(2)
               for k = 1:voldim(3)
                   
                   if v1(i,j,k) ~= 0
                      if xflip == 1
                         vol(i+shift(1),j+shift(2),k+shift(3)) = v1(i,j,k); 
                      elseif xflip == -1
                         vol(i-shift(1),j+shift(2),k+shift(3)) = v1(i,j,k);  
                      end
                   end   
               end
           end
       end
   
   elseif method == 20
      xmax = voldim(1);
      volL = zeros(voldim(1),voldim(2),voldim(3)); % empty
      volR = zeros(voldim(1),voldim(2),voldim(3)); % empty
      
      volF = ones(voldim(1),voldim(2),voldim(3));
      
      % manually do this for L(eft), R(ight), and F(ull)
     
      if xflip == 1
        volL(1:origin(1)-1,:,:) = 1;
        volR(origin(1)+1:xmax,:,:) = 1;
      
      elseif xflip == -1 % the opposite
        volL(origin(1)+1:xmax,:,:) = 1;
        volR(1:origin(1)-1,:,:) = 1;
      end
     
            
      % need to mask to template if one was selected
      volL = volL .* vm;
      volR = volR .* vm;
      volF = volF .* vm;
      
   elseif method == 21
      % note: doing this at the matrix level only works if the origin is at
      % the center of the x dimension (it should be for SPM data and masks)
      % This was already checked above.
      
      % copy of the original volume, for clarity
      v1o = v1;
      
      % flip the volume in the x-dimension
      v1x = flipdim(v1,1);
      
      % do the appropriate calculation
      
      if op == 1
         vol = v1o - v1x; 
      elseif op == 2
         vol = v1x - v1o;  
      elseif op == 3
         vol = min(v1o,v1x);  
      elseif op == 4
         vol = max(v1o,v1x);   
      elseif op == 5
         vol = (v1o+v1x) / 2; 
      elseif op == 6
          % note: this only makes sense for integer labeled masks (e.g., the AAL)
          % 1. compare the original and the flipped image at every voxel
          % 2. write a binarized mask that has a "1" if same, and "0" if different
          % 3. vol is then the truth mask times the original image
          tmask = v1o == v1x;
          tmask = double(tmask);
          vol = v1o .* tmask;
      end
      
   elseif method == 22
       if rpwith == 3 || rpwith == 4
           % need to get the max value
           rwval = max(max(max(v1)));
       end
           
       if rpwith == 1 || rpwith == 3 || rpwith == 5
           v1(v1 == 0) = rwval;
       
       elseif rpwith == 2 || rpwith == 4 || rpwith == 6
           v1(isnan(v1)) = rwval;
       end
       
       op1 = '_repl';
       vol = v1;

   elseif method == 23
       
       reg = regs(f);
       v1n = spm_vol(mask);     
       vol = zeros(dimm(1),dimm(2),dimm(3)); % reset for each reg value (dimm is from the mask)
       
       for s = 1:nsubs
          file1 = files1(s,:); 
          v1  = spm_read_vols(spm_vol(file1));  % the actual volume 
          rvol = v1 == reg;
          vol = vol + rvol; % read through each subject map each time, writing to this volume
       end
       
       vol = vol / nsubs; % to scale from 0 to 1
       vol = vol .* vm; % the implicit mask

   elseif method == 24
       vol = zeros(voldim(1)+padvox(1)*2,voldim(2)+padvox(2)*2,voldim(3)+padvox(3)*2); % if 3 x N values, 0 means no padding
       
       vol(1+padvox(1):padvox(1)+voldim(1),1+padvox(2):padvox(2)+voldim(2),1+padvox(3):padvox(3)+voldim(3)) = v1; 
       
       % move the origin; easy to do with the inverted matrix
       mat_inv(1:3,4) = mat_inv(1:3,4)+padvox(:); % just shift all values
       
       % now, reinvert the matrix, and put it back in the file header
       v1n.mat = inv(mat_inv);
       
   end % method choice


     if isempty(findstr(file1,',1'))
        % OK 
     else
        % get rid of the suffix
        file1 = file1(1:findstr(file1,',1')-1);
     end
     
     % volume header - if not already defined
     if method == 0
        if grpmat == 2 
           v1n.('fname') = strcat(file1, '_', exprname,'.nii');
        elseif grpmat == 1
           v1n.('fname') = strcat(newdir, '4D_matrix_calc_', exprname, '.nii'); 
        end
     elseif method == 1
        v1n.('fname') = strcat(file1, '_bin.nii');    % will write to the "file1" directory
     elseif method == 10 
        v1n.('fname') = strcat(file1,'_',zmeth,'.nii');    % will write to the "file1" directory
     elseif method == 11
        v1n.('fname') = strcat(file1, '_wins.nii');    % will write to the "file1" directory
     elseif method == 13
        v1n.('fname') = strcat(file1, '_prc.nii');  
     elseif method == 14
        v1n.('fname') = strcat(file1, '_std', num2str(nmult), '.nii'); 
     elseif method == 15
        if numim == 2
            v1n.('fname') = strcat(file1, '_ctr_', num2str(f),'_vox_', num2str(ocoor(1)), '_', num2str(ocoor(2)), '_', num2str(ocoor(3)), '.nii');
        elseif numim == 1
            v1n.('fname') = strcat(file1,'_all_coords', imsuf, '.nii');
        end
     elseif method == 16
        v1n.('fname') = strcat(newdir, 'cluster_image_', newname, '.nii'); 
     elseif method == 17
        if split_type == 2
            v1n.('fname') = strcat(file1, '_reg_', num2str(reg), '.nii');  
        elseif split_type == 1
            % code has been moved up
        end
     elseif method == 18
        v1n.('fname') = strcat(file1, '_shift_', num2str(shift(1)), '_', num2str(shift(2)), '_', num2str(shift(3)), '.nii'); 
     elseif method == 21
        v1n.('fname') = strcat(file1,op1,'.nii'); 
     elseif method == 22
        v1n.('fname') = strcat(file1,op1,'.nii');
     elseif method == 23
        v1n.('fname') = strcat(newdir,'prop_reg_',num2str(reg),'.nii');
     elseif method == 24
        if var(padvox) == 0
            % values are all the same
            v1n.('fname') = strcat(file1, '_pad_xyz_', num2str(padvox(1)), '_vox.nii');
        else
            % different x y z padding, so write that to the file
            v1n.('fname') = strcat(file1, '_pad_', num2str(padvox(1)), '_', num2str(padvox(2)), '_', num2str(padvox(3)), '_vox.nii');
        end
         
     end
     
     
     % ============================================
     %% write the volumes
     if method == 16 && f == numf1
         
         % check for appropriate dtype
         dtype = dtype_check(vol(:));
         v1n.('dt') = dtype;
         
         v1n = spm_write_vol(v1n,vol);
         
     elseif method == 20
         % write three images     
         if isempty(findstr(file1,',1'))
            % OK 
         else
            % get rid of the suffix
            file1 = file1(1:findstr(file1,',1')-1);
         end

         v1n.('pinfo') = [0;0;0]; % so that binarized  == 1
         v1n.('fname') = strcat(file1, '_mask_left.nii'); 
         v1n.('dt') = [2 0]; % since a binarized image
         spm_write_vol(v1n,volL);

         v1n.('fname') = strcat(file1, '_mask_right.nii'); 
         v1n.('dt') = [2 0]; % since a binarized image
         spm_write_vol(v1n,volR);

         v1n.('fname') = strcat(file1, '_mask_full.nii'); 
         v1n.('dt') = [2 0]; % since a binarized image
         spm_write_vol(v1n,volF);
     
     % need to reset pinfo for all binarized images
     elseif method == 1 || method == 2
         
         v1n.('dt') = [2 0]; % since a binarized image
         v1n.('pinfo') = [0;0;0]; % so that binarized  == 1
         
         v1n = spm_write_vol(v1n,vol);  
     
     elseif method == 3
        % check for appropriate dtype
        dtype = dtype_check(vol(:));
        v1n.('dt') = dtype;
        v1n = spm_write_vol(v1n,vol); 
        
     elseif method == 4 || method == 5
         % deprecated 
         
     elseif method == 17
         if split_type == 2
            % check for appropriate dtype
            dtype = dtype_check(vol(:));
            v1n.('dt') = dtype; 
            
            v1n = spm_write_vol(v1n,vol);
         elseif split_type == 1
             % code has been moved up
         end
          
     elseif method == 24
         v1n.dim = [voldim(1)+padvox(1)*2 voldim(2)+padvox(2)*2 voldim(3)+padvox(3)*2];
         
         % check for appropriate dtype
         dtype = dtype_check(vol(:));
         v1n.('dt') = dtype;
         
         v1n = spm_write_vol(v1n,vol);        
 
     else % all other cases
        if grpmat == 2 % default; individual images are written
            % check for appropriate dtype
            dtype = dtype_check(vol(:));
            v1n.('dt') = dtype;
            
            v1n = spm_write_vol(v1n,vol);  
        end
     end
    
   end % for f = 1:numf1
   
   % for method 0 using the matrix input
   if grpmat == 1
       
      % check for appropriate dtype
      dtype = dtype_check(vol(:));
      v1n.('dt') = dtype;
    
      v1n = spm_write_vol(v1n,vol);  
   end
   
end % method list


%% below are distinct loops with special operations and volume writing
% ----------------
if method == 2
   % a unique value to binarize each image at
   
   for f=1:numf1
      
        if numf1 > 1
           progressbar(f/numf1)
        end

   file1 = files1(f,:);
   thr = vec(f);
   v1n = spm_vol(file1);
   v1  = spm_read_vols(spm_vol(file1));  % the actual 3D volume
   
   % pay attention to (1) manual vs percentile and (2) binarize vs preserve values       
   if thrit == 1
       vnew = v1; % copy of original values
   else
       vnew = ones(size(v1)); % simple mask of ones
   end

    if UorD == 1
        vol = vnew .* (v1 > thr);  
        suf = '_binA_';            
    elseif UorD == 2
        vol = vnew .* (v1 >= thr); 
        suf = '_binA_'; 
    elseif UorD == 3
        vol = vnew .* (v1 < thr);  
        suf = '_binB_';
    elseif UorD == 4
        vol = vnew .* (v1 <= thr);  
        suf = '_binB_';
    elseif UorD == 5
        vol = vnew .* (v1 == thr);  
        suf = '_binE_';
    end      


     % mask (with mask image or itself)
     vol = vol .* vm; % important, so we don't retain non-brain voxels
   
     % write the volume
          if isempty(findstr(file1,',1'))
            % OK 
         else
            % get rid of the suffix
            file1 = file1(1:findstr(file1,',1')-1);
         end
     v1n.('fname') = strcat(file1, suf ,num2str(thr),'.nii');    % will write to the "file1" directory

     % check for appropriate dtype
     dtype = dtype_check(vol(:));
     v1n.('dt') = dtype;
     
     v1n = spm_write_vol(v1n,vol);
 
   end % file loop
   
   
end

% ----------------
if method == 6
 
% to tally retained percentage for v1 and v2
retained_vox = zeros(numf1,2);
    
    for f=1:numf1
        if numf1 > 1
           progressbar(f/numf1)
        end
        
        file1 = files1(f,:);
        file2 = files2(f,:);
        v1n = spm_vol(file1);

        v1  = spm_read_vols(spm_vol(file1));  % the actual volume
        v1(isnan(v1)) = 0;   % replace NaN if they are present

        v2  = spm_read_vols(spm_vol(file2));
        v2(isnan(v2)) = 0;   % replace NaN if they are present


        % the file dimensions need to be correct

        v1num = size(v1);
        v2num = size(v2);

        if sum(v1num - v2num) == 0
           % OK
        else
           fprintf('\n Warning: The two files do not have the same dimensions. Terminating.\n')
           return
        end


        dim1 = size(v1);   
        vol = zeros(dim1(1), dim1(2), dim1(3));

        % first, binarize the images, so we only deal with voxels common to both 

        v1bin = v1~=0;
        v2bin = v2~=0;

        v3bin = min(v1bin,v2bin); % will be 1 only if both v1 and v2 are non-zero

        % out of curiosity, what percentage of voxels are retained?
        v1val = sum(sum(sum(v1bin)));
        v2val = sum(sum(sum(v2bin)));

        v3val = sum(sum(sum(v3bin)));
        retained_vox(f,1) = (v1val / v3val) * 100;
        retained_vox(f,2) = (v2val / v3val) * 100;


            if op == 1
               vol = v1 + v2;
               
            elseif op == 2
               vol = v1 - v2;
               
            elseif op == 3
               vol = power(v1,pow(1)) .* power(v2,pow(2));
               
            elseif op == 4
               vol = power(v1,pow(1)) ./ power(v2,pow(2));
               
            elseif op == 5
                % the minimum statistic across the two volumes, at every
                % voxel
                vol = min(v1,v2);
                
            elseif op == 6
                % maximum statistic
                vol = max(v1,v2);
                
            elseif op == 7
               % the "null statistic" - take the value that is closer to zero
                vol = ((abs(v1)<abs(v2)) .* v1) + ((abs(v2)<abs(v1)) .* v2);
               
               % now, zero values where the sign of v1 and v2 does not match

                signchk = sign(v1) == sign(v2); % this will be zero if the signs do not agree
               
                % now, add the random eps vector to signchk - no, don't do this
                % signchk = signchk + (eps * sign(randn(v1num,1)));
                vol = vol .* signchk; % will turn voxels to zero if they don't have the same sign
                
            elseif op == 8
               % average of cond1 and cond2
               vol = (v1 + v2) / 2;
               
            end
            
            % replace NaN with zero (for voxels where 0/0)
            vol(isnan(vol)) = 0;
            
            % now, mask to the appropriate template - only if using
            % explicit mask!
            if use_mask == 1
               if maskch == 1 % inclusive
                    vol = vol .* v3bin;
               elseif maskch == 2 % exclusive
                    vol = vol;
               end
            end
            
            % mask (with mask image or itself) - vm is bounding box if no
            % mask is included
            vol = vol .* vm;
            
            
        % write the volume
        
     if isempty(findstr(file1,',1'))
        % OK 
     else
        % get rid of the suffix
        file1 = file1(1:findstr(file1,',1')-1);
     end
     
         v1n.('fname') = strcat(file1, op1,'.nii');    % will write to the "file1" directory
     % don't need to change pinfo; keep same s the input file

         % check for appropriate dtype
         dtype = dtype_check(vol(:));
         v1n.('dt') = dtype;
         
         v1n = spm_write_vol(v1n,vol);

    end % file loop

end % method 6

% % ----------------
% % for method 4
% if method == 4 || method == 5
% 
%    for f=1:numf1
%         if numf1 > 1
%            progressbar(f/numf1)
%         end
% 
%    file1 = files1(f,:);
%    v1n = spm_vol(file1);
%    v1  = spm_read_vols(spm_vol(file1));  % the actual volume
%    v1(isnan(v1)) = 0;   % replace NaN if they are present
%    sizev1 = size(v1);
%    
%    if f == 1 % get the template
%        vol = zeros(sizev1(1),sizev1(2),sizev1(3));
%    end
%    
%    if binit == 1
%       vol2 = v1~=0;  % binarize all non-zero voxels to 0
%    elseif binit == 2
%       % don't binarize
%       vol2 = v1;
%    end
%    
%    vol = vol + vol2; % add all the binarized images together consecutively
%     
%    end % file loop
%    
% 
%    % mask (with mask image or itself)
%    vol = vol .* vm;
%  
%      if method == 4  
%          % now, write the summed volume
%          v1n.('fname') = strcat(newdir, 'Group_BinSum.nii');    % will write to the "file1" directory
%          v1n.('dt') = [2 0]; % since binarized image
%          
%          v1n = spm_write_vol(v1n,vol);
% 
%      elseif method == 5
%          % now, rebinarize and write
%          vol = vol ~= 0;
% 
%          v1n.('fname') = strcat(newdir, 'Group_BinSumBin.nii');    % will write to the "file1" directory
%          v1n.('dt') = [2 0]; % since binarized image
%         
%         v1n = spm_write_vol(v1n,vol);
%      end
%  
% end % method == 4 || method == 5

% ----------------
% for method 3
if method == 3

   for f=1:numf1
       
        if numf1 > 1
           progressbar(f/numf1)
        end

       file1 = files1(f,:);
       thr = vec(f); % either manual or percentile

       v1n = spm_vol(file1);
       v1  = spm_read_vols(spm_vol(file1));  % the actual volume
       v1(isnan(v1)) = 0;   % replace NaN if they are present
       sizev1 = size(v1);

       if f == 1 % get the template
          vol = zeros(sizev1(1),sizev1(2),sizev1(3));
       end

       % pay attention to (1) manual vs percentile and (2) binarize vs preserve values       
       if thrit == 1
           vnew = v1; % copy of original values
       else
           vnew = ones(size(v1)); % simple mask of ones
       end
       
       if use_prc == 2
           % must turn the percentile into a reaal threshold
           v = v1(v1~=0); % all non-zero values (so percentiles are not biased); will be a vector
           
           thr = prctile_nist(v,thr); % doesn't require stats toolbox
       end
       
        % vnew handles whether we keep binarize or restore original values
        if UorD == 1
            vol2 = vnew .* (v1 > thr);              
        elseif UorD == 2
            vol2 = vnew .* (v1 >= thr); 
        elseif UorD == 3
            vol2 = vnew .* (v1 < thr);  
        elseif UorD == 4
            vol2 = vnew .* (v1 <= thr);  
        elseif UorD == 5
            vol2 = vnew .* (v1 == thr);  
        end      


       % add all the images together **cumulatively**
       vol = vol + vol2; 
       
       % optional rescaling step
       if f == numf1
           if rescale_it == 1
               vol = vol ./ numf1; % so that values go from 0 to 1 instead of 0 to N
               final_name = 'Group_ThrBinSumRS.nii';
           else
               % values remain from 0 to N
               final_name = 'Group_ThrBinSum.nii';
           end
       end
   
   end % file loop
   
   % mask (with mask image or itself)
   vol = vol .* vm;
    
     % check for appropriate dtype
     dtype = dtype_check(vol(:));
     v1n.('dt') = dtype;
     v1n.('fname') = strcat(newdir, final_name);    % will write to the "file1" directory
     
     spm_write_vol(v1n,vol); % volume written here
 
end % method == 3

% --------------
% for method 7: 

if method == 7

     for f=1:numf1
       
        if numf1 > 1
           progressbar(f/numf1)
        end

       file1 = files1(f,:);
       v1n = spm_vol(file1);
       v1  = spm_read_vols(spm_vol(file1));  % the actual volume
       v1(isnan(v1)) = 0;   % replace NaN if they are present


       vol = -1 * v1;  % flip sign

       % mask (with mask image or itself)
       vol = vol .* vm;

     % write the volume
     if isempty(findstr(file1,',1'))
        % OK 
     else
        % get rid of the suffix
        file1 = file1(1:findstr(file1,',1')-1);
     end     
     v1n.('fname') = strcat(file1, '_sflip','.nii');    % will write to the "file1" directory
     
     % check for appropriate dtype
     dtype = dtype_check(vol(:));
     v1n.('dt') = dtype;
     
     spm_write_vol(v1n,vol); % volume written here
   
     end % file loop
     
end % method == 7

% ------
% for method == 8
if method == 8
    for f=1:numf1
        
        if numf1 > 1
           progressbar(f/numf1)
        end
       
       file1 = files1(f,:);
       v1n = spm_vol(file1);
       v1  = spm_read_vols(spm_vol(file1));  % the actual volume
       v1(isnan(v1)) = 0;   % replace NaN if they are present


       vol = flipdim(v1,1);  % flip along the x-dimension

       % mask (with mask image or itself)
       vol = vol .* vm;

     % write the volume
     if isempty(findstr(file1,',1'))
        % OK 
     else
        % get rid of the suffix
        file1 = file1(1:findstr(file1,',1')-1);
     end
     
     v1n.('fname') = strcat(file1, '_xflip','.nii');    % will write to the "file1" directory

     % check for appropriate dtype
     dtype = dtype_check(vol(:));
     v1n.('dt') = dtype;
     
     spm_write_vol(v1n,vol); % volume written here
   
     end % file loop
end

% ------
% for method == 9
if method == 9
    for f=1:numf1
        
        if numf1 > 1
           progressbar(f/numf1)
        end

       file1 = files1(f,:);
       v1n = spm_vol(file1);
       v1  = spm_read_vols(spm_vol(file1));  % the actual volume
       v1(isnan(v1)) = 0; % replace NaN if they are present

       % mask (with mask image or itself)
       vol = v1 .* vm;

         % write the volume
         if isempty(findstr(file1,',1'))
            % OK 
         else
            % get rid of the suffix
            file1 = file1(1:findstr(file1,',1')-1);
         end     
         v1n.('fname') = strcat(file1, '_masked','.nii'); % will write to the "file1" directory
         
         % check for appropriate dtype
         dtype = dtype_check(vol(:));
         v1n.('dt') = dtype;
     
         spm_write_vol(v1n,vol); % volume written here
   
     end % f loop    
    
end

if numf1 > 1
    progressbar(1) % make sure it's closed
end

%% end
if method == 17
    if tot_files == 1
        fprintf('\n Finished. Output images written to the directory of the source image.\n\n')
    else
        fprintf('\n Finished. Output images written to the directories of the source images.\n\n')
    end
    
elseif method == 0 || method == 1 || method == 2 || method >= 7 
  
   if numf1 == 1
       fprintf('\n Finished. Output image written to the directory of the source image.\n\n');
   else
       fprintf(['\n Finished. ' num2str(numf1) ' output images written to the directories of the source images.\n\n']);
   end
   
   
   if method == 0
      % report out if any images have 0 voxels
      if max(zero_vox) > 0
          fprintf([' *Note*: images(s) [' num2str(zero_vox) '] have zero active voxels.\n\n'])
      else
          fprintf('\n');
      end
   end
   
elseif method == 3
    if same_dir == 1
        fprintf('\n Finished. Output image written to the directory of the input files. \n\n')
    else
        fprintf('\n Finished. Output image written to the directory specified. \n\n');
    end
% elseif method == 4 || method == 5
%     if numf1 == 1
%         fprintf('\n Finished. Summed binarized image written to the directory specified. \n\n');
%     else
%         fprintf(['\n Finished. ' num2str(numf1) ' summed binarized images written to the directory specified. \n\n']);
%     end
    
elseif method == 6
    if numf1 == 1
        fprintf('\n Finished. Output image written to the directories of the *Group 1* file.\n\n');
    else
        fprintf(['\n Finished. ' num2str(numf1) ' output image(s) written to the directories of the *Group 1* file(s).\n\n']);
    end
   
end

% remove copies of SPM functions from the path to avoid errors
remove_subfolder('imcalc.m','spm_fxns')

% finally, delete the resliced mask(s)
cd(maskdir)
delete('rs_*')
cd(curdir)