function [same direc] = dir_check(files)

% check the directories in "files" and return a 
%   1 - if they are all the same
%   0 - otherwise

same = [];
direc  = [];

% how many files?
nfiles = size(files,1);

if nfiles >= 2
    file1 = strtrim(files(1,:));
    path1 = fileparts(file1);
    
    for n = 2:nfiles
        % we do this using simple logic; just test adjacent pairs of files
        file2 = strtrim(files(n,:));
        path2 = fileparts(file2);
        
        % are they the same?
        if strcmp(path1,path2)
            % keep going ...
            same = 1;
            direc = path1;
        else
            % mismatch
            same = 0;
            break
        end

    end    
    
else
    % pointless to continue
    return
end

