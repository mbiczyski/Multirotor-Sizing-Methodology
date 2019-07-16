%% ------------------------------------------------------------------------
% Multirotor Sizing Methodology
% with Flight Time Estimation
%
% M. Biczyski, R. Sehab, G.Krebs, J.F. Whidborne, P. Luk
%
% load_propList.m - implements function load_propList() that reads basic
% propeller information from the database; requires file 'WEBLIST.xlsx' 
% in 'PERFILES_WEB' folder and file 'PROP-DATA-FILE_201812.xlsx' from 
% APC Performance Database Summary
%% ------------------------------------------------------------------------

% output = name, file, diameter (in), pitch (in), mass (g)
function output = load_propList()
    [~, ~, everything] = xlsread('PERFILES_WEB/WEBLIST.xlsx'); % read from the data file

    fileList = [{} {}];
    for ii = 6:size(everything,1) % skip first 5 lines and save propeller names and respective filenames
        fileList(ii-5,:) = [everything(ii,3) everything(ii,1)];
    end

    clear everything; % clear unused data

    [~, ~, everything] = xlsread('PROP-DATA-FILE_201812.xlsx', 2); % read from the data file

    for ii = 1:size(fileList,1) % look for each propeller specified by name and save its parameters
        for jj = 2:size(everything,1)
            if startsWith(everything{jj, 1}, fileList{ii,1})
                fileList(ii, 3:5) = [str2num(everything{jj, 11})  % diameter
                                     everything(jj, 12)  % pitch
                                     everything(jj, 16)]; % mass
                break;
            end
        end
    end

    fileList(cellfun('isempty', fileList(:,3)), :) = []; %remove entries without parameters
    
    % RPM limits
    % E   - 145000/D
    % F   - 145000/D
    % MR  - 105000/D
    % SF  - 65000/D
    % E-3 - 270000/D
    % E-4 - 270000/D
    
    for ii = 1:size(fileList, 1) % depending on propeller series and diameter add speed limit parameter
        temp1 = fileList{ii,1};
        temp2 = temp1(isstrprop(temp1,'alpha'));
        switch temp2(2:end)
            case 'PN'
                SpeedLimitSet = 190000;
            case {'W' 'N' 'P'}
                SpeedLimitSet = 270000;
            case {'E' 'F' 'WE' 'EPN' 'ECD'}
                SpeedLimitSet = 145000;
            case 'MR'
                SpeedLimitSet = 105000;
            case 'SF'
                SpeedLimitSet = 65000;
            case {'E-3' 'E-4' 'C'}
                SpeedLimitSet = 270000;
            otherwise
                SpeedLimitSet = 225000;
        end
        fileList{ii,6} = SpeedLimitSet/fileList{ii,3};
    end

    output = fileList;
end