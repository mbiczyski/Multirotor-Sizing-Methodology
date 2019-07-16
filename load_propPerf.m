%% ------------------------------------------------------------------------
% Multirotor Sizing Methodology
% with Flight Time Estimation
%
% M. Biczyski, R. Sehab, G.Krebs, J.F. Whidborne, P. Luk
%
% load_propPerf.m - implements function load_propPerf() that reads static 
% propeller performance from the database; requires file 
% 'PER2_STATIC-2.DAT' from APC Performance Database Summary in
% 'PERFILES_WEB' folder
%% ------------------------------------------------------------------------

% output = RPM, Thrust (g), Power (W), Torque (Nm), Cp, Ct
function output = load_propPerf(name, bPlot)
    fileID = fopen('PERFILES_WEB/PER2_STATIC-2.DAT', 'rt'); % access the data file
    
    temp1 = fgetl(fileID);
    while ~strcmp(sscanf(temp1, '%s'), name(6:end)) % look for the propeller specified by name
        temp1 = fgetl(fileID);
    end
    
    for ii = 1:4 % skip blank space & headings
        fgetl(fileID);
    end
    
    currentLine = sscanf(fgetl(fileID), '%f %f %f %f %f %f'); % read first line of data
    staticData = [];
    while ~isempty(currentLine) % read next lines of data until blank space or end of file
        staticData(end+1, :) = currentLine;
        temp2 = fgetl(fileID);
        if temp2 == -1
            temp2 = '';
        end
        currentLine = sscanf(temp2, '%f %f %f %f %f %f');
    end
    fclose(fileID); % close data file
    
    staticData(:,2) = staticData(:,2)*4.4482216/9.81*1000; % convert thrust values to gram-force
    staticData(:,3) = staticData(:,3)*746; % convert power values to Watt
    staticData(:,4) = staticData(:,4)*0.112985; % convert torque values to Newton-meter
    
    if bPlot % (optional) plot static propeller characteristics
        figure;
        subplot(2,1,1);
        yyaxis left;
        plot(staticData(:,1), staticData(:,2));
        xlabel('RPM');
        ylabel('Thrust (g)');
        ylim([0 1.1*max(staticData(:,2))]);
        title([name(6:end-4) ' - Static Thrust / Torque']);
        yyaxis right;
        plot(staticData(:,1), staticData(:,4));
        ylabel('Torque (Nm)');
        ylim([0 1.1*max(staticData(:,4))]);
        grid;
        
        subplot(2,1,2);
        yyaxis left;
        plot(staticData(:,1), staticData(:,6));
        xlabel('RPM');
        ylabel('Thrust Coefficient');
        ylim([0 1.1*max(staticData(:,6))]);
        title([name(6:end-4) ' - Thrust / Power Coefficients']);
        yyaxis right;
        plot(staticData(:,1), staticData(:,5));
        ylabel('Power Coefficient');
        ylim([0 1.1*max(staticData(:,5))]);
        grid;
    end
    
    output = staticData;
end