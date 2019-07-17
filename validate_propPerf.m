%% ------------------------------------------------------------------------
% Multirotor Sizing Methodology
% with Flight Time Estimation
%
% M. Biczyski, R. Sehab, G.Krebs, J.F. Whidborne, P. Luk
%
% validate_propPerf.m - top-level script that calculates and compares 
% simulated and measured performance of APC propellers; requires 
% 'PROP-DATA-EXPERIMENTAL.xlsx' file based on:
% J.B. Brandt, R.W. Deters, G.K. Ananda, and M.S. Selig (21.05.2019), 
% UIUC Propeller Database, University of Illinois at Urbana-Champaign, 
% retrieved from http://m-selig.ae.illinois.edu/props/propDB.html.
%% ------------------------------------------------------------------------

close all; clear; clc;

%% load experimental propellers list
[num, text, everything] = xlsread('PROP-DATA-EXPERIMENTAL.xlsx'); % reads experimental data

expList = {text{2,1}, 3, 0}; % prepare first entry
for ii = 3:size(text,1) % fill the list with names and data ranges of propellers
    if ~isempty(text{ii,1})
        expList(end,3) = {ii-1};
        expList(end+1,:) = {text{ii,1}, ii+1, 0};
    end
end
expList(end,3) = {size(everything,1)};

%% load experimental propellers performance
expPerf = {};
for ii = 1:size(expList, 1) % read thrust and power coefficients varying with speed
    temp_String = expList{ii,1}(end-5:end-4);
    propGroup = temp_String(isstrprop(temp_String,'alpha')); % read propeller series or set it as Sport
    if strcmp(propGroup,'x') == true || isempty(propGroup)
        propGroup = 'S';
    end
    
    expPerf(ii,:) = {[everything{expList{ii,2}:expList{ii,3},1}]',...
                     [everything{expList{ii,2}:expList{ii,3},2}]',...
                     [everything{expList{ii,2}:expList{ii,3},3}]',...
                     propGroup};
end

%% load simulated propellers performance
simPerf = {};
for ii = 1:size(expList, 1) % read propeller performance provided by the manufacturer
    temp_simPerf = load_propPerf(expList{ii,1}, false);
    simPerf(ii,:) = {temp_simPerf(:,1),...
                     temp_simPerf(:,6),...
                     temp_simPerf(:,5)};
end

%% calculate average errors
errPerf = [];
for ii = 1:size(expList, 1) % calculate mean absolute and relative errors, standard deviation and mean values for C_T and C_P
    errPerf(ii,1) = mean(expPerf{ii,2} - interp1(simPerf{ii,1}, simPerf{ii,2}, expPerf{ii,1}), 'omitnan'); % Mean error C_T
    errPerf(ii,2) = mean((expPerf{ii,2} - interp1(simPerf{ii,1}, simPerf{ii,2}, expPerf{ii,1}))./expPerf{ii,2}, 'omitnan').*100; % Mean relative error C_T
    errPerf(ii,3) = std(expPerf{ii,2} - interp1(simPerf{ii,1}, simPerf{ii,2}, expPerf{ii,1}), 'omitnan'); % Standard deviation C_T
    errPerf(ii,4) = mean(expPerf{ii,3} - interp1(simPerf{ii,1}, simPerf{ii,3}, expPerf{ii,1}), 'omitnan'); % Mean error C_P
    errPerf(ii,5) = mean((expPerf{ii,3} - interp1(simPerf{ii,1}, simPerf{ii,3}, expPerf{ii,1}))./expPerf{ii,3}, 'omitnan').*100; % Mean relative error C_P
    errPerf(ii,6) = std(expPerf{ii,3} - interp1(simPerf{ii,1}, simPerf{ii,3}, expPerf{ii,1}), 'omitnan'); % Standard deviation C_P
    avgPerf(ii,1) = mean(expPerf{ii,2});
    avgPerf(ii,2) = mean(simPerf{ii,2});
    avgPerf(ii,3) = mean(expPerf{ii,3});
    avgPerf(ii,4) = mean(simPerf{ii,3});
end

% thrust coefficient is overpredicted
% torque is underpredicted
averageError = [mean(errPerf(:,1)), mean(errPerf(:,2)), mean(errPerf(:,3)), mean(errPerf(:,4)), mean(errPerf(:,5)), mean(errPerf(:,6))]; % average the errors across all samples

%% plot results
figure;
subplot(2,1,1);
bar(errPerf(:,2), 1);
title('Mean relative thrust coefficient (C_T) error');
ylabel('Error (%)');
xlabel('Propeller sample number');
subplot(2,1,2);
bar(errPerf(:,5), 1);
title('Mean relative power coefficient (C_P) error');
ylabel('Error (%)');
xlabel('Propeller sample number');

figure;
cat = categorical({'Carbon', 'Sport (IC)', 'Electrical', 'Slow Flyer'});
subplot(1,2,1);
bar(cat, [mean(avgPerf(strcmp(expPerf(:,4),'C'),1)) mean(avgPerf(strcmp(expPerf(:,4),'S'),1)) mean(avgPerf(strcmp(expPerf(:,4),'E'),1)) mean(avgPerf(strcmp(expPerf(:,4),'SF'),1))]);
ylabel('Thrust coefficient');
subplot(1,2,2);
bar(cat, [mean(avgPerf(strcmp(expPerf(:,4),'C'),3)) mean(avgPerf(strcmp(expPerf(:,4),'S'),3)) mean(avgPerf(strcmp(expPerf(:,4),'E'),3)) mean(avgPerf(strcmp(expPerf(:,4),'SF'),3))]);
ylabel('Power coefficient');


