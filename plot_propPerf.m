%% ------------------------------------------------------------------------
% Multirotor Sizing Methodology
% with Flight Time Estimation
%
% M. Biczyski, R. Sehab, G.Krebs, J.F. Whidborne, P. Luk
%
% plot_propPerf.m - plots data related to propeller's performance; requires
% main.m or example_DJI_Phantom_4_v2.m to be run before
%% ------------------------------------------------------------------------

plot_scaler = 1.1;

%% Thrust-speed charcateristic
figure; hold on;
temp_plot_points = [];
temp_RPMRange = 1000:100:max([propList_considered{:,6}]);
for ii = 1:consideredNo
    plot(temp_RPMRange, interp1(propPerf{ii}(:,1), propPerf{ii}(:,2), temp_RPMRange));
end
colorOrder = get(gca, 'ColorOrder');
for ii = 1:consideredNo
    scatter(operatingPoints{ii,1}(1), operatingPoints{ii,1}(2), 30, colorOrder(mod(ii-1,7)+1,:), 'filled');
    scatter(operatingPoints{ii,2}(1), operatingPoints{ii,2}(2), 30, colorOrder(mod(ii-1,7)+1,:), 'filled');
    scatter(operatingPoints{ii,3}(1), operatingPoints{ii,3}(2), 30, colorOrder(mod(ii-1,7)+1,:), 'filled');
    temp_plot_points(end+1) = operatingPoints{ii,2}(2);
end
hold off; grid;
legend (propList_considered{:,1}, 'Location', 'NorthWest');
ylim([0 plot_scaler*max(temp_plot_points)]);
ylabel('Thrust (g)');
xlabel('Speed (RPM)');

%% Torque-speed characteristic
figure; hold on;
temp_plot_points = [];
for ii = 1:consideredNo
    plot(temp_RPMRange, interp1(propPerf{ii}(:,1), propPerf{ii}(:,4), temp_RPMRange));
end
for ii = 1:consideredNo
    scatter(operatingPoints{ii,1}(1), operatingPoints{ii,1}(3), 30, colorOrder(mod(ii-1,7)+1,:), 'filled');
    scatter(operatingPoints{ii,2}(1), operatingPoints{ii,2}(3), 30, colorOrder(mod(ii-1,7)+1,:), 'filled');
    scatter(operatingPoints{ii,3}(1), operatingPoints{ii,3}(3), 30, colorOrder(mod(ii-1,7)+1,:), 'filled');
    temp_plot_points(end+1) = operatingPoints{ii,2}(3);
end
hold off; grid;
legend (propList_considered{:,1}, 'Location', 'NorthWest');
ylim([0 plot_scaler*max(temp_plot_points)]);
ylabel('Torque (Nm)');
xlabel('Speed (RPM)');

%% Power-speed characteristic
figure; hold on;
temp_plot_points = [];
for ii = 1:consideredNo
    plot(temp_RPMRange, interp1(propPerf{ii}(:,1), propPerf{ii}(:,3), temp_RPMRange));
end
for ii = 1:consideredNo
    scatter(operatingPoints{ii,1}(1), operatingPoints{ii,1}(4), 30, colorOrder(mod(ii-1,7)+1,:), 'filled');
    scatter(operatingPoints{ii,2}(1), operatingPoints{ii,2}(4), 30, colorOrder(mod(ii-1,7)+1,:), 'filled');
    scatter(operatingPoints{ii,3}(1), operatingPoints{ii,3}(4), 30, colorOrder(mod(ii-1,7)+1,:), 'filled');
    temp_plot_points(end+1) = operatingPoints{ii,2}(4);
end
hold off; grid;
legend (propList_considered{:,1}, 'Location', 'NorthWest');
ylim([0 plot_scaler*max(temp_plot_points)]);
ylabel('Power (W)');
xlabel('Speed (RPM)');

%% Motor load versus ideal motor torque
figure; hold on;
plot(temp_RPMRange, interp1(propPerf{temp_propChosen_pos}(:,1), propPerf{temp_propChosen_pos}(:,4), temp_RPMRange));
plot(temp_RPMRange(15:end), motorSpecification{5}./(temp_RPMRange(15:end)*2*pi/60));
scatter(operatingPoints{temp_propChosen_pos,1}(1), operatingPoints{temp_propChosen_pos,1}(3), 50, colorOrder(1,:), 'filled');
scatter(operatingPoints{temp_propChosen_pos,2}(1), operatingPoints{temp_propChosen_pos,2}(3), 50, colorOrder(1,:), 'filled');
scatter(operatingPoints{temp_propChosen_pos,3}(1), operatingPoints{temp_propChosen_pos,3}(3), 50, colorOrder(1,:), 'filled');
hold off; grid;
legend ({[propSpecification{1} ' propeller'] ['Constant power - ' num2str(round(motorSpecification{5})) ' W']}, 'Location', 'NorthWest');
ylabel('Torque (Nm)');
xlabel('Speed (RPM)');
title('Torque characteristics of the propulsion system');

%% Battery discharge at hover
figure;
yyaxis left;
plot(time_hover*60, voltage_hover); hold on;
plot(time_hover*60, current_hover);
ylabel('--- Voltage (V) / - - Current (A)'); ylim([0 1.1*max(current_max)]);
yyaxis right;
plot(time_hover*60, capacity_hover*1000);
ylabel('Available capacity (mAh)'); ylim([0 1.1*max(capacity_hover*1000)]);
xlabel('Time (min)'); xlim([0 time_hover(end)*60]); grid;
title([num2str(baterrySpecification(1)) 'S ' num2str(ceil(baterrySpecification(2))) 'C ' num2str(ceil(baterrySpecification(3)/100)*100) 'mAh - energy usage at hover']);

%% Battery discharge at WOT
figure;
yyaxis left;
plot(time_max*60, voltage_max); hold on;
plot(time_max*60, current_max);
ylabel('--- Voltage (V) / - - Current (A)'); ylim([0 1.1*max(current_max)]);
yyaxis right;
plot(time_max*60, capacity_max*1000);
ylabel('Available capacity (mAh)'); ylim([0 1.1*max(capacity_max*1000)]);
xlabel('Time (min)'); xlim([0 time_max(end)*60]); grid;
title([num2str(baterrySpecification(1)) 'S ' num2str(ceil(baterrySpecification(2))) 'C ' num2str(ceil(baterrySpecification(3)/100)*100) 'mAh - energy usage at 100% throttle']);