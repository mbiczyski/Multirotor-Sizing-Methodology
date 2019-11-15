%% ------------------------------------------------------------------------
% Multirotor Sizing Methodology
% with Flight Time Estimation
%
% M. Biczyski, R. Sehab, G.Krebs, J.F. Whidborne, P. Luk
%
% plot_motorPerf.m - plots data related to motors's performance; requires
% main.m or example_DJI_Phantom_4_v2.m to be run before
%
% uses function addaxis version 1.1.0.0 by Harry Lee
%% ------------------------------------------------------------------------

conn = sqlite('DCbase.dcd'); % connect to the database
motors = fetch(conn, ['SELECT * FROM Motors WHERE myid = ' num2str(motorChosen{1})]); % run SQL query to obtain motor performance data based on specified ID
mdata = fetch(conn, 'SELECT * FROM MData'); % run SQL query to obtain full 'MData' table
esc = fetch(conn, 'SELECT * FROM ESC'); % run SQL query to obtain full 'ESC' table
close(conn); % close connection to database

esc_id = motors{6}; % obtain ID of ESC recommended for the motor
performance = mdata([mdata{:,3}] == motorChosen{1},:); % look for test data for the motor with given ID

current_maximum = motors{13}; % obtain maximum motor current
if current_maximum < 1 % if no maximum current data is available, obtain it from ESC's maximum current or guess based on voltage
    if esc_id > 0 && esc{esc_id,5} > 0
        current_maximum = esc{esc_id,5};
    else
        current_maximum = 2*BattCellNo*BattCellVoltage;
    end
end
current = 1:current_maximum; % prepare current domain
[~, noLoad_pos] = min(abs([performance{[performance{:,7}]==1,4}]-BattCellNo*BattCellVoltage)); % find no-load test data close to target voltage
Rm = 2*motors{18}; % read motor windings resitance
kV = motors{17}; % read motor KV rating

coppLoss = Rm*current.^2; %calclulate copper losses based on windings resistance
ironLoss = BattCellNo*BattCellVoltage*performance{noLoad_pos,6}; % calculate iron losses based on no-load current
loss = coppLoss+ironLoss; % calculate total losses
speed = (BattCellNo*BattCellVoltage-coppLoss./current)*kV; % calculate maximum speed
powerMax = BattCellNo*BattCellVoltage*current; % calculate electrical power
power = powerMax - loss; % calculate motor power
efficiency = ((powerMax-loss)./powerMax); % calculate efficiency

%% Motor performance characteristics
figure;
plot(current, speed); 
ylim([0 ceil(max(speed)*1.1/1000)*1000]);
ax = gca;
ax.YColor = [0, 0.4470, 0.7410];

scaler = ceil(max(power)*1.1/100)*100;
addaxis(current, power, [0 scaler]);
addaxisplot(current, powerMax, 2, 'k-.');

addaxis(current(2:end), efficiency(2:end)*100, [40 100], 'Color', [0.4660, 0.6740, 0.1880]);

addaxisplot(motorChosen{7}, operatingPoints{temp_propChosen_pos,2}(3)*operatingPoints{temp_propChosen_pos,2}(1)*RPM2RAD, 2, '.', 'MarkerSize', 20);
addaxisplot(motorChosen{10}, operatingPoints{temp_propChosen_pos,1}(3)*SafetyFactor*operatingPoints{temp_propChosen_pos,1}(1)*RPM2RAD, 2, '.', 'MarkerSize', 20);
addaxisplot(motorChosen{7}, motorChosen{9}, 3, '.', 'MarkerSize', 20);
addaxisplot(motorChosen{10}, motorChosen{12}, 3, '.', 'MarkerSize', 20);

grid on; grid minor;
xlim([0 ceil(current_maximum/10)*10]);
xlabel('Current (A)');
title([motorChosen{2} ' (' num2str(round(motorSpecification{2}/10)*10) ' kV)']);
legend('Speed (RPM)', 'Mech. power (W)', 'El. Power (W)', 'Efficiency (%)', 'Location', 'southeast');