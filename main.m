%% ------------------------------------------------------------------------
% Multirotor Sizing Methodology
% with Flight Time Estimation
%
% M. Biczyski, R. Sehab, G.Krebs, J.F. Whidborne, P. Luk
%
% main.m - base file of the script realising input preparation, propeller 
% and motor selection, battery simulation and results generation
%% ------------------------------------------------------------------------

close all; clear; clc;
format compact; format shortG;

RPM2RAD = 2*pi/60;
RAD2RPM = 60/2/pi;

%% User parameters
RotorNo = 4; % Number of rotors
OptimisationGoal = 'hover'; % selection criteria
                    % hover - best specific thrust (g/W) at hover
                    % max - best specific thrust (g/W) at 100% throttle
                    % utilisation - maximum usable power range of propeller
ThrustWeightRatio = 3; % estimator for maximum paerformance
                    % 2 - minimum
                    % 3 - payload transport
                    % 4 - survaillence
                    % 5+ - aerobatics / hi-speed video
                    % 7+ - racing
PropDiameter_Min = 12; % inch, min. propeller diameter
PropDiameter_Max = 12; % inch, max. propeller diameter
SafetyFactor = 1.1; % [1-2], arbitrary safety parameter
AcceptedTypes = {'MR' 'E-3' 'E-4'}; % preferred propeller series
                  % E	Electric
                  % F	Folding Blade (electric only)
                  % MR	Multi-Rotor (electric only)
                  % SF	Slow Fly (electric only)
                  % E-3	3-Blade
                  % E-4	4-Blade

BattCellNo = 4; %S 1P, battery cell count
BattCellVoltage = 3.7; % V per cell, battery cell voltage
BattCapacity = 5200; % mAh, battery capacity
BattPeukertConstant = 1.3; % for LiPo, Peukert's constant for selected type of battery
BattVoltageSagConstant = 0.5/0.8*BattCellNo; % 0.5V decrease per cell in resting volatage for 80% DoD
BattHourRating = 1; % h

%% Mass data [g]
mass_Frame = 543; % Lumenier QAV500 V2 with 540mm arms
mass_FC = 21; % Vector Flight Controller + OSD
mass_FC_GPS = 13;
mass_FC_CurrentSensor = 15;
mass_Receiver = 2; % FrSky R-XSR 2.4GHz 16CH ACCST
mass_Motor_Est = 184; % FXC4006-13 740kv - 92 g
mass_ESC_Est = 14; % Lumenier 35A BLHeli_S ESC OPTO - 7 g
mass_Propeller_Est = 36; % HQProp 12x4.5 Props - 18g
% mass_Payload = 977;
mass_Payload = 0;
mass_Battery = 473; % Lumenier 5200mAh 4s 35c
mass_Other_Est = 20; % cabling, straps, standoffs, etc.

mass_NoDrive_Est = mass_Frame + mass_FC + mass_FC_GPS + mass_FC_CurrentSensor + mass_Receiver + mass_Payload + mass_Other_Est;
mass_Total_Est = mass_NoDrive_Est + RotorNo*(mass_Motor_Est + mass_ESC_Est + mass_Propeller_Est) + mass_Battery;

%% Filter propeller set
% propList = name, file, diameter (in), pitch (in), mass (g), speedLimit (RPM)
propList = load_propList(); % loading propeller set

propList_considered_bySize = propList(cellfun(@(x) x>=PropDiameter_Min && x<=PropDiameter_Max, propList(:,3)), :); % filtering propellers based on diameter
propList_considered_byMass = propList(cellfun(@(x) x<=mass_Propeller_Est, propList(:,5)), :); % filtering propellers based on mass
propList_considered_byType = propList(cellfun(@(x) endsWith(x, AcceptedTypes), propList(:,1)), :); % filtering propellers based on series

temp_propList_considered_names = intersect(intersect(propList_considered_bySize(:,1), propList_considered_byMass(:,1)), propList_considered_byType(:,1));

propList_considered = propList(ismember(propList(:,1), temp_propList_considered_names),:);  % intersections of filtered sets
consideredNo = size(propList_considered,1); % size of filtered set

if consideredNo < 1
    error('ERROR! No matching propeller found!');
end

%% Load propeller performance
% propPerf = RPM, Thrust (g), Power (W), Torque (Nm), Cp, Ct
% operatingPoints = {hover, max, limit}
%                   [speed, thrust, torque, power]

propPerf = {};
for ii = 1:consideredNo
    % TRUE/FALSE for plot
    propPerf(ii) = {load_propPerf(propList_considered{ii,2}, false)}; % loading propeller static performance data
end

%% Calculate operating points
thrustHover_Est = mass_Total_Est/RotorNo; % calcculate thrust required for hover
thrustMax_Est = thrustHover_Est*ThrustWeightRatio; % calculate estimated thrust at WOT
for ii = 1:consideredNo
    speedHover = interp1(propPerf{ii}(2:end,2), propPerf{ii}(2:end,1), thrustHover_Est); % obtaining propeller speed at hover from required thrust for hover
    speedMax = interp1(propPerf{ii}(2:end,2), propPerf{ii}(2:end,1), thrustMax_Est); % obtaining propeller speed at WOT from estimated thrust at WOT 
    speedLimit = propList_considered{ii,6}; % obtaining propeller's limiting speed specified by the manufacturer
    operatingPoints(ii,1) = {[speedHover thrustHover_Est interp1(propPerf{ii}(:,1), propPerf{ii}(:,4), speedHover) interp1(propPerf{ii}(:,1), propPerf{ii}(:,3), speedHover)]}; % obtaining hover operating point
    operatingPoints(ii,2) = {[speedMax thrustMax_Est interp1(propPerf{ii}(:,1), propPerf{ii}(:,4), speedMax) interp1(propPerf{ii}(:,1), propPerf{ii}(:,3), speedMax)]}; % obtaining WOT operating point
    operatingPoints(ii,3) = {[speedLimit interp1(propPerf{ii}(:,1), propPerf{ii}(:,2), speedLimit) interp1(propPerf{ii}(:,1), propPerf{ii}(:,4), speedLimit)...
        interp1(propPerf{ii}(:,1), propPerf{ii}(:,3), speedLimit)]}; % obtaining speed limit operating point
end

%% Select propeller
for ii = 1:consideredNo
    switch OptimisationGoal % selection of approperiate criteria based on user's choice
        case 'hover'
            selectionCriterion(ii,1) = (operatingPoints{ii,1}(1)*2*pi/60)*operatingPoints{ii,1}(3);
            selectionCriterion(ii,2) = operatingPoints{ii,1}(4);  % power at hover
        case 'max'
            selectionCriterion(ii,1) = (operatingPoints{ii,2}(1)*2*pi/60)*operatingPoints{ii,2}(3);
            selectionCriterion(ii,2) = operatingPoints{ii,2}(4); % power at WOT
        case 'utilisation'
            selectionCriterion(ii,1) = (operatingPoints{ii,3}(1)*2*pi/60)*operatingPoints{ii,3}(3) - (operatingPoints{ii,2}(1)*2*pi/60)*operatingPoints{ii,2}(3);
            selectionCriterion(ii,2) = operatingPoints{ii,3}(4) - operatingPoints{ii,2}(4); % best usage of propeller's speed range
        otherwise
            error('ERROR! Wrong optimisation criteria!');
    end
end

methodError(:,1) = abs(selectionCriterion(:,2) - selectionCriterion(:,1)); % absolute error between power and the product of speed and torque due to interpolation
methodError(:,2) = abs(selectionCriterion(:,2) - selectionCriterion(:,1))./abs(selectionCriterion(:,2)); % relative interpolation error
for ii = 1:consideredNo
    if operatingPoints{ii,3}(4) < operatingPoints{ii,2}(4) || isnan(operatingPoints{ii,2}(1))
        selectionCriterion(ii,:) = inf; % rejecting propellers with numerical errors and with WOT speed over limit speed
    end
end

[~, temp_propChosen_pos] = min(mean(selectionCriterion,2)); % selection of best propeller for the application

if selectionCriterion(temp_propChosen_pos,2) == inf
    error('ERROR! No matching propeller found!');
end

%% Load & filter motor data
% motorList = ID, name, ILimit (A), mass (g), kV, Rm (Ohm),
            % op_Imax (A), op_powerMaxEl(W), op_effMax (%), op_IHover (A), op_powerHoverEl (W), op_effHover (%)
            
motorList = load_motorList(BattCellNo*BattCellVoltage, operatingPoints{temp_propChosen_pos,2}(1), operatingPoints{temp_propChosen_pos,2}(3),...
    operatingPoints{temp_propChosen_pos,1}(1), operatingPoints{temp_propChosen_pos,1}(3)*SafetyFactor,...
    mass_Motor_Est); % loading motor set with operating points

if size(motorList,1) < 1
    error('ERROR! No matching motor found!');
end

%% Select motor
switch OptimisationGoal % selection of approperiate criteria based on user's choice
    case 'hover'
        [~, temp_motorChosen_pos] = min([motorList{:,11}]); % power at hover
    case 'max'
        [~, temp_motorChosen_pos] = min([motorList{:,8}]); % power at WOT
    case 'utilisation'
        [~, temp_motorChosen_pos] = min(abs([motorList{:,3}]-[motorList{:,7}]));  % best usage of motor's current range
    otherwise
        error('ERROR! Wrong optimisation criteria!');
end

motorChosen = motorList(temp_motorChosen_pos,:); % selection of best motor for the application

%% Determine drive specification
% propSpecification = name, diameter (in), pitch (in)
% motorSpecification = name, speedMax (RPM), torqueMax (Nm), powerMax (W), powerMaxEl (W), EfficiencyMax(%), voltageNominal (V)
% escSpecification = currentMax (A)
% baterrySpecification = NoCells, C-rating, minCapacity (mAh)

propSpecification = propList_considered(temp_propChosen_pos,[1 3:4]);
motorSpecification = {motorChosen{2}, motorChosen{5}, operatingPoints{temp_propChosen_pos,2}(1), operatingPoints{temp_propChosen_pos,2}(3)*SafetyFactor,...
    operatingPoints{temp_propChosen_pos,2}(4)*SafetyFactor, motorChosen{8}, motorChosen{9}, BattCellNo*BattCellVoltage};
escSpecification = motorChosen{7};

BattCapacity_Ah = BattCapacity/1000;
minBattRating = escSpecification*RotorNo*SafetyFactor/BattCapacity_Ah; % calculate min. battery C-rating required to supply enough current to motors
baterrySpecification = [BattCellNo, minBattRating, BattCapacity];

mass_Propeller = propList_considered{temp_propChosen_pos,5};
mass_Motor = motorChosen{4};
mass_ESC = mass_ESC_Est;
mass_Total = mass_NoDrive_Est + RotorNo*(mass_Motor + mass_ESC + mass_Propeller) + mass_Battery; % recalculate total mass of multirotor using real component weights

%% Calculate initioal battery state
voltage_hover(1) = (BattCellNo*(BattCellVoltage+0.5)); % 4.2 V per cell times number of cells
voltage_max(1) = (BattCellNo*(BattCellVoltage+0.5));
current_hover(1) = motorChosen{11}*RotorNo/voltage_hover(1); % calculate total current at hover
current_max(1) = motorChosen{8}*RotorNo/voltage_max(1); % calculate total current at WOT
capacity_hover(1) = (current_hover(1)^(1-BattPeukertConstant))*(BattHourRating^(1-BattPeukertConstant))*(BattCapacity_Ah^BattPeukertConstant); % from modified Peukert's equation calculate available capacity at hover
capacity_max(1) = (current_max(1)^(1-BattPeukertConstant))*(BattHourRating^(1-BattPeukertConstant))*(BattCapacity_Ah^BattPeukertConstant); % from modified Peukert's equation calculate available capacity at WOT

%% Calculate next flight iterations
timeStep = 1/60/60; % set timestep as 1 s
ii = 1;
while voltage_hover(ii) > BattCellVoltage*BattCellNo && ii*timeStep < 2
    voltage_hover(ii+1) = voltage_hover(1) - (BattVoltageSagConstant/capacity_hover(1))*(capacity_hover(1) - capacity_hover(ii)); % calculate instantaneus voltage including voltage sag
    current_hover(ii+1) = motorChosen{11}*RotorNo/voltage_hover(ii+1); % calculate instantaneus current based on required power for hover
    capacity_hover(ii+1) = (current_hover(ii+1)^(1-BattPeukertConstant))*(BattHourRating^(1-BattPeukertConstant))*(BattCapacity_Ah^BattPeukertConstant) - sum(current_hover(2:end)*timeStep); % calculate remaining available capacity according to Paeukert's effect
    ii = ii+1;
end
time_hover = (0:ii-1)*timeStep; % calculate time spent in hover

ii = 1;
while voltage_max(ii) > BattCellVoltage*BattCellNo && ii*timeStep < 2
    voltage_max(ii+1) = voltage_max(1) - (BattVoltageSagConstant/capacity_max(1))*(capacity_max(1) - capacity_max(ii)); % calculate instantaneus voltage including voltage sag
    current_max(ii+1) = motorChosen{8}*RotorNo/voltage_max(ii+1); % calculate instantaneus current based on estimated power at WOT
    capacity_max(ii+1) = (current_max(ii+1)^(1-BattPeukertConstant))*(BattHourRating^(1-BattPeukertConstant))*(BattCapacity_Ah^BattPeukertConstant) - sum(current_max(2:end)*timeStep); % calculate remaining available capacity according to Paeukert's effect
    ii = ii+1;
end
time_max = (0:ii-1)*timeStep; % calculate time spent at WOT

%% Display results and plot characteristics
disp(['For a ' num2str(RotorNo) '-rotor drone with estimated AUM of ' num2str(round(mass_Total_Est)) ' g (calculated TOM of ' num2str(round(mass_Total)) ' g):']);

switch OptimisationGoal
    case 'hover'
        textOptimisation = ['the highest specific thrust of ' num2str(round(operatingPoints{temp_propChosen_pos,1}(2)/motorChosen{11}*100)/100)  ' gf/W per motor at hover.'];
    case 'max'
        textOptimisation = ['the highest specific thrust of ' num2str(round(operatingPoints{temp_propChosen_pos,2}(2)/motorChosen{8}*100)/100)  ' gf/W per motor at WOT.'];
    case 'utilisation'
        textOptimisation = 'maximum usable power range of propeller';
    otherwise
        error('ERROR! Wrong optimisation criteria!');
end
        
disp(['APC ' propSpecification{1} ' propeller should be chosen for ' textOptimisation]);
disp([motorSpecification{1} ' (' num2str(round(motorSpecification{2}/10)*10) ' KV) motor should be selected with '...
    num2str(round(motorSpecification{4}*100)/100) ' Nm torque at maximum speed of ' num2str(round(motorSpecification{3}/100)*100) ' RPM.']);
disp(['The motor uses ' num2str(round(motorChosen{11})) ' W of electrical power at hover and ' num2str(round(motorChosen{8})) ' W of electrical power at WOT.']);
disp(['The drive should be controlled by a ' num2str(ceil(escSpecification)) ' A ESC per motor.']);
disp(['The whole system should be powered by a ' num2str(baterrySpecification(1)) 'S ' num2str(ceil(baterrySpecification(2))) 'C LiPo battery of '...
    num2str(baterrySpecification(3)) ' mAh.']);
disp('---------');
disp(['Hovering flight requires ' num2str(round(operatingPoints{temp_propChosen_pos,1}(4))) ' W of mechanical power (' num2str(round(operatingPoints{temp_propChosen_pos,1}(3)*100)/100)...
    ' Nm at ' num2str(round(operatingPoints{temp_propChosen_pos,1}(1)/100)*100) ' RPM) to achieve ' num2str(round(operatingPoints{temp_propChosen_pos,1}(2)*RotorNo)) ' gf of total thrust.']);
disp(['WOT flight requires ' num2str(round(operatingPoints{temp_propChosen_pos,2}(4))) ' W of mechanical power (' num2str(round(operatingPoints{temp_propChosen_pos,2}(3)*100)/100)...
    ' Nm at ' num2str(round(operatingPoints{temp_propChosen_pos,2}(1)/100)*100) ' RPM) to achieve ' num2str(round(operatingPoints{temp_propChosen_pos,2}(2)*RotorNo)) ' gf of total thrust.']);
disp(['This configuration should achieve around ' num2str(round(time_hover(end)*60)) ' min of hover and around ' num2str(round(time_max(end)*60)) ' min of flight at WOT.']);

plot_propPerf; % plot propeller performance & battery simulation results 
plot_motorPerf; % plot motor performance
