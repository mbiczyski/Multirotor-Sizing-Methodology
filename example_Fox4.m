%% ------------------------------------------------------------------------
% Multirotor Sizing Methodology
% with Flight Time Estimation
%
% M. Biczyski, R. Sehab, G.Krebs, J.F. Whidborne, P. Luk
%
% example_Fox4.m - modified main.m file with data from
% Heliceo Fox4 quadrotor; for documentation refer to main.m
%% ------------------------------------------------------------------------

close all; clear; clc;
format compact; format shortG;

RPM2RAD = 2*pi/60;
RAD2RPM = 60/2/pi;

%% User parameters
RotorNo = 4;
OptimisationGoal = 'hover'; % hover - best specific thrust (g/W) at hover
                            % max - best specific thrust (g/W) at 100% throttle
                            % utilisation - maximum usable power range of propeller
ThrustWeightRatio = 2.5; % 2 - minimum
                       % 3 - payload transport
                       % 4 - survaillence
                       % 5+ - aerobatics / hi-speed video
                       % 7+ - racing
PropDiameter_Min = 10; % inch
PropDiameter_Max = 15; % inch
SafetyFactor = 1.0; % [1-2]
AcceptedTypes = {'MR' 'E' 'E-3' 'E-4'};       
                              % E	Electric
                              % F	Folding Blade (electric only)
                              % MR	Multi-Rotor (electric only)
                              % SF	Slow Fly (electric only)
                              % E-3	3-Blade
                              % E-4	4-Blade

BattCellNo = 6; %S 1P
BattCellVoltage = 3.7; % V per cell
BattCapacity = 2*5000; % mAh
BattPeukertConstant = 1.3; % for LiPo
BattVoltageSagConstant = 0.5/0.8*BattCellNo; % 0.5V decrease per cell in resting volatage for 80% DoD
BattHourRating = 1; % h

%% Mass data [g]
mass_Motor_Est = 200; % assumption
mass_ESC_Est = 100; % assumption
mass_Propeller_Est = 150;
mass_Battery = 700;

mass_NoDrive_Est = 0;
mass_Total_Est = 4000;

%% Filter propeller set
% propList = name, file, diameter (in), pitch (in), mass (g), speedLimit (RPM)
propList = load_propList();

propList_considered_bySize = propList(cellfun(@(x) x>=PropDiameter_Min && x<=PropDiameter_Max, propList(:,3)), :);
propList_considered_byMass = propList(cellfun(@(x) x<=mass_Propeller_Est, propList(:,5)), :);
propList_considered_byType = propList(cellfun(@(x) endsWith(x, AcceptedTypes), propList(:,1)), :);

temp_propList_considered_names = intersect(intersect(propList_considered_bySize(:,1), propList_considered_byMass(:,1)), propList_considered_byType(:,1));

propList_considered = propList(ismember(propList(:,1), temp_propList_considered_names),:);
consideredNo = size(propList_considered,1);

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
    propPerf(ii) = {load_propPerf(propList_considered{ii,2}, false)};
end

%% Calculate operating points
thrustHover_Est = mass_Total_Est/RotorNo;
thrustMax_Est = thrustHover_Est*ThrustWeightRatio;
for ii = 1:consideredNo
    speedHover = interp1(propPerf{ii}(2:end,2), propPerf{ii}(2:end,1), thrustHover_Est);
    speedMax = interp1(propPerf{ii}(2:end,2), propPerf{ii}(2:end,1), thrustMax_Est);
    speedLimit = propList_considered{ii,6};
    operatingPoints(ii,1) = {[speedHover thrustHover_Est interp1(propPerf{ii}(:,1), propPerf{ii}(:,4), speedHover) interp1(propPerf{ii}(:,1), propPerf{ii}(:,3), speedHover)]};
    operatingPoints(ii,2) = {[speedMax thrustMax_Est interp1(propPerf{ii}(:,1), propPerf{ii}(:,4), speedMax) interp1(propPerf{ii}(:,1), propPerf{ii}(:,3), speedMax)]};
    operatingPoints(ii,3) = {[speedLimit interp1(propPerf{ii}(:,1), propPerf{ii}(:,2), speedLimit) interp1(propPerf{ii}(:,1), propPerf{ii}(:,4), speedLimit)...
        interp1(propPerf{ii}(:,1), propPerf{ii}(:,3), speedLimit)]};
end

%% Select propeller
for ii = 1:consideredNo
    switch OptimisationGoal
        case 'hover'
            selectionCriterion(ii,1) = (operatingPoints{ii,1}(1)*2*pi/60)*operatingPoints{ii,1}(3);
            selectionCriterion(ii,2) = operatingPoints{ii,1}(4);
        case 'max'
            selectionCriterion(ii,1) = (operatingPoints{ii,2}(1)*2*pi/60)*operatingPoints{ii,2}(3);
            selectionCriterion(ii,2) = operatingPoints{ii,2}(4);
        case 'utilisation'
            selectionCriterion(ii,1) = (operatingPoints{ii,3}(1)*2*pi/60)*operatingPoints{ii,3}(3) - (operatingPoints{ii,2}(1)*2*pi/60)*operatingPoints{ii,2}(3);
            selectionCriterion(ii,2) = operatingPoints{ii,3}(4) - operatingPoints{ii,2}(4);
        otherwise
            error('ERROR! Wrong optimisation criteria!');
    end
end

methodError(:,1) = abs(selectionCriterion(:,2) - selectionCriterion(:,1));
methodError(:,2) = abs(selectionCriterion(:,2) - selectionCriterion(:,1))./abs(selectionCriterion(:,2));
for ii = 1:consideredNo
    if operatingPoints{ii,3}(4) < operatingPoints{ii,2}(4) || isnan(operatingPoints{ii,2}(1))
        selectionCriterion(ii,:) = inf;
    end
end

[~, temp_propChosen_pos] = min(mean(selectionCriterion,2));

if selectionCriterion(temp_propChosen_pos,2) == inf
    error('ERROR! No matching propeller found!');
end

%% Load & filter motor data
% motorList = ID, name, ILimit (A), mass (g), kV, Rm (Ohm),
            % op_Imax (A), op_powerMaxEl(W), op_effMax (%), op_IHover (A), op_powerHoverEl (W), op_effHover (%)
            
motorList = load_motorList(BattCellNo*BattCellVoltage, operatingPoints{temp_propChosen_pos,2}(1), operatingPoints{temp_propChosen_pos,2}(3),...
    operatingPoints{temp_propChosen_pos,1}(1), operatingPoints{temp_propChosen_pos,1}(3)*SafetyFactor,...
    mass_Motor_Est);

if size(motorList,1) < 1
    error('ERROR! No matching motor found!');
end

%% Select motor
switch OptimisationGoal
    case 'hover'
        [~, temp_motorChosen_pos] = min([motorList{:,11}]);
    case 'max'
        [~, temp_motorChosen_pos] = min([motorList{:,8}]); 
    case 'utilisation'
        [~, temp_motorChosen_pos] = min(abs([motorList{:,3}]-[motorList{:,7}])); 
    otherwise
        error('ERROR! Wrong optimisation criteria!');
end

motorChosen = motorList(temp_motorChosen_pos,:);

%% Determine drive specification
% propSpecification = name, diameter (in), pitch (in)
% motorSpecification = name, kV, speedMax (RPM), torqueMax (Nm), powerMax (W), powerMaxEl (W), EfficiencyMax(%), voltageNominal (V)
% escSpecification = currentMax (A)
% baterrySpecification = NoCells, C-rating, minCapacity (mAh)

propSpecification = propList_considered(temp_propChosen_pos,[1 3:4]);
motorSpecification = {motorChosen{2}, motorChosen{5}, operatingPoints{temp_propChosen_pos,2}(1), operatingPoints{temp_propChosen_pos,2}(3)*SafetyFactor,...
    operatingPoints{temp_propChosen_pos,2}(4)*SafetyFactor, motorChosen{8}, motorChosen{9}, BattCellNo*BattCellVoltage};
escSpecification = motorChosen{7};

BattCapacity_Ah = BattCapacity/1000;
minBattRating = escSpecification*RotorNo*SafetyFactor/BattCapacity_Ah;
baterrySpecification = [BattCellNo, minBattRating, BattCapacity];

mass_Propeller = propList_considered{temp_propChosen_pos,5};
mass_Motor = motorChosen{4};
mass_ESC = mass_ESC_Est;
mass_Total = mass_NoDrive_Est + RotorNo*(mass_Motor + mass_ESC + mass_Propeller) + mass_Battery;

%% Calculate initioal battery state
voltage_hover(1) = (BattCellNo*(BattCellVoltage+0.5));
voltage_max(1) = (BattCellNo*(BattCellVoltage+0.5));
current_hover(1) = motorChosen{11}*RotorNo/voltage_hover(1);
current_max(1) = motorChosen{8}*RotorNo/voltage_max(1);
capacity_hover(1) = (current_hover(1)^(1-BattPeukertConstant))*(BattHourRating^(1-BattPeukertConstant))*(BattCapacity_Ah^BattPeukertConstant);
capacity_max(1) = (current_max(1)^(1-BattPeukertConstant))*(BattHourRating^(1-BattPeukertConstant))*(BattCapacity_Ah^BattPeukertConstant);

%% Calculate next flight iterations
timeStep = 1/60/60; % 1 s
ii = 1;
while voltage_hover(ii) > BattCellVoltage*BattCellNo && ii*timeStep < 2
    voltage_hover(ii+1) = voltage_hover(1) - (BattVoltageSagConstant/capacity_hover(1))*(capacity_hover(1) - capacity_hover(ii));
    current_hover(ii+1) = motorChosen{11}*RotorNo/voltage_hover(ii+1);
    capacity_hover(ii+1) = (current_hover(ii+1)^(1-BattPeukertConstant))*(BattHourRating^(1-BattPeukertConstant))*(BattCapacity_Ah^BattPeukertConstant) - sum(current_hover(2:end)*timeStep);
    ii = ii+1;
end
time_hover = (0:ii-1)*timeStep;

ii = 1;
while voltage_max(ii) > BattCellVoltage*BattCellNo && ii*timeStep < 2
    voltage_max(ii+1) = voltage_max(1) - (BattVoltageSagConstant/capacity_max(1))*(capacity_max(1) - capacity_max(ii));
    current_max(ii+1) = motorChosen{8}*RotorNo/voltage_max(ii+1);
    capacity_max(ii+1) = (current_max(ii+1)^(1-BattPeukertConstant))*(BattHourRating^(1-BattPeukertConstant))*(BattCapacity_Ah^BattPeukertConstant) - sum(current_max(2:end)*timeStep);
    ii = ii+1;
end
time_max = (0:ii-1)*timeStep;

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
disp(['One motor uses ' num2str(round(motorChosen{11})) ' W of electrical power at hover and ' num2str(round(motorChosen{8})) ' W of electrical power at WOT.']);
disp(['The drive should be controlled by a ' num2str(ceil(escSpecification)) ' A ESC per motor.']);
disp(['The whole system should be powered by a ' num2str(baterrySpecification(1)) 'S ' num2str(ceil(baterrySpecification(2))) 'C LiPo battery of '...
    num2str(baterrySpecification(3)) ' mAh.']);
disp('---------');
disp(['Hovering flight requires ' num2str(round(RotorNo*operatingPoints{temp_propChosen_pos,1}(4))) ' W of mechanical power (' num2str(round(operatingPoints{temp_propChosen_pos,1}(3)*100)/100)...
    ' Nm at ' num2str(round(operatingPoints{temp_propChosen_pos,1}(1)/100)*100) ' RPM) to achieve ' num2str(round(operatingPoints{temp_propChosen_pos,1}(2)*RotorNo)) ' gf of total thrust.']);
disp(['WOT flight requires ' num2str(round(RotorNo*operatingPoints{temp_propChosen_pos,2}(4))) ' W of mechanical power (' num2str(round(operatingPoints{temp_propChosen_pos,2}(3)*100)/100)...
    ' Nm at ' num2str(round(operatingPoints{temp_propChosen_pos,2}(1)/100)*100) ' RPM) to achieve ' num2str(round(operatingPoints{temp_propChosen_pos,2}(2)*RotorNo)) ' gf of total thrust.']);
disp(['This configuration should achieve around ' num2str(round(time_hover(end)*600)/10) ' min of hover and around ' num2str(round(time_max(end)*600)/10) ' min of flight at WOT.']);

plot_propPerf;
plot_motorPerf;
