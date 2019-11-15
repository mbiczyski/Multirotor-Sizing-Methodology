%% ------------------------------------------------------------------------
% Multirotor Sizing Methodology
% with Flight Time Estimation
%
% M. Biczyski, R. Sehab, G.Krebs, J.F. Whidborne, P. Luk
%
% load_motorList.m - implements function load_motorList() that reads basic
% motor information from the database and calculates its operating points; 
% requires file ''DCbase.dcd'' from Drive Calculator software
%% ------------------------------------------------------------------------

function motorList = load_motorList(voltage, prop_speedMax, prop_torqueMax, prop_speedHover, prop_torqueHover, spec_mass)
    motorList = {};
    conn = sqlite('DCbase.dcd'); % connect to the database
    motors = fetch(conn, 'SELECT * FROM Motors'); % run SQL query to obtain full 'Motors' table
    mdata = fetch(conn, 'SELECT * FROM MData'); % run SQL query to obtain full 'MData' table
    esc = fetch(conn, 'SELECT * FROM ESC'); % run SQL query to obtain full 'ESC' table
    close(conn); % close connection to database

    %% ---------------------
    for ii = 1:size(motors,1) % filter motors in the database based on maximum current, mass and maximum speed
        motor_id = motors{ii,1}; % obtain motor ID
        esc_id = motors{ii,16}; % obtain ID of ESC recommended for the motor
        performance = mdata([mdata{:,3}]==motor_id,:); % look for test data for the motor with given ID
        if isempty(performance) || motors{ii,5} > 0 || motors{ii,15} > 1 % if no test data is found, skip motor
            continue;
        end
        
        current_max = motors{ii,13}; % obtain maximum motor current
        if current_max < 1 % if no maximum current data is available, obtain it from ESC's maximum current or guess based on voltage
            if esc_id > 0 && esc{esc_id,5} > 0
                current_max = esc{esc_id,5};
            else
                current_max = 2*voltage;
            end
        end
        current = 1:current_max; % prepare current domain
        [~, noLoad_pos] = min(abs([performance{[performance{:,7}]==1,4}]-voltage)); % find no-load test data close to target voltage
        Rm = 2*motors{ii,18}; % read motor windings resitance
        kV = motors{ii,17}; % read motor KV rating
        mass = motors{ii,14}; % read motor mass
        
        ironLoss = voltage*performance{noLoad_pos,6}; % calculate iron losses based on no-load current
        voltageMax = max([performance{:,4}]); % obtain max. voltage used for tests

        %% ---------------------
        prop_powerMax = prop_torqueMax*prop_speedMax/60*2*pi; % calculate required propeller power at WOT
        prop_powerHover = prop_torqueHover*prop_speedHover/60*2*pi; % calculate required propeller power at hover
        motor_currentMax = (voltage - sqrt(voltage^2-4*Rm*(ironLoss+prop_powerMax)))/(2*Rm); % calculate motor current at WOT
        motor_currentHover = (voltage - sqrt(voltage^2-4*Rm*(ironLoss+prop_powerHover)))/(2*Rm); % calculate motor current at hover
        
%         disp([motor_id motor_currentMax current_max mass voltage voltageMax*1.15 0.8*voltage*kV prop_speedMax]);
        if isreal(motor_currentMax) && isreal(motor_currentHover) && motor_currentMax > 0 && motor_currentHover > 0 && motor_currentMax <= current_max &&...
                mass <= spec_mass && mass > 0 && 0.8*voltage*kV > prop_speedMax % filter motor
            motor_powerMaxEl = voltage*motor_currentMax; % calculate electrical power at WOT
            motor_effMax = prop_powerMax/(voltage*motor_currentMax)*100; % calculate efficiency at WOT
            motor_powerHoverEl = voltage*motor_currentHover; % calculate electrical power at WOT
            motor_effHover = prop_powerHover/(voltage*motor_currentHover)*100; % calculate efficiency at WOT
           
            % motorList = ID, name, ILimit (A), mass (g), kV, Rm (Ohm),
            % op_Imax (A), op_powerMaxEl (W), op_effMax (%), op_IHover (A), op_powerHoverEl (W), op_effHover (%)
            motorList(end+1,:) = {motor_id, motors{ii,3}, motors{ii,13}, mass, kV, Rm,...
                                  motor_currentMax, motor_powerMaxEl, motor_effMax, motor_currentHover, motor_powerHoverEl, motor_effHover};
                              
        end
    end
end