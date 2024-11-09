clear all; close all; clc;

load Tram.mat

Pload = [T.pelec]/1000; % Convert power to kW
time = t';

% Calculate the time interval, assuming uniform time steps
dt = mean(diff(time));

% Compute energy by cumulative integration
Energy_IGBT = cumtrapz(Pload) / 3600; % Converting from Ws to kWh;
for i = 1:length(Pload)
    if Pload(i) > 0
        EnergyD(i) = Pload(i);
    else
        EnergyD(i) = 0;
    end
end
Energy_Diode = cumtrapz(EnergyD)/3600;

%% PPH and EPH
[PPH,PESS,EPH,Es,Eu]= PPH_EPH(Pload,time);

%% Frequency Spectrum
[S,k,fssdim,freelle,S2,k2,fece]=spectrece(PESS,1);

%% Energy Management System 

fc = freelle(find(S2 == max(S2)));
[Grid, ONESS] = EMS(Pload,fc,EPH); % Assuming that Cutoff Frequency in Low-Pass Filter is 0.5Hz
% [S,k,fssdim,freelle,S2,k2,fece]=spectrece(ONESS,1);
max(ONESS)
%% Sizing of the ONESS
%% a)

conversion_efficiency = 0.9; % Assume that the convert efficiency is 90%
ONESS_Actual = ONESS*conversion_efficiency; % Calculate actual power provided by ONESS considering efficiency
% Define the maximum capacity of ONESS in kWh
max_capacity = 100; % Example value; set this to your actual maximum capacity

% Initialize cumulative energy stored in kWh
E_ONESS = zeros(size(ONESS_Actual)); 

E_ONESS(1) = ONESS_Actual(1)* dt/ 3600; % Set initial stored energy

for t = 2:length(ONESS_Actual)
    % Calculate energy change based on whether ONESS is charging or discharging
    E_ONESS(t) = E_ONESS(t-1) + ONESS_Actual(t) * dt/ 3600; % Charge the ONESS
    if E_ONESS(t) > max_capacity
        E_ONESS(t) = max_capacity; % Cap at max capacity
    elseif E_ONESS(t) <= 0
        E_ONESS(t) = 0; % Prevent negative values
    end
    if E_ONESS(t) > 0
        En_ONESS(t) = E_ONESS(t); % Cap the energy at maximum capacity
    elseif E_ONESS(t) <= 0
        En_ONESS(t) = 0; % Prevent negative energy
    end
end
Total_EONESS = cumtrapz(En_ONESS);
%% c
Total_Capcity = round(max(Total_EONESS));
disp(['Total Capacity of ONESS: ', num2str(Total_Capcity), ' kWh']);

% Li-ion Battery
DOD_Li = 0.8; % Assume that the depth of discharge in the system is 80% rate of Li-ion battery DOD
Useful_energy_Li = max(Total_EONESS)/DOD_Li;
Total_Capacity_ONESS_Li_Dis = round(max(Useful_energy_Li));
disp(['Total Capacity of ONESS Discharge of Li-ion Battery: ', num2str(Total_Capacity_ONESS_Li_Dis), ' kWh']);

% NiMH Battery
DOD_NiMH = 0.95;
Useful_energy_NiMH = max(Total_EONESS)./DOD_NiMH;
Total_Capacity_ONESS_NiMH_Dis = round(max(Useful_energy_NiMH));
disp(['Total Capacity of ONESS Discharge of NiMH Battery: ', num2str(Total_Capacity_ONESS_NiMH_Dis), ' kWh']);

% Lead Fluid Battery
DOD_Lead = 0.5;
Useful_energy_Lead = max(Total_EONESS)/DOD_Lead;
Total_Capacity_ONESS_Lead_Dis = round(max(Useful_energy_Lead));
disp(['Total Capacity of ONESS Discharge of Lead Fluid Battery: ', num2str(Total_Capacity_ONESS_Lead_Dis), ' kWh']);

%% Battery or Supercapacitor
% Battery 
discharge_capacities = [Useful_energy_Li, Useful_energy_NiMH, Useful_energy_Lead, Total_Capcity];
bat_capacities_price = [Useful_energy_Li*130, Useful_energy_NiMH*350, Useful_energy_Lead*150];

Battery_Energy_Density = 120;  % Battery Energy Density: Typically around 120 Wh/kg of Nickel-Metal Hydride (NiMH) 
[~, min_index] = min(bat_capacities_price); % Get index of minimum price
Battery_mass = (discharge_capacities(min_index) * 1000) / Battery_Energy_Density;
% Assume that voltage cell is 3.7voltages; current cell is 2500mA
volt_cell_batt = 3.7; current_cell_batt = 2.5;
Num_Cells_batt= ((Total_Capcity*1000/DOD_NiMH))/(volt_cell_batt*current_cell_batt);

% Supercapacitor0
[Grid1, ONESS1] = EMS(Pload, fc, EPH);
conversion_efficiency = 0.9; 
ONESS_Actual1 = ONESS1 * conversion_efficiency; 
max_capacity = 100;
% Initialize cumulative energy storage in kWh
E_ONESS1 = zeros(size(ONESS_Actual1)); 
En_ONESS1 = zeros(size(ONESS_Actual1)); 
E_ONESS1(1) = ONESS_Actual1(1) * dt / 3600; 

for t = 2:length(ONESS_Actual1)
    % Calculate energy change based on charging or discharging
    E_ONESS1(t) = E_ONESS1(t-1) + ONESS_Actual1(t) * dt / 3600;
    % Cap energy storage to max capacity and prevent negatives
    if E_ONESS1(t) > max_capacity
        E_ONESS1(t) = max_capacity;
    elseif E_ONESS1(t) <= 0
        E_ONESS1(t) = 0;
    end
    % Update En_ONESS1 with capped energy values
    En_ONESS1(t) = E_ONESS1(t);
end

% Calculate cumulative energy over time (total energy processed)
Total_EONESS1 = cumsum(En_ONESS1);
Total_Capacity1 = round(max(Total_EONESS1));

% Calculate supercapacitor requirements
Capacitor_Energy_density = 10; % Energy Density in Wh/kg
Supercapacitor_mass = (Total_Capacity1 * 1000) / Capacitor_Energy_density;

Super_cells = 5; % Individual cell energy rating in Wh
Num_Cells_Super = (Total_Capacity1 * 1000) / Super_cells;


%% Sizing DC/DC Converter 
%% Size the wheel and the associated motor

RPM = 10000; % Assume the RPM of the motor

omega = RPM*(2*pi)/60; % 60 come from 60s/mn , Angular Velocity of the Flywheel
Energy_in_Joule = abs(min(ONESS_Actual))*3.6*1e6;
Inertia = (2*Energy_in_Joule)/(omega^2); %alculate Required Moment of Inertia  I kg⋅m 2
r = 0.5;  % r based on feasible design parameters and solve for m
m = (2*Inertia)/(r^2); % Flywheel Mass and Radius: Assuming a solid disk shape (kg) 
% a flywheel with a radius of 0.5 meters would require a mass of approximately m kg to store max(Energy_IGBT) kWh.
%% Sizing Motor
alpha = omega / length(time);
Troque = Inertia * alpha;
Max_T = max(Troque);
P_motor = Max_T *omega/1000;


%% Plot Power Consumption
figure();
plot(time, Pload, 'r', 'LineWidth', 1);
hold on;
time_max = find(Pload == max(Pload));
time_min = find(Pload == min(Pload));
plot(time(time_max), max(Pload), 'go', 'LineWidth', 1);
plot(time(time_min), min(Pload), 'mo', 'LineWidth', 1);
text(time(time_min), min(Pload), num2str(min(Pload)), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
text(time(time_max), max(Pload), num2str(max(Pload)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
title('Power of Electricity Consumption from the Trams');
xlim([0, time(end)]);
xlabel('Time (s)');
ylabel('Power Consumption (kW)');
legend('Power Consumption for 1 Vehicle');
box on;
%% Plot Energy Consumption
figure();
hold on;
plot(time, Energy_IGBT, 'y-', 'LineWidth', 1);
plot(time, Energy_Diode, 'g-.', 'LineWidth',1);
title('Cumulative Energy Consumption');
xlabel('Time (s)');
ylabel('Energy (kWh)');
xlim([0, time(end)]);
text(length(time), max(Energy_Diode), num2str(max(Energy_Diode)), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
text(length(time), max(Energy_IGBT), num2str(max(Energy_IGBT)), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
legend('Energy with IGBT Rectifier', 'Energy with Diode Rectifier', 'Location', 'Best');
box on

%% Plot EMS power sharing
figure()
plot(time, Grid, 'b', time, ONESS, 'r',LineWidth=1);
xlabel('Time (s)');
ylabel('Power (kW)');
xlim([0, time(end)]);
title('EMS Power Sharing between Grid and ONESS');
legend('Grid Power', 'ONESS Power');
box on;

%% Plot Frequency Spectrum
figure();
plot(freelle, S2, 'm', 'LineWidth', 1);
xlabel('Frequency (Hz)'); ylabel('Magnitude PESS');
title('Frequency Spectrum of PESS Power');
box on;

%% Plot PPH and EPH as Bar Plots with Correct Legend and Text
figure()
h1 = bar(1, PPH, 'FaceColor', 'b', 'EdgeColor', 'k'); hold on;
h2 = bar(2, EPH, 'FaceColor', 'r', 'EdgeColor', 'k');
set(gca, 'XTick', [1 2], 'XTickLabel', {'PPH', 'EPH'});
title('PPH and EPH Values');
text(1, PPH(1), sprintf('%.4f', PPH(1)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(2, EPH(1), sprintf('%.4f', EPH(1)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
box on;

%% Plot Energy Stored in ONESS
figure()
plot(time, En_ONESS, 'c', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Energy Stored in ONESS (kWh)');
xlim([0, time(end)]);
title('Energy in ONESS Over Time');
box on

%% Plot Energy Stored in ONESS
figure();
plot(time, Total_EONESS, 'g-.', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Energy Stored in ONESS (kWh)');
xlim([0, time(end)]);
title('Cumulative Energy Stored in ONESS Over Time');
text(length(time), max(Total_EONESS), num2str(max(Total_EONESS)), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
legend('Cumulative Energy in ONESS (kWh)', 'Location', 'Best');
box on

%% Plot DOD of Each battery, Cell and Weight
% Plot DOD of Each Battery, Cell, Mass, and Price
figure();

% Plot 1: Total Capacity of ONESS Discharge for different battery types
subplot(4, 1, 1); % Changed to 4 rows
b1 = bar(discharge_capacities);
b1.FaceColor = 'flat';
b1.CData(1,:) = [0.2 0.6 1.0]; % Li-ion color
b1.CData(2,:) = [1.0 0.5 0.0]; % NiMH color
b1.CData(3,:) = [0.3 0.7 0.4]; % Lead Fluid color
b1.CData(4,:) = [0.5 0.5 0.5]; % Total Capacity color
ylim([0, max(discharge_capacities)*1.1]);
set(gca, 'XTickLabel', {'Li-ion', 'NiMH', 'Lead Fluid', 'Total Capacity'});
ylabel('Capacity (kWh)');
title('Capacity for Different Battery Types');
text(1:length(discharge_capacities), discharge_capacities, ...
    num2str(discharge_capacities', '%.1f'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% Plot 2: Mass of Battery and Supercapacitor
subplot(4, 1, 2);
mass_values = [Battery_mass, Supercapacitor_mass];
b2 = bar(mass_values);
b2.FaceColor = 'flat';
b2.CData(1,:) = [0.6 0.2 0.7]; % Battery mass color
b2.CData(2,:) = [0.0 0.7 0.9]; % Supercapacitor mass color
ylim([0, max(mass_values) * 1.2]);
set(gca, 'XTickLabel', {'Battery', 'Supercapacitor'});
ylabel('Mass (kg)');
title('Mass Required for Battery and Supercapacitor');
text(1:length(mass_values), mass_values, ...
    num2str(mass_values', '%.1f'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% Plot 3: Number of Cells Required for Battery and Supercapacitor
subplot(4, 1, 3);
cell_counts = [Num_Cells_batt, Num_Cells_Super];
b3 = bar(cell_counts);
b3.FaceColor = 'flat';
b3.CData(1,:) = [0.8 0.3 0.3]; % Battery cells color
b3.CData(2,:) = [1.0 0.5 0.3]; % Supercapacitor cells color
ylim([0, max(cell_counts) * 1.2]);
set(gca, 'XTickLabel', {'Battery Cells', 'Supercapacitor Cells'});
ylabel('Number of Cells');
title('Number of Cells Required for Battery and Supercapacitor');
text(1:length(cell_counts), cell_counts, ...
    num2str(cell_counts', '%.1f'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% Plot 4: Price of Battery and Supercapacitor
subplot(4, 1, 4);
b4 = bar(bat_capacities_price);
b4.FaceColor = 'flat';
b4.CData(1,:) = [0.5 0.7 0.2]; % Li-ion price color
b4.CData(2,:) = [0.7 0.3 0.1]; % NiMH price color
b4.CData(3,:) = [0.2 0.5 0.3]; % Lead Fluid price color
ylim([0, max(bat_capacities_price) * 1.2]);
set(gca, 'XTickLabel', {'Li-ion Price', 'NiMH Price', 'Lead Fluid Price'});
ylabel('Price ($)');
title('Price for Battery Types');
text(1:length(bat_capacities_price), bat_capacities_price, ...
    num2str(bat_capacities_price', '%.1f'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% Adjust layout
sgtitle('ONESS System Analysis');
