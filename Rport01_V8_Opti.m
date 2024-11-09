clear all; close all; clc;

load Tram.mat
% Convert power to kW
Pload = [T.pelec] / 1000; 
time = t';
dt = mean(diff(time));

% Optimization options for PSO
options = optimoptions('particleswarm', ...
                       'SwarmSize', 50, ...
                       'MaxIterations', 100, ...
                       'Display', 'iter', ...
                       'FunctionTolerance', 1e-6);

% Battery types, parameters, and specific cutoff frequencies
battery_types = {'Li-ion', 'NiMH', 'Lead Fluid'};
cost_per_type = []; % To store costs for each battery type
optimal_params_per_type = {}; % To store optimal parameters for each type

% Define specific cutoff frequencies for each battery type
cutoff_frequencies = [0.01, 0.01, 0.01]; % Example cutoff frequencies for Li-ion, NiMH, and Lead Fluid
a = 0.001; b = 0.5;

for i = 1:length(battery_types)
    battery_type = battery_types{i};
    
    % Set bounds and initial parameters based on battery type
    switch battery_type
        case 'Li-ion'
            lb = [a, 0.2, 0.1]; % Lower bounds for Li-ion
            ub = [b, 0.8, 1]; % Upper bounds for Li-ion
            price_per_kWh = 130; % in USD
        case 'NiMH'
            lb = [a, 0.8, 0.1]; % Lower bounds for NiMH
            ub = [b, 1, 1]; % Upper bounds for NiMH
            price_per_kWh = 350; % average between $200 and $500
        case 'Lead Fluid'
            lb = [a, 0.2, 0.1]; % Lower bounds for Lead Fluid
            ub = [b, 0.5, 1]; % Upper bounds for Lead Fluid
            price_per_kWh = 150; % average between $100 and $200
    end
    
    % Run PSO optimization for the current battery type
    [optimal_params, optimal_cost] = particleswarm(...
        @(x) Objective_Function(x, Pload, time, price_per_kWh), ...
        3, lb, ub, options);
    
    % Store results
    result{i} = optimal_params;
    [PPH, PESS] = PPH1(Pload);
    [S, k, fssdim, freelle, S2, k2, fece] = spectrece(PESS, optimal_params(1));
    freq_cutoff = freelle(find(S2 == max(S2)));
    conversion_efficiency = 0.9;
    [grid_power, oness_power] = EMS(Pload, freq_cutoff, optimal_params(3)); % EMS function call
    oness_power_actual = oness_power * conversion_efficiency;
    max_capacity = 100; % Example value; set this to your actual maximum capacity

    % Initialize cumulative energy stored in kWh
    E_ONESS = zeros(size(oness_power_actual)); 
    
    E_ONESS(1) = oness_power_actual(1) * dt / 3600; % Set initial stored energy
    
    for t = 2:length(oness_power_actual)
        % Calculate energy change based on whether ONESS is charging or discharging
        E_ONESS(t) = E_ONESS(t-1) + oness_power_actual(t) * dt / 3600; % Charge the ONESS
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

    % Calculate cumulative energy in ONESS and maximum capacity required
    total_energy_oness = cumtrapz(En_ONESS); % kWh    
    total_capacity_oness = max(total_energy_oness) / optimal_params(2); % kWh
    
    % Calculate energy drawn from the grid
    total_energy_from_grid = cumtrapz(time, grid_power) / 3600; % kWh
    total_capacity{i} = [total_capacity_oness, max(total_energy_from_grid), max(total_energy_oness)];
    
    % Cost calculation for the ONESS capacity
    oness_cost = total_capacity_oness * price_per_kWh;
    cost_per_type(i) = oness_cost;
    EPH_in(i) = optimal_params(3);
    Battery_Capacity(i) = total_capacity_oness;
    Grid{i} = [grid_power];
    ONESS{i} = [En_ONESS];
    optimal_params_per_type{i} = [freq_cutoff, optimal_params(2), optimal_params(3)];
    fprintf('Battery Type: %s\n', battery_type);
    fprintf('Optimal Cutoff Frequency (Hz): %.8f\n', freq_cutoff);
    fprintf('Optimal Depth of Discharge (DoD): %.8f\n', optimal_params(2));
    fprintf('Optimal Energy Potential Hybridization (EPH): %.8f\n', optimal_params(3));
    fprintf('Cost for this configuration: %.4f\n\n', cost_per_type(i));
end

%% Select the Best Battery Type
[~, idx_min_cost] = min(cost_per_type);

best_battery_type = battery_types{idx_min_cost};
best_params = optimal_params_per_type{idx_min_cost};
best_grid = Grid{idx_min_cost};
best_ONESS = ONESS{idx_min_cost};

fprintf('Best Battery Type: %s\n', best_battery_type);
fprintf('Optimal Parameters:\n');
fprintf(' - Cutoff Frequency (Hz): %.4f\n', best_params(1));
fprintf(' - Depth of Discharge (DoD): %.4f\n', best_params(2));
fprintf(' - Energy Potential Hybridization (EPH): %.4f\n', best_params(3));
fprintf(' - Minimum Cost: %.2f\n', cost_per_type(idx_min_cost));

%%
EPH = best_params(3);
[Grid, ONESS] = EMS(Pload,best_params(1),EPH);
figure()
plot(time, Grid, 'r', time, ONESS, 'c',LineWidth=1);
xlabel('Time (s)');
ylabel('Power (kW)');
xlim([0, time(end)]);
title('EMS Power Sharing between Grid and ONESS');
legend('Grid Power', 'ONESS Power');
box on;
%%

figure();
plot(freelle, S2, 'b', 'LineWidth', 1);
xlabel('Frequency (Hz)'); ylabel('Magnitude PESS');
title('Frequency Spectrum of PESS Power');
box on;

%%
%%
[S,k,fssdim,freelle2,S22,k2,fece]=spectrece(Grid,optimal_params(1));
figure();
plot(freelle2, S22, 'r', 'LineWidth', 1);
xlabel('Frequency (Hz)'); ylabel('Magnitude PESS');
title('Frequency Spectrum of Grid Power');
box on;

%% Plot PPH and EPH as Bar Plots with Correct Legend and Text
% Plot EPH of each battery type
figure;
bar(EPH_in, 'FaceColor', 'flat'); % Use the stored EPH values for plotting
xlabel('Battery Type');
ylabel('Energy Potential Hybridization (EPH)');
title('EPH of Each Battery Type');
set(gca, 'XTickLabel', {'Li-ion', 'NiMH', 'Lead Fluid'});
ylim([0, max(EPH_in) * 1.2]); % Set y-limit for better visualization
for i = 1:length(EPH_in)
    text(i, EPH_in(i), sprintf('%.4f', EPH_in(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end
box on;

%% Plot Energy Stored in ONESS
figure()
plot(time, En_ONESS, 'c', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Energy Stored in ONESS (kWh)');
xlim([0, time(end)]);
title('Energy in ONESS Over Time');
box on

%%
DOD = [optimal_params_per_type{1}(2)*100, optimal_params_per_type{2}(2)*100, optimal_params_per_type{3}(2)*100];
Capacity_at_DOD = Battery_Capacity .* (DOD / 100);

figure;

subplot(1, 4, 1);
b1 = bar(cost_per_type);
b1.FaceColor = 'flat';
b1.CData(1, :) = [0.2 0.6 1.0];
b1.CData(2, :) = [1.0 0.5 0.0];
b1.CData(3, :) = [0.3 0.7 0.4];
title('Cost per Battery Type');
xlabel('Battery Type');
ylabel('Cost (USD)');
set(gca, 'XTickLabel', {'Li-ion', 'NiMH', 'Lead Fluid'});
ylim([0, max(cost_per_type) * 1.1]);
text(1:length(cost_per_type), cost_per_type, num2str(cost_per_type', '%.1f'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

subplot(1, 4, 2);
b3 = bar(DOD);
b3.FaceColor = 'flat';
b3.CData(1, :) = [0.2 0.6 1.0];
b3.CData(2, :) = [1.0 0.5 0.0];
b3.CData(3, :) = [0.3 0.7 0.4];
title('Depth of Discharge');
xlabel('Battery Type');
ylabel('DOD (%)');
set(gca, 'XTickLabel', {'Li-ion', 'NiMH', 'Lead Fluid'});
ylim([0, 110]);
text(1:length(DOD), DOD, num2str(DOD', '%.1f'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
[max_dod_value, idx_max_dod] = max(DOD);
text(idx_max_dod, max_dod_value + 2, sprintf('Max DOD: %.1f%%', max_dod_value), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold', 'Color', 'black');

subplot(1, 4, 3);
b2 = bar(Battery_Capacity);
b2.FaceColor = 'flat';
b2.CData(1, :) = [0.2 0.6 1.0];
b2.CData(2, :) = [1.0 0.5 0.0];
b2.CData(3, :) = [0.3 0.7 0.4];
hold on; % Hold on to overlay the Battery Capacity line
h = plot(1:3, Capacity_at_DOD, 'k-o', 'MarkerFaceColor', 'red'); % Battery Capacity line
hold off; % Release hold
title('Battery Capacity');
xlabel('Battery Type');
ylabel('Capacity (kWh)');
set(gca, 'XTickLabel', {'Li-ion', 'NiMH', 'Lead Fluid'});
ylim([0, max(Battery_Capacity) * 1.1]);
text(1:length(Battery_Capacity), Battery_Capacity, num2str(Battery_Capacity', '%.1f'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% Add label for the Battery Capacity line
legend(h, 'Capacity of DOD', 'Location', 'northwest');

subplot(1, 4, 4);
b4 = bar(Capacity_at_DOD);
b4.FaceColor = 'flat';
b4.CData(1, :) = [0.2 0.6 1.0];
b4.CData(2, :) = [1.0 0.5 0.0];
b4.CData(3, :) = [0.3 0.7 0.4];
title('DOD Capacity');
xlabel('Battery Type');
ylabel('Capacity (kWh)');
set(gca, 'XTickLabel', {'Li-ion', 'NiMH', 'Lead Fluid'});
ylim([0, max(Battery_Capacity) * 1.1]);
% Annotating DOD Capacity values
text(1:length(Capacity_at_DOD), Capacity_at_DOD, num2str(Capacity_at_DOD', '%.1f'), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
%%
function obj = Objective_Function(x, Pload, time, price_per_kWh)
    % Extract parameters
    cutoff_freq = x(1); % Cutoff Frequency
    DOD = x(2);         % Depth of Discharge (fraction)
    EPH = x(3);         % Energy Potential Hybridization
    
    % Constants
    conversion_efficiency = 0.9; % Assume a 90% efficiency
    dt = mean(diff(time));
    max_capacity = 100; % Example max capacity in kWh
    
    % Energy Management System
    [grid_power, oness_power] = EMS(Pload, cutoff_freq, EPH); % EMS function
    oness_power_actual = oness_power * conversion_efficiency;
    
    % Initialize cumulative energy stored in kWh
    E_ONESS = zeros(size(oness_power_actual)); 
    E_ONESS(1) = oness_power_actual(1) * dt / 3600;
    
    % Compute energy storage changes
    for t = 2:length(oness_power_actual)
        E_ONESS(t) = E_ONESS(t-1) + oness_power_actual(t) * dt / 3600;
        % Apply constraints on ONESS capacity
        if E_ONESS(t) > max_capacity
            E_ONESS(t) = max_capacity; % Cap at max capacity
        elseif E_ONESS(t) <= 0
            E_ONESS(t) = 0; % Prevent negative energy
        end
    end
    
    % Calculate metrics
    total_energy_oness = cumtrapz(E_ONESS); % Cumulative ONESS energy (kWh)
    total_capacity_oness = max(total_energy_oness) / DOD; % Capacity (kWh)
    total_energy_from_grid = trapz(time, grid_power) / 3600; % Grid energy (kWh)
    
    % ONESS cost
    oness_cost = total_capacity_oness * price_per_kWh;
    
    % Objective Function
    c1 = 1;    % Weight for ONESS cost
    c2 = 100;  % Weight for grid energy cost
    c3 = 1000;
    obj = c1 * (oness_cost ) + c2 * (total_energy_from_grid) + c3*total_capacity_oness ; 
    
    % Add penalty for exceeding capacity
    penalty = sum(E_ONESS > max_capacity) * 1000; % Penalize infeasibility
    obj = obj + penalty;
end

%% Energy Management System Function
function [grid_power, oness_power] = EMS(Pload, freq_cutoff, PPH)
    grid_power = lowpass(Pload, freq_cutoff, 1); % Assuming 1 H sampling rate
    % High-frequency component for ONESS
    oness_power = Pload - grid_power;
    max_ONESS_power = max(PPH * max(Pload));
    
    % Allow ONESS to capture more energy
    oness_power(oness_power > max_ONESS_power) = max_ONESS_power; 
    oness_power(oness_power < -max_ONESS_power) = -max_ONESS_power; 
    grid_power = Pload - oness_power;
    
    grid_power = double(grid_power); % Ensure output is of type double
    oness_power = double(oness_power); % Ensure output is of type double
end

%%
function [PPH,PESS] = PPH1(x)
%% PPH
    for i = 1:length(x)
        if max(x) > 0
           PPH(i) = 1-(mean(x)/max(x));
        else 
           PPH(i) = 1;
        end
    end

    for i = 1:length(x)
        PESS(i) = x(i) - mean(x);
    end
end

function [S,k,fssdim,freelle,S2,k2,fece]=spectrece(xx,fece)

    % [S,k]=spectre(ce,fe);
    %
    % Matlab function to estimate the spectrum of a function xx
    % 
    % Inputs :
    %  - xx  : sampled function [N x 1]
    %  - fece : sampling frequency (in Hertz)
    %
    % Outputs :
    %  - S : spectrum value    [Nfft x 1]
    %  - k : Sampling elements of the spectrum [Nfft x 1]
    
    
    xx=xx(:);
    N=length(xx);
    g=hamming(N);% Hamming window
    
    xx=xx.*g;
    Nfft=8*length(xx);  
    TFD=fft(xx,Nfft); 
    S=abs(TFD);      
    
    %% Calibration of the abscissa
    S2=S(1:Nfft/2+1);
    k=0:Nfft-1; 
    k2=k(1:Nfft/2+1);
    fssdim=k2/Nfft; 
    freelle=k2*fece/Nfft; 
    
    S2=2*S2/abs(sum(g));
end
