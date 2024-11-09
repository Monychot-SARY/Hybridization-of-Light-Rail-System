function obj = Objective_Function(x, Pload, time)
    freq_cutoff = x(1); % Cutoff Frequency
    DOD = x(2);          % Depth of Discharge 
    max_capacity = 100;  % Maximum capacity for ONESS (kWh) - define based on your system requirements
    dt = mean(diff(time)); % Time step based on input `time` array
    %% Energy Management System 
    [Grid, ONESS] = EMS(Pload, freq_cutoff); % Assuming Cutoff Frequency in Low-Pass Filter is freq_cutoff
    conversion_efficiency = 0.9; % Assume that the conversion efficiency is 90%
    ONESS_Actual = ONESS * conversion_efficiency; % Adjust power for efficiency
    Total_Energy_ONESS = cumtrapz(ONESS_Actual)/DOD;
    % Total_Capacity_EONESS = max(Total_Energy_ONESS)/DOD;
    total_energy_from_grid = trapz(time, Grid) / 3600;
    alpha = 1000;
    obj = alpha * Total_Energy_ONESS + total_energy_from_grid; % Overall cost metric
end
