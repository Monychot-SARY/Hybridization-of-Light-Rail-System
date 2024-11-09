function [Grid, ONESS] = EMS(Pload, Freq_cutoff,PPH)
    Grid = lowpass(Pload, Freq_cutoff); % Assuming 1 Hz sampling rate
    % High-frequency component for ONESS
    ONESS = Pload - Grid;
    max_ONESS_power = max(PPH*max(Pload));
    ONESS(ONESS > max_ONESS_power) = max_ONESS_power; 
    ONESS(ONESS < -max_ONESS_power) = -max_ONESS_power; 
    Grid = Pload - ONESS;
end