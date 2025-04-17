function filtered_data = apply_notch_filter(eeg_data, sampling_rate, line_freq, quality_factor, harmonics)
    % Apply a notch filter to the EEG dataset to remove line noise.
    % 
    % Parameters:
    % eeg_data: EEG data of shape (n_samples x n_channels)
    % sampling_rate: Sampling rate of the EEG data in Hz
    % line_freq: Line noise frequency to filter (default: 50 Hz)
    % quality_factor: Quality factor of the notch filter (default: 30)
    % harmonics: Number of harmonics to filter (default: 1, filters only the line frequency)
    %
    % Returns:
    % filtered_data: Filtered EEG data of the same shape as input
    
    if nargin < 3
        line_freq = 60; % Default line noise frequency (60 Hz)
    end
    if nargin < 4
        quality_factor = 30; % Default quality factor
    end
    if nargin < 5
        harmonics = 1; % Default number of harmonics to filter
    end
    
    % Initialize the filtered data
    filtered_data = eeg_data;
    
    % Apply the notch filter for each harmonic
    for harmonic = 1:harmonics
        target_freq = line_freq * harmonic;
        
        % Calculate the notch filter coefficients using the MATLAB function 'iirnotch'
        [b, a] = iirnotch(target_freq / (sampling_rate / 2), target_freq / (sampling_rate / 2) / quality_factor);
        
        % Apply the notch filter along the time dimension (rows)
        filtered_data = filtfilt(b, a, filtered_data);
    end
end
