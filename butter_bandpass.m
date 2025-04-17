function eeg = butter_bandpass(data,lowcut, highcut, fs)
    % Design a Butterworth bandpass filter
    %
    % Parameters:
    % lowcut: Lower cutoff frequency (Hz)
    % highcut: Higher cutoff frequency (Hz)
    % fs: Sampling frequency (Hz)
    % order: The order of the filter (default: 4)
    %
    % Returns:
    % b, a: Filter coefficients for the bandpass filter
    
    order = 4; % Default filter order
    
    nyq = 0.5 * fs;  % Nyquist frequency
    low = lowcut / nyq;
    high = highcut / nyq;
    
    % Design the Butterworth bandpass filter using MATLAB's butter function
    [b, a] = butter(order, [low, high], 'bandpass');
    eeg = filtfilt(b,a,data);
end
