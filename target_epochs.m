function [epoch,epoch_r,rt,start_indices_nontarget] = target_epochs(eeg, marker_values ,marker_timestamps, fs, time)

% Epoch extraction (-200ms to 1000ms after stim onset)

format long g;

start_indices_target = [];
rt = [];start_indices_nontarget=[];
baseline_corrected = [];e1=[];
t = [];
index = [];

for i=1:1:size(marker_values,2)
    if marker_values(1,i) == 11 % target trial

        if marker_values(1,i+2) == 100 || marker_values(1,i+2) == 200
            format long g;
            start_indices_target = [start_indices_target marker_timestamps(1,i+1)];
        end
    elseif marker_values(1,i) == 12
        if marker_values(1,i+2) == 100 || marker_values(1,i+2) == 200
            format long g;
            start_indices_nontarget = [start_indices_nontarget marker_timestamps(1,i+1)];
            rt = [rt marker_timestamps(1,i+2)];

        end
    end

end

   epoch_r = []; 
epoch = []; 
b = round(0.2 * fs);           % 200 ms baseline
post_stim = round(2.0 * fs);   % 2000 ms after stimulus

for i = 1:length(start_indices_target)
    stim_time = start_indices_target(i);
    
    [~, idx] = min(abs(time - stim_time));
    
    if idx - b >= 1 && idx + post_stim - 1 <= size(eeg, 1)
        baseline = 0; full_segment=[];corrected=[];
        % Get baseline + post-stimulus segment
        full_segment = eeg(idx - b : idx + post_stim - 1, :);  % size: (b + post_stim) x channels
        
        % Baseline: first b samples
        baseline = mean(full_segment(1:b, :), 1);  % Check if the first b samples are clean
        
        % Correct only post-stimulus portion
        corrected = full_segment(b + 1:end, :) - baseline;  % size: post_stim x channels
        
        % Save only the post-stimulus corrected data
        epoch = cat(3, epoch, corrected);  % size: channels x time x trials
    end
end


b=0.2*fs;
for i = 1:size(start_indices_nontarget, 2)
    stim_timenon = start_indices_nontarget(i);
    
    % Find closest sample index corresponding to the stimulus time
    idxnon = find(time <= stim_timenon, 1, 'last');  % time in seconds
    
    % Ensure you don’t go out of bounds (start and end time checks)
    if idxnon - b >= 1 && idxnon + fs <= size(eeg, 1)
        % Extract the segment (epoch) around the stimulus time
        segmentnon = eeg(idxnon - b : idxnon + fs, :);  % shape: [epoch_len x channels]
        
        % Baseline correction: subtract mean of pre-stimulus period (-200 to 0 ms)
        baselinenon = mean(segmentnon(1:b, :), 1);  % mean across baseline time window
        correctednon  = segmentnon - baselinenon;      % subtract baseline from the entire segment
   
        % Save the corrected epoch to the 3D matrix (trials x channels x timepoints)
        epoch_r = cat(3, epoch_r, correctednon);
    end
end
end