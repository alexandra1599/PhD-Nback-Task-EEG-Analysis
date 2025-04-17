function laplacian_data = apply_laplacian_filter(eeg_data, channel_labels)
    % Apply Laplacian filter to Cz, CPz, and Pz channels
    %
    % Parameters:
    % eeg_data: 2D matrix (time x channels)
    % channel_labels: Cell array of channel names (1 x channels)
    %
    % Returns:
    % laplacian_data: 2D matrix (time x 3) with Laplacian-filtered Cz, CPz, and Pz

    % Define the target channels and their surrounding neighbors
    target_channels = {'CZ', 'CPZ', 'PZ'};
    neighbors = {
        {'C1', 'C2', 'FCZ', 'PZ'},       % Neighbors of Cz
        {'CP1', 'CP2', 'CZ', 'PZ'},      % Neighbors of CPz
        {'P1', 'P2', 'CPZ', 'POZ'}       % Neighbors of Pz
    };

    % Initialize the Laplacian data array
    [num_samples, ~] = size(eeg_data);
    laplacian_data = zeros(num_samples, length(target_channels));

    % Loop through each target channel (Cz, CPz, Pz)
    for i = 1:length(target_channels)
        target_channel = target_channels{i};
        neighbor_channels = neighbors{i};

        % Find indices of the target and neighboring channels
        target_idx = find(strcmp(channel_labels, target_channel));
        neighbor_indices = find(ismember(channel_labels, neighbor_channels));

        % Compute the Laplacian as the difference between the target and the mean of its neighbors
        if ~isempty(target_idx) && ~isempty(neighbor_indices)
            laplacian_data(:, i) = eeg_data(:, target_idx) - mean(eeg_data(:, neighbor_indices), 2);
        else
            warning(['Channel ' target_channel ' or its neighbors not found in data.']);
        end
    end
end
