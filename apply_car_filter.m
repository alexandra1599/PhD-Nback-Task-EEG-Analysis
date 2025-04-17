function filtered_data = apply_car_filter(data)
    % Apply Common Average Reference (CAR) filter to EEG data
    %
    % Parameters:
    % data: EEG data matrix (time x channels)
    %
    % Returns:
    % filtered_data: EEG data with CAR applied (time x channels)
    
    % Calculate the average of the data across channels (axis 2)
    avg = mean(data, 2);
    
    % Subtract the average from the data (Common Average Reference)
    filtered_data = data - avg;
end
