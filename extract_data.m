function [eeg,time,ts,fs,m_j] = extract_data(streams)

fs = cell(1,size(streams,2)); 
ts = cell(1,size(streams,2)); 
eeg1 = cell(1,size(streams,2)); 
eeg = cell(1,size(streams,2)); 
time1 = cell(1,size(streams,2)); 
time = cell(1,size(streams,2)); 
m_j = cell(1,size(streams,2)); 


for j = 1:1:size(streams,2)

for i = 1:1:size(streams{1,j},2)
    data_type = string(streams{1,j}{1,i}.info.type);
    format long g

        if strcmp(data_type, 'Markers')  % Markers data stream
            if size(streams{1,j}{1,1}.time_series,2) > 10
            
             % all the markers are in streams{1,a}{1,1}
                m_j{j} = streams{1,j}{1,1}.time_series;
            
             end
        end

        if strcmp(data_type, 'EEG')  % Markers data stream
        ts{j} = str2double(streams{1,j}{1,i}.info.first_timestamp); %first time stamp
        fs{j} = streams{1,j}{1,i}.info.effective_srate; %sampling rate
        
        time1{j} = (streams{1,j}{1,i}.time_stamps); % time vector for EEG 
        eeg1{j} = (streams{1,j}{1,i}.time_series); % EEG data

       end
end
% we have the data as channels x time and we want them as time x channels
eeg{j} = eeg1{j}';
time{j} = time1{j}';
end



end
