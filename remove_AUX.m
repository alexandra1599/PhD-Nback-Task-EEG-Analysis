function [eeg,EOG] = remove_AUX(eeg,ch)

% 32 channel cap : M1 : 13, M2 : 19 + last 7 AUX
% 64 channel cap : M1 : 13, M2 : 19, EOG : 32, AUX : 65 - 71 
r = 7; % Number of AUX columns to remove


for j =1:1:size(eeg,2)
    
    VEOG{1,j} = eeg{1,j}(:,36); %AUX 7
    HEOG{1,j} = eeg{1,j}(:,37); %AUX8
    EOG{1,j} = cat(2,VEOG{1,j} , HEOG{1,j});
    eeg{1,j}(:,end-r+1:end) = []; % remove AUX
    
end

end