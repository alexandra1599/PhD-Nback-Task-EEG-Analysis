function [epoch_i,epoch_c,t,rt_correct,rt_incorrect] = correct_epochs(eeg, marker_values ,marker_timestamps, fs, time)

% Epoch extraction (-200ms to 1000ms after stim onset)

format long g;

start_correct = []; start_incorrect = [];
rt_correct = [];rt_incorrect = [];
e = [];e1=[];
t = [];
index = [];

    for i=1:1:size(marker_values,2)
        if marker_values(1,i) == 11 && marker_values(1,i+3) == 1 % target + correct trial
            
            if marker_values(1,i+2) == 100 || marker_values(1,i+2) == 200
                format long g;
                start_correct = [start_correct marker_timestamps(1,i+1)];
            
                rt_correct = [rt_correct marker_timestamps(1,i+2)];
            end

        elseif marker_values(1,i) == 11 && marker_values(1,i+3) == 2 % incorrect trial
            
            if marker_values(1,i+2) == 100 || marker_values(1,i+2) == 200
                format long g;
                start_incorrect = [start_incorrect marker_timestamps(1,i+1)];
            
                rt_incorrect = [rt_incorrect marker_timestamps(1,i+2)];
            end
        end

    end

epoch_i = []; epoch_c=[];
b= 0.2*512;
for i = 1:1:size(start_correct,2)
    for j = 1:1:size(time,1)
        if time(j,1) <= start_correct(1,i) && time(j+1,1) >= start_correct(1,i)
            t = [t time(j,1)];
            index= [index j];
            e = eeg(j-b:j+2.5*fs,:);
            epoch_c = cat(3,epoch_c,e(b:end,:));
             
        end
    end
end
for i = 1:1:size(start_incorrect,2)
    for j = 1:1:size(time,1)
        if time(j,1) <= start_incorrect(1,i) && time(j+1,1) >= start_incorrect(1,i)

            e1 = eeg(j-b:j+2.5*fs,:);
            epoch_i = cat(3,epoch_i,e1);

        
        end
    end
end


end