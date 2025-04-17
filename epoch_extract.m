function [epoch,epoch_r,start_indices,r_t] = epoch_extract(eeg, marker_values ,marker_timestamps, fs, time)

% Epoch extraction (-200ms to 1000ms after stim onset)

format long g;



start_indices = [];
s=0;
rt = []; r_t = [];
e = []; e1=[];


epoch = []; epoch_r=[];
b= 0.2*512;
for i=11:1:size(marker_values,2)
    if(marker_values(1,i) == 0)
    format long g
    s = marker_timestamps(1,i);
    start_indices = [start_indices,marker_timestamps(1,i)];
    b = 0.2*512;
    for j = 1:1:size(time,1)
            if time(j,1) <= s && time(j+1,1) >= s
                e = eeg(j-b:j+2*fs,:);
             
                epoch = cat(3,epoch,e);
            end
    end

    elseif((marker_values(1,i) == 100 || marker_values(1,i) == 200))
        format long g

        rt =  [rt, marker_timestamps(1,i)];
    end
end
  

for i = 1:1:size(rt,2)
    r_t = [r_t, (rt(1,i)-start_indices(1,i))*100];
end
end