function filtered_signal = filtering(eeg,session,labels)
ds = 6;     % downsample to 40 Hz
fs = 512;
%notch filter to remove line noise
Wo = 60/(fs/2);  
BW = Wo/35;
[bn,an] = iirnotch(Wo,BW); 
fprintf("Notch Applied\n");
filtnotch = filter(bn,an,eeg);

if (strcmp(session,'tACS'))
    % Remove M1,M2
    fprintf("Removing M1,M2\n");

    remove_M1_M2 = filtnotch;
    remove_M1_M2(:,19) = [];
    remove_M1_M2(:,13) = [];
    labels{1,19} = []; labels{1,13}=[];
else 
    remove_M1_M2 = filtnotch;
end

%car
avg_channel = mean(remove_M1_M2,2);
car_channels = [];

for i=1:1:size(filtnotch,2)
        car_channels(:,i) = filtnotch(:,i)-avg_channel;
end
fprintf("CAR applied\n");
low = 1;
high = 30;
nyquist = 0.5 * fs;


[b,a] = butter(4,[low high]/nyquist,'bandpass');
bp = filtfilt(b,a,car_channels);

filtered_signal = bp;
fprintf(sprintf("BP applied in %3d - %3d Hz : ",low,high));

end 