function bp = nback_proc(eeg,fs)

%notch filter to remove line noise
Wo = 60/(fs/2);  
BW = Wo/35;
[bn,an] = iirnotch(Wo,BW); 

filtnotch = filter(bn,an,eeg);

car_channels = [];
avg_channel = mean(filtnotch,2);
car_channels = filtnotch-avg_channel;


%Step 1 : filter the signal, P300 falls in [0.1-30] Hz CPP falls in
%[0.1-45]Hz

low = 1;
high = 40;
nyquist = 0.5 * fs;

[b,a] = butter(4,[low/nyquist,high/nyquist],'bandpass');
bp = filtfilt(b,a,car_channels);

end