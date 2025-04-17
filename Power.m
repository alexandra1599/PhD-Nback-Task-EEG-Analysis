prompt = "Subject ID : ";
ID =  input(prompt);


prompt1 = "Session ID : ";
session =  input(prompt1);
prompt2 = "Day # : "; 
day =  input(prompt2);
prompt3 = "Number of runs in the session : ";
runs = input(prompt3);
prompt4 = "Nback N (input 9 if it is not Nback data) : ";
N = input(prompt4);

load('ErrP_cap_chan_file.mat');
load('chanlocs64.mat');


streams = load_data(ID,session,runs,day,N);


[eeg,time,ts,fs,m] = extract_data(streams);

% Remove EOG and AUX 

eeg1 = remove_AUX(eeg,32);

%% Relaxation Analysis 

%Alpha = 8-12 so 7.5-12.5
%Beta = 13-30 so 12.5-30.5
%Theta = 4-8 so 3.5-8.5
fl = 7.5;
fh = 12.5;
[BLP_power, PSD_norm ] = processing(eeg1{1,1},chan,fs_r{1,1},fl,fh,session);


%[BLP_power_ec, PSD_norm_ec ] = processing(eeg1{1,1},chan,fs{1,1},fl,fh,session);

