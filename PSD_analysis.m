%% PSD analysis EC Relax

prompt = "Subject ID : ";
ID =  input(prompt);
prompt1 = "Session ID : ";
session = input(prompt1);
prompt2 = "Day # : "; 
day =  input(prompt2);
prompt3 = "Number of runs in the session : ";
runs = input(prompt3);
prompt4 = "Nback N (input 9 if it is not Nback data) : ";
N = input(prompt4);

streams_ec = load_data(ID,session,runs,day,N);

[eeg,time_nr,ts_nr,fs_nr,m_nr] = extract_data(streams_ec);

% Remove EOG and AUX 

eeg1 = remove_AUX(eeg,32);
eeg_ec = eeg1{1,1};
% eeg_ecr = eeg1{1,2};
% eeg_find = eeg1{1,3};

%%
prompt1 = "Session ID : ";
session = input(prompt1);
prompt2 = "Day # : "; 
day =  input(prompt2);
prompt3 = "Number of runs in the session : ";
runs = input(prompt3);
prompt4 = "Nback N (input 9 if it is not Nback data) : ";
N = input(prompt4);

streams_r1 = load_data(ID,session,runs,day,N);

[eeg,time_r,ts_r,fs_r,m_r] = extract_data(streams_r1);

% Remove EOG and AUX 

eeg2 = remove_AUX(eeg,32);
eeg_r1 = eeg2{1,1};
% eeg_r2 = eeg2{1,2};


%%

[A_norm,A,f] = psd_plot(eeg_t_ecpostn{1,1},512);
%[A1_norm,A1,f1] = psd_plot(eeg_r1,512);
% [A2_norm,A2,f2] = psd_plot(eeg_r2,512);
%[A3_norm,A3,f3] = psd_plot(eeg_ecr,512);
% [A1_norm,A1,f1] = psd_plot(eeg_find,512);


%labels = get_label(streams_ec);

 Aec_f_fc = mean(A(:,1:12),2);
 %Ar1_f_fc = mean(A1(:,1:12),2);
% Ar2_f_fc = mean(A2(:,1:12),2);
% Aecr_f_fc = mean(A3(:,1:12),2);
 % Aectacs_f_fc = mean(A1(:,1:12),2);

 Aec_c_cp = mean(A(:,13:22),2);
% Ar1_c_cp = mean(A1(:,13:22),2);
% Ar2_c_cp = mean(A2(:,13:22),2);
% Aecr_c_cp = mean(A3(:,13:22),2);
 % Aectacs_c_cp = mean(A1(:,13:22),2);

 Aec_p = mean(A(:,23:32),2);
 %Ar1_p = mean(A1(:,23:32),2);
% Ar2_p = mean(A2(:,23:32),2);
% Aecr_p = mean(A3(:,23:32),2);
 % Aectacs_p = mean(A1(:,23:32),2);


% Aec_f_fc = mean(A_norm(:,1:12),2);
% Ar1_f_fc = mean(A1_norm(:,1:12),2);
% Ar2_f_fc = mean(A2_norm(:,1:12),2);
% Aecr_f_fc = mean(A3_norm(:,1:12),2);
% 
% Aec_c_cp = mean(A_norm(:,13:22),2);
% Ar1_c_cp = mean(A1_norm(:,13:22),2);
% Ar2_c_cp = mean(A2_norm(:,13:22),2);
% Aecr_c_cp = mean(A3_norm(:,13:22),2);
% 
%Aec_p = mean(A_norm(:,23:32),2);
% Ar1_p = mean(A1_norm(:,23:32),2);
% Ar2_p = mean(A2_norm(:,23:32),2);
%Aecr_p = mean(A3_norm(:,23:32),2);

%%

figure();
loglog(f,smoothdata(Aec_f_fc),'Color','b','LineWidth',2)
hold on 
% loglog(f1,smoothdata(Ar1_f_fc),'Color','r','LineWidth',2)
% hold on
% loglog(f2,smoothdata(Ar2_f_fc),'Color','g','LineWidth',2)
% hold on
%loglog(f1,smoothdata(Aectacs_f_fc),'Color','r','LineWidth',2)
%hold on
%loglog(f3,smoothdata(Aecr_f_fc),'Color','black','LineWidth',2)
legend("Eyes Closed pre-tACS",'find tACS','Eyes Closed post-tACS')

% legend("Eyes Closed pre-relax",'Relax 1','Relax 2','Eyes Closed post-relax')
xlabel('Frequency (Hz)')
ylabel('dB/Hz')
grid on
%title("Normalized PSD for F/FC electrodes")

title("PSD for F/FC electrodes")
xlim([8 30])


figure();
loglog(f,smoothdata(Aec_c_cp),'Color','b','LineWidth',2)
hold on 
% loglog(f1,smoothdata(Ar1_c_cp),'Color','r','LineWidth',2)
% hold on
% loglog(f2,smoothdata(Ar2_c_cp),'Color','g','LineWidth',2)
% hold on
%loglog(f1,smoothdata(Aectacs_c_cp),'Color','r','LineWidth',2)
%hold on
%loglog(f3,smoothdata(Aecr_c_cp),'Color','black','LineWidth',2)
legend("Eyes Closed pre-tACS",'find tACS','Eyes Closed post-tACS')

% legend("Eyes Closed pre-relax",'Relax 1','Relax 2','Eyes Closed post-relax')
xlabel('Frequency (Hz)')
ylabel('dB/Hz')
grid on

%title("Normalized PSD for C/CP electrodes")

 title("PSD for C/CP electrodes")

xlim([8 30])


figure();
loglog(f,smoothdata(Aec_p),'Color','b','LineWidth',2)
hold on 
% loglog(f1,smoothdata(Ar1_p),'Color','r','LineWidth',2)
% hold on
% loglog(f2,smoothdata(Ar2_p),'Color','g','LineWidth',2)
% hold on
% loglog(f1,smoothdata(Aectacs_p),'Color','r','LineWidth',2)
% hold on
% loglog(f3,smoothdata(Aecr_p),'Color','black','LineWidth',2)
% % legend("Eyes Closed pre-relax",'Relax 1','Relax 2','Eyes Closed post-relax')
legend("Eyes Closed pre-tACS",'find tACS','Eyes Closed post-tACS')

xlabel('Frequency (Hz)')
ylabel('dB/Hz')
grid on
%title("Normalized PSD for P/PO electrodes")

title("PSD for P/PO electrodes")

xlim([8 30])
