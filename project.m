% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         ||       -----      -----       ----- \
%         ||     ||     ||  ||     ||   ||      ||
%         ||     ||     ||  ||_____||   ||      ||
%         ||     ||     ||  ||     ||   ||      ||
%          ----    -----    ||     ||     ----- /
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prompt = "Subject ID : ";
ID =  input(prompt);
prompt1 = "Session (Relaxation / tACS) : ";
session = input(prompt1);
prompt2 = "Type of the session : ";
type = input(prompt2);
prompt3 = "Number of runs in the session : ";
runs = input(prompt3);
prompt4 = "Nback N (input 9 if it is not Nback data) : ";
N = input(prompt4);

if strcmp(session,'Relaxation')
    if strcmp(type, 'Eyes Closed post')

        streams_r_ecpost = load_project_data(ID,session,runs,type,N);

        [eeg,time_r_ecpost,ts_r_ecpost,fs_r_ecpost,m_r_ecpost] = extract_data(streams_r_ecpost);

        % Remove EOG and AUX
        [eeg_r_ecpost, eog_r_ecpost] = remove_AUX(eeg,32);

    elseif strcmp(type, 'Eyes Closed pre')

        streams_r_ecpre = load_project_data(ID,session,runs,type,N);

        [eeg,time_r_ecpre,ts_r_ecpre,fs_r_ecpre,m_r_ecpre] = extract_data(streams_r_ecpre);

        % Remove EOG and AUX
        [eeg_r_ecpre, eog_r_ecpre] = remove_AUX(eeg,32);

    elseif strcmp(type,'Nback')

        streams_r_nbpre = load_project_data(ID,session,runs,type,N);

        [eeg,time_r_nbpre,ts_r_nbpre,fs_r_nbpre,m_r_nbpre] = extract_data(streams_r_nbpre);

        % Remove EOG and AUX
        [eeg_r_nbpre, eog_r_nbpre] = remove_AUX(eeg,32);

    elseif strcmp(type,'Nback + relax')

        streams_r_nbpost = load_project_data(ID,session,runs,type,N);

        [eeg,time_r_nbpost,ts_r_nbpost,fs_r_nbpost,m_r_nbpost] = extract_data(streams_r_nbpost);

        % Remove EOG and AUX
        [eeg_r_nbpost,eog_r_nbpost] = remove_AUX(eeg,32);

    elseif strcmp(type,'Relax')

        streams_relax = load_project_data(ID,session,runs,type,N);

        [eeg,time_relax,ts_relax,fs_relax,m_relax] = extract_data(streams_relax);

        % Remove EOG and AUX
        [eeg_relax, eog_relax] = remove_AUX(eeg,32);


    end


elseif strcmp(session, 'tACS')
    if strcmp(type, 'Eyes Closed post tACS')

        streams_t_ecpost = load_project_data(ID,session,runs,type,N);

        [eeg,time_t_ecpost,ts_t_ecpost,fs_t_ecpost,m_t_ecpost] = extract_data(streams_t_ecpost);

        % Remove EOG and AUX
        [eeg_t_ecpost, eog_t_ecpost] = remove_AUX(eeg,32);

    elseif strcmp(type, 'Eyes Closed pre tACS')

        streams_t_ecpre = load_project_data(ID,session,runs,type,N);

        [eeg,time_t_ecpre,ts_t_ecpre,fs_t_ecpre,m_t_ecpre] = extract_data(streams_t_ecpre);

        % Remove EOG and AUX
        eeg_t_ecpre = remove_AUX(eeg,32);


    elseif strcmp(type, 'Eyes Closed find tACS')

        streams_find = load_project_data(ID,session,runs,type,N);

        [eeg,time_find,ts_find,fs_find,m_find] = extract_data(streams_find);

        % Remove EOG and AUX
        eeg_t_find = remove_AUX(eeg,32);

    elseif strcmp(type, 'Eyes Closed pre nothing')

        streams_t_ecpren = load_project_data(ID,session,runs,type,N);

        [eeg,time_t_ecpren,ts_t_ecpren,fs_t_ecpren,m_t_ecpren] = extract_data(streams_t_ecpren);

        % Remove EOG and AUX
        eeg_t_ecpren = remove_AUX(eeg,32);

    elseif strcmp(type, 'Eyes Closed post nothing')

        streams_t_ecpostn = load_project_data(ID,session,runs,type,N);

        [eeg,time_t_ecpostn,ts_t_ecpostn,fs_t_ecpostn,m_t_ecpostn] = extract_data(streams_t_ecpostn);

        % Remove EOG and AUX
        eeg_t_ecpostn = remove_AUX(eeg,32);

    elseif strcmp(type,'Nback')

        streams_t_nbpre = load_project_data(ID,session,runs,type,N);

        [eeg,time_t_nbpre,ts_t_nbpre,fs_t_nbpre,m_t_nbpre] = extract_data(streams_t_nbpre);

        % Remove EOG and AUX
        [eeg_t_nbpre, eog_t_nbpre] = remove_AUX(eeg,32);

    elseif strcmp(type,'Nback + tACS')

        streams_t_nbpost = load_project_data(ID,session,runs,type,N);

        [eeg,time_t_nbpost,ts_t_nbpost,fs_t_nbpost,m_t_nbpost] = extract_data(streams_t_nbpost);

        % Remove EOG and AUX
        [eeg_t_nbpost, eog_t_nbpost] = remove_AUX(eeg,32);

    elseif strcmp(type,'Nothing')

        streams_n = load_project_data(ID,session,runs,type,N);

        [eeg,time_n,ts_n,fs_n,m_n] = extract_data(streams_n);

        % Remove EOG and AUX
        eeg_n = remove_AUX(eeg,32);


    end

end

load("ErrP_cap_chan_file.mat");

%%
n=n+1;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          -----      -----    \\                  //   -------    -----
%         ||    ||  ||     ||   \\      //\\      //   ||        ||     ||
%         ||----    ||     ||    \\    //  \\    //    ||-----   || ----
%         ||        ||     ||     \\  //    \\  //     ||        ||    \\
%         ||          -----        \\//      \\//       -------  ||     \\
%          _____  ____   ____   ____   ____        ____  _____
%            |   |    | |____| |    | |____| |    |    |   |
%            |   |____| |      |____| |      |___ |____|   |
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('ch32Locations.mat');
fprintf("Frequency band specification for topoplots\n ");
prompt1 = "Please input the Low freq cutoff : ";
fl = input(prompt1);

prompt2 = "Please input the High freq cutoff : ";
fh = input(prompt2);

% ALPHA 7.5-12.5

[BLP_power, PSD_norm ] = processing(eeg_r_ecpre{1,1},ch32Locations,fs_r_ecpre{1,1},fl,fh,type);


%%

f_fc_all_mixed_pre = cell(1,15);f_fc_all_mixed_post = cell(1,15);
c_cp_all_mixed_pre = cell(1,15);c_cp_all_mixed_post = cell(1,15);
f_fc_all_fract_pre = cell(1,15);f_fc_all_fract_post = cell(1,15);
c_cp_all_fract_pre =cell(1,15);c_cp_all_fract_post =cell(1,15);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          -----      -----    \\                  //   -------    -----
%         ||    ||  ||     ||   \\      //\\      //   ||        ||     ||
%         ||----    ||     ||    \\    //  \\    //    ||-----   || ----
%         ||        ||     ||     \\  //    \\  //     ||        ||    \\
%         ||          -----        \\//      \\//       -------  ||     \\
%                        ____   ____   ____
%                       |____| |____  |     \
%                       |       ____| |_____/
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FRACTAL
srate = fs_t_ecpre{1,1}; % sampling frequency
movingwin = [3 1]; % [window size, sliding step] %sliding window: common
frange = [1 50];
win = int32(movingwin(1)*srate);
step = movingwin(2)*srate;

% load ECoG data from one sensor recorded in the left occipital of one
% macaque in eyes-closed awake state, totally 5 mins
% load('ECoG_data.mat');
[b,a] = butter(2, [1 100]./(srate/2), 'bandpass');

prompt = "What channels do you want ? (C/CP [1] or F/FC [2]) ?";
c = input(prompt);

if(strcmp(session,'tACS'))
    if c==1
        channel_eegpre =  cat(2,eeg_t_ecpre{1,1}(:,15:17),eeg_t_ecpre{1,1}(:,20:23));
        channel_eegpost =  cat(2,eeg_t_ecpost{1,1}(:,15:17),eeg_t_ecpost{1,1}(:,20:23));
    elseif c==2
        channel_eegpre =  eeg_t_ecpre{1,1}(:,4:12);
        channel_eegpost =  eeg_t_ecpost{1,1}(:,4:12);

    end
elseif (strcmp(session,'Relaxation'))

    if c==1
        channel_eegpre =  eeg_r_ecpre{1,1}(:,13:22);
        channel_eegpost =  eeg_r_ecpost{1,1}(:,13:22);
        %channel_eegrelax = eeg_relax{1,1}(:,13:22);
    elseif c==2
        channel_eegpre =  eeg_r_ecpre{1,1}(:,3:12);
        channel_eegpost =  eeg_r_ecpost{1,1}(:,3:12);
        %channel_eegrelax = eeg_relax{1,1}(:,3:12);

    end
end
data_pre = filtfilt(b, a, channel_eegpre);
data_post = filtfilt(b, a, channel_eegpost);
%data_relax = filtfilt(b, a, channel_eegrelax);

% separate fractal and oscillatory components using sliding window
nwinpre = floor((length(data_pre) - win)/step);
nwinpost = floor((length(data_post) - win)/step);
%nwinrelax = floor((length(data_relax) - win)/step);

sigpre = zeros(win,nwinpre);
sigpost = zeros(win,nwinpost);
%sigrelax = zeros(win,nwinrelax);

for i = 1 : nwinpre
    sigpre(:,i) = data_pre(ceil((i-1)*step)+1 : ceil((i-1)*step)+win);
end
for i = 1 : nwinpost
    sigpost(:,i) = data_post(ceil((i-1)*step)+1 : ceil((i-1)*step)+win);
end
%  for i = 1 : nwinrelax
%     sigrelax(:,i) = data_relax(ceil((i-1)*step)+1 : ceil((i-1)*step)+win);
% end
tic

Frac_pre = amri_sig_fractal(sigpre,srate,'detrend',1,'frange',frange);
Frac_pre.time = (0:step/srate:step*(nwinpre-1)/srate)';
Frac_post = amri_sig_fractal(sigpost,srate,'detrend',1,'frange',frange);
Frac_post.time = (0:step/srate:step*(nwinpost-1)/srate)';
% Frac_relax = amri_sig_fractal(sigrelax,srate,'detrend',1,'frange',frange);
%  Frac_relax.time = (0:step/srate:step*(nwinrelax-1)/srate)';

toc

% fitting power-law function to the fractal power spectra
Frange = [1,30]; % define frequency range for power-law fitting
Frac_pre = amri_sig_plawfit(Frac_pre,Frange);
Frac_post = amri_sig_plawfit(Frac_post,Frange);
%Frac_relax = amri_sig_plawfit(Frac_relax,Frange);


% show averaged fractal and oscillatory power spectrum
figure;

%subplot(2,1,1);
loglog(Frac_pre.freq,mean(Frac_pre.mixd,2),'b','LineWidth',2);
hold on
loglog(Frac_pre.freq,mean(Frac_pre.frac,2),'--b','LineWidth',1);
hold on
loglog(Frac_post.freq,mean(Frac_post.mixd,2),'r','LineWidth',2);
hold on
loglog(Frac_post.freq,mean(Frac_post.frac,2),'--r','LineWidth',1);
% hold on
% loglog(Frac_relax.freq,mean(Frac_relax.mixd,2),'g','LineWidth',2);
%  hold on
%  loglog(Frac_relax.freq,mean(Frac_relax.frac,2),'--g','LineWidth',1);

if(strcmp(session,'tACS') && ( strcmp(type,'Eyes Closed post nothing')))
    legend ('pre control mixed signal', 'pre control fractal signal', 'post control mixed signal' ...
        , ' post control fractal signal');

elseif(strcmp(session,'Relaxation'))
    legend ('pre relaxation mixed signal', 'pre relaxation fractal signal', 'post relaxation mixed signal' ...
        , ' post relaxation fractal signal'); %,'relaxation mixed signal','relaxation fractal signal');
elseif(strcmp(session,'tACS'))

    legend ('pre tACS mixed signal', 'pre tACS fractal signal', 'post tACS mixed signal' ...
        , ' post tACS fractal signal');

end

% % % subplot(2,1,2);
% plot(Frac_pre.freq, mean(Frac_pre.osci,2));
% hold on
% plot(Frac_post.freq, mean(Frac_post.osci,2));
% if(strcmp(session,'tACS'))
%     legend('pre tACS oscillations','post tACS oscillations')
% elseif(strcmp(session,'Relaxation'))
%     legend('pre relaxation oscillations','post relaxation oscillations')
% end
%legend('pre control oscillations','post control oscillations')
%legend('pre relaxation oscillations','post relaxation oscillations')
%title('averaged fractal and oscillatory power spectrum for F/FC channels ')


%% 

if c==1

    c_cp_all_mixed_pre{1,n} = mean(Frac_pre.mixd,2);
    c_cp_all_mixed_post{1,n} = mean(Frac_post.mixd,2);
    c_cp_all_fract_pre{1,n} = mean(Frac_pre.frac,2);
    c_cp_all_fract_post{1,n} = mean(Frac_post.frac,2);

elseif c==2

    f_fc_all_mixed_pre{1,n} = mean(Frac_pre.mixd,2);
    f_fc_all_mixed_post{1,n} = mean(Frac_post.mixd,2);
    f_fc_all_fract_pre{1,n} = mean(Frac_pre.frac,2);
    f_fc_all_fract_post{1,n} = mean(Frac_post.frac,2);
end


%%
channel_eegpre_c_cp_mixed=[];channel_eegpost_c_cp_mixed=[];
channel_eegpre_f_fc_mixed=[];channel_eegpost_f_fc_mixed=[];
channel_eegpre_c_cp_frac=[];channel_eegpost_c_cp_frac=[];
channel_eegpre_f_fc_frac=[];channel_eegpost_f_fc_frac=[];
srate = fs_r_ecpre{1,1}; % sampling frequency
movingwin = [3 1]; % [window size, sliding step] %sliding window: common
frange = [1 50];
win = int32(movingwin(1)*srate);
step = movingwin(2)*srate;
[b,a] = butter(2, [1 100]./(srate/2), 'bandpass');

for k=1:1:14
channel_eegpre_c_cp_mixed = cat(2,channel_eegpre_c_cp_mixed,c_cp_all_mixed_pre{1,k}(1:196,:));
channel_eegpost_c_cp_mixed = cat(2,channel_eegpost_c_cp_mixed,c_cp_all_mixed_post{1,k}(1:196,:));
channel_eegpre_f_fc_mixed = cat(2,channel_eegpre_f_fc_mixed,f_fc_all_mixed_pre{1,k}(1:196,:));
channel_eegpost_f_fc_mixed = cat(2,channel_eegpost_f_fc_mixed,f_fc_all_mixed_post{1,k}(1:196,:));

channel_eegpre_c_cp_frac = cat(2,channel_eegpre_c_cp_frac,c_cp_all_fract_pre{1,k}(1:196,:));
channel_eegpost_c_cp_frac = cat(2,channel_eegpost_c_cp_frac,c_cp_all_fract_post{1,k}(1:196,:));
channel_eegpre_f_fc_frac = cat(2,channel_eegpre_f_fc_frac,f_fc_all_fract_pre{1,k}(1:196,:));
channel_eegpost_f_fc_frac = cat(2,channel_eegpost_f_fc_frac,f_fc_all_fract_post{1,k}(1:196,:));

end


% show averaged fractal and oscillatory power spectrum
figure;

%subplot(2,1,1);
loglog(Frac_pre.freq,mean(channel_eegpre_c_cp_mixed,2),'b','LineWidth',2);
hold on
loglog(Frac_pre.freq,mean(channel_eegpre_c_cp_frac,2),'--b','LineWidth',1);
hold on
loglog(Frac_post.freq,mean(channel_eegpost_c_cp_mixed,2),'r','LineWidth',2);
hold on
loglog(Frac_post.freq,mean(channel_eegpost_c_cp_frac,2),'--r','LineWidth',1);

if(strcmp(session,'tACS') && ( strcmp(type,'Eyes Closed post nothing')))
    legend ('pre control mixed signal', 'pre control fractal signal', 'post control mixed signal' ...
        , ' post control fractal signal');

elseif(strcmp(session,'Relaxation'))
    legend ('pre relaxation mixed signal', 'pre relaxation fractal signal', 'post relaxation mixed signal' ...
        , ' post relaxation fractal signal'); %,'relaxation mixed signal','relaxation fractal signal');
elseif(strcmp(session,'tACS'))

    legend ('pre tACS mixed signal', 'pre tACS fractal signal', 'post tACS mixed signal' ...
        , ' post tACS fractal signal');

end
title("C/CP Electrodes")

figure;

%subplot(2,1,1);
loglog(Frac_pre.freq,mean(channel_eegpre_f_fc_mixed,2),'b','LineWidth',2);
hold on
loglog(Frac_pre.freq,mean(channel_eegpre_f_fc_frac,2),'--b','LineWidth',1);
hold on
loglog(Frac_post.freq,mean(channel_eegpost_f_fc_mixed,2),'r','LineWidth',2);
hold on
loglog(Frac_post.freq,mean(channel_eegpost_f_fc_frac,2),'--r','LineWidth',1);

if(strcmp(session,'tACS') && ( strcmp(type,'Eyes Closed post nothing')))
    legend ('pre control mixed signal', 'pre control fractal signal', 'post control mixed signal' ...
        , ' post control fractal signal');

elseif(strcmp(session,'Relaxation'))
    legend ('pre relaxation mixed signal', 'pre relaxation fractal signal', 'post relaxation mixed signal' ...
        , ' post relaxation fractal signal'); %,'relaxation mixed signal','relaxation fractal signal');
elseif(strcmp(session,'tACS'))

    legend ('pre tACS mixed signal', 'pre tACS fractal signal', 'post tACS mixed signal' ...
        , ' post tACS fractal signal');

end

title("F/FC Electrodes")




%% NBACK TASK ANALYSIS
if (strcmp(session, 'tACS'))
    epoch_p300_pre_tACStarget = cell(1,size(streams_t_nbpre,2));
    epoch_p300_post_tACStarget = cell(1,size(streams_t_nbpost,2));
    epoch_p300_pre_tACSnontarget = cell(1,size(streams_t_nbpre,2));
    epoch_p300_post_tACSnontarget = cell(1,size(streams_t_nbpost,2));
    reaction_time_pre = cell(1,size(streams_t_nbpre,2));
    reaction_time_post = cell(1,size(streams_t_nbpost,2));
    start_time_pre = cell(1,size(streams_t_nbpre,2));
    start_time_post = cell(1,size(streams_t_nbpost,2));

elseif (strcmp(session, 'Relaxation'))
    epoch_p300_pre_rtarget = cell(1,size(streams_r_nbpre,2));
    epoch_p300_post_rtarget = cell(1,size(streams_r_nbpost,2));
    epoch_p300_pre_rnontarget = cell(1,size(streams_r_nbpre,2));
    epoch_p300_post_rnontarget = cell(1,size(streams_r_nbpost,2));
    reaction_time_pre = cell(1,size(streams_r_nbpre,2));
    reaction_time_post = cell(1,size(streams_r_nbpost,2));
    start_time_pre = cell(1,size(streams_r_nbpre,2));
    start_time_post = cell(1,size(streams_r_nbpost,2));
end

k=1;
for r =1:1:5
    % STEP 1 : Extract epochs

    if (strcmp(session, 'tACS'))
        labels = get_label(streams_t_nbpost);


        % eeg_pre = detrend(eeg_t_nbpre{1,r}, 'constant');
        % eeg_post = detrend(eeg_t_nbpost{1,r}, 'constant');
        % eog_pre = detrend(eog_t_nbpre{1,r}(:,1:2), 'constant');  % only use VEOG & HEOG
        % eog_post = detrend(eog_t_nbpost{1,r}(:,1:2), 'constant');  % only use VEOG & HEOG
        % 
        % eeg_clean_pre = remove_eog_artifacts(eeg_pre, eog_pre);
        % eeg_clean_post = remove_eog_artifacts(eeg_post, eog_post);

        filtered_eeg_pre = filtering(eeg_t_nbpre{1,r},'tACS',labels);
        filtered_eeg_post = filtering(eeg_t_nbpost{1,r},'tACS',labels);

        [epochs_pre_target,epoch_pre_nontarget,t_pre,st_pre] = target_epochs(filtered_eeg_pre,m_t_nbpre{1,r}(1,:), m_t_nbpre{1,r}(2,:), 512, time_t_nbpre{1,r});
        [epochs_post_target,epoch_post_nontarget,t_post,st_post] = target_epochs(filtered_eeg_post,m_t_nbpost{1,r}(1,:), m_t_nbpost{1,r}(2,:), 512, time_t_nbpost{1,r});
        reaction_time_pre{1,k} = t_pre ;
        reaction_time_post{1,k} = t_post ;
        start_time_pre{1,k} = st_pre;
        start_time_post{1,k} = st_post;
        epoch_p300_pre_tACStarget{1,k} = epochs_pre_target;
        epoch_p300_post_tACStarget{1,k} = epochs_post_target;
        epoch_p300_pre_tACSnontarget{1,k} = epoch_pre_nontarget;
        epoch_p300_post_tACSnontarget{1,k} = epoch_post_nontarget;
        k = k+1;
    elseif (strcmp(session, 'Relaxation'))
        labels = get_label(streams_r_nbpost);

        % eeg_pre = detrend(eeg_r_nbpre{1,r}, 'constant');
        % eeg_post = detrend(eeg_r_nbpost{1,r}, 'constant');
        % eog_pre = detrend(eog_r_nbpre{1,r}(:,1:2), 'constant');  % only use VEOG & HEOG
        % eog_post = detrend(eog_r_nbpost{1,r}(:,1:2), 'constant');  % only use VEOG & HEOG
        %
        % eeg_clean_pre = remove_eog_artifacts(eeg_pre, eog_pre);
        % eeg_clean_post = remove_eog_artifacts(eeg_post, eog_post);

        filtered_eeg_pre = filtering(eeg_r_nbpre{1,r},'Relaxation',labels);
        filtered_eeg_post = filtering(eeg_r_nbpost{1,r},'Relaxation',labels);

        [epochs_pre_target,epoch_pre_nontarget,t_pre,st_pre] = target_epochs(filtered_eeg_pre,m_r_nbpre{1,r}(1,:), m_r_nbpre{1,r}(2,:), 512, time_r_nbpre{1,r});
        [epochs_post_target,epoch_post_nontarget,t_post,st_post] = target_epochs(filtered_eeg_post,m_r_nbpost{1,r}(1,:), m_r_nbpost{1,r}(2,:), 512, time_r_nbpost{1,r});
        reaction_time_pre{1,k} = t_pre;
        reaction_time_post{1,k} = t_post;
        start_time_pre{1,k} = st_pre;
        start_time_post{1,k} = st_post;

        epoch_p300_pre_rtarget{1,k} = epochs_pre_target;
        epoch_p300_post_rtarget{1,k} = epochs_post_target;
        epoch_p300_pre_rnontarget{1,k} = epoch_pre_nontarget;
        epoch_p300_post_rnontarget{1,k} = epoch_post_nontarget;
        k = k+1;
    end



end
%%
if (strcmp(session, 'tACS'))
    p300_pre = cat(3,epoch_p300_pre_tACStarget{1,1},epoch_p300_pre_tACStarget{1,2},epoch_p300_pre_tACStarget{1,3},epoch_p300_pre_tACStarget{1,4},epoch_p300_pre_tACStarget{1,5});
    p300_post = cat(3,epoch_p300_post_tACStarget{1,1},epoch_p300_post_tACStarget{1,2},epoch_p300_post_tACStarget{1,3},epoch_p300_post_tACStarget{1,4},epoch_p300_post_tACStarget{1,5});
    nop300_pre = cat(3,epoch_p300_pre_tACSnontarget{1,1},epoch_p300_pre_tACSnontarget{1,2},epoch_p300_pre_tACSnontarget{1,3},epoch_p300_pre_tACSnontarget{1,4},epoch_p300_pre_tACSnontarget{1,5});
    nop300_post = cat(3,epoch_p300_post_tACSnontarget{1,1},epoch_p300_post_tACSnontarget{1,2},epoch_p300_post_tACSnontarget{1,3},epoch_p300_post_tACSnontarget{1,4},epoch_p300_post_tACSnontarget{1,5});
elseif (strcmp(session, 'Relaxation'))
    p300_pre = cat(3,epoch_p300_pre_rtarget{1,1},epoch_p300_pre_rtarget{1,2},epoch_p300_pre_rtarget{1,3},epoch_p300_pre_rtarget{1,4});
    p300_post = cat(3,epoch_p300_post_rtarget{1,1},epoch_p300_post_rtarget{1,2},epoch_p300_post_rtarget{1,3},epoch_p300_post_rtarget{1,4});
    nop300_pre = cat(3,epoch_p300_pre_rnontarget{1,1},epoch_p300_pre_rnontarget{1,2},epoch_p300_pre_rnontarget{1,3},epoch_p300_pre_rnontarget{1,4},epoch_p300_pre_rnontarget{1,5});
    nop300_post = cat(3,epoch_p300_post_rnontarget{1,1},epoch_p300_post_rnontarget{1,2},epoch_p300_post_rnontarget{1,3},epoch_p300_post_rnontarget{1,4},epoch_p300_post_rnontarget{1,5});
end

%% INITIALIZATION

cpp_all_pre = cell(1,15);cpp_all_post = cell(1,15);
p300_all_pre = cell(1,15);p300_all_post = cell(1,15);
n=1;

%%
prompt = "CPP[1] or P300[2] analysis";
c = input(prompt);

if(strcmp(session,'tACS'))
    pz = 24;cz=15;c1=19; c2=21;
    labels = get_label(streams_t_nbpost);

elseif(strcmp(session,'Relaxation'))
    pz = 25;cz=15;cpz = 20; c1=19; c2=21;
    labels = get_label(streams_r_nbpost);

end

if (c==1)
    %rt_pre1 = []; rt_post1 = [];rt_pre2=[];rt_pre_all=[];rt_pre_avg=0;
    % nop300_pre=[];nop300_post=[];
    % nop300_pre = cat(3,epoch_p300_pre_tACSnontarget{1,3},epoch_p300_pre_tACSnontarget{1,4});
    % nop300_post = cat(3,epoch_p300_post_tACSnontarget{1,3},epoch_p300_post_tACSnontarget{1,4});

    cpp_pre = mean(mean(cat(2,nop300_pre(:,c1,:),nop300_pre(:,c2,:)),3),2);
    stdpre = std(cat(2,nop300_pre(:,c1,:),nop300_pre(:,c2,:)))/size(nop300_pre,3);
    stdpre = mean(stdpre,3);
    stdpre = mean(stdpre,2);

    cpp_post = mean(mean(cat(2,nop300_post(:,c1,:),nop300_post(:,c2,:)),3),2);
    stdpost = std(cat(2,nop300_post(:,c1,:),nop300_post(:,c2,:)))/size(nop300_post,3);
    stdpost = mean(stdpost,3);
    stdpost = mean(stdpost,2);

    rt_pre1 = (reaction_time_pre{1,3} - start_time_pre{1,3})*1000;
    rt_pre1 = rt_pre1(rt_pre1 <= 700);
    rt_pre2 = (reaction_time_pre{1,4} - start_time_pre{1,4})*1000;
    rt_pre2 = rt_pre2(rt_pre2 <= 700);
    rt_pre_all = cat(2,rt_pre1,rt_pre2);
    rt_pre_avg = mean(rt_pre_all,2);

    rt_post1 = (reaction_time_post{1,3} - start_time_post{1,3})*1000;
    rt_post1 = rt_post1(rt_post1 <= 800);
    rt_post2 = (reaction_time_post{1,4} - start_time_post{1,4})*1000;
    rt_post2 = rt_post2(rt_post2 <= 800);
    rt_post_all = cat(2,rt_post1,rt_post2);
    rt_post_avg = mean(rt_post_all,2);

    figure();
    plot(((smoothdata(cpp_pre))),'Color','b','LineWidth',2);
    hold on
    plot(((smoothdata(cpp_post))),'Color','r','LineWidth',2);
    hold on
    plot(((smoothdata(cpp_pre+mean(stdpre,3)))),'Color','b','LineWidth',0.5);
    hold on
    plot(((smoothdata(cpp_pre-mean(stdpre,3)))),'Color','b','LineWidth',0.5);
    hold on
    plot(((smoothdata(cpp_post+mean(stdpost,3)))),'Color','r','LineWidth',0.5);
    hold on
    plot(((smoothdata(cpp_post-mean(stdpost,3)))),'Color','r','LineWidth',0.5);
    hold on 
    % xline(rt_pre_avg,'Color','b','LineWidth',1.5);
    % hold on 
    % xline(397,'Color','r','LineWidth',1.5);
    grid on
    xlim([0 600])
    xlabel('Time (ms)');ylabel("Amplitude (uV)")
    legend('pre tACS target','post tACS target'); %,'pre relax non-target','post relax non-target');
    title('CPP Analysis with standard error - Avg CP1/CP2')


elseif(c==2)
    % % p300_pre=[];p300_post=[];pz_pre_target=[];stdpre=0;post_stim=[];
    % p300_pre = cat(3,epoch_p300_pre_tACStarget{1,4},epoch_p300_pre_tACStarget{1,5});
    % p300_post = cat(3,epoch_p300_post_tACStarget{1,4},epoch_p300_post_tACStarget{1,5});

    pz_pre_target = mean(p300_pre(:,pz,:),3);
    stdpre = std(p300_pre(:,pz,:))/size(p300_pre,3);

    pz_post_target = mean(p300_post(:,pz,:),3);
    stdpost = std(p300_post(:,pz,:))/size(p300_post,3);

    post_stim = round(2.0 * 512);   % 2000 ms after stimulus
   time_vector = (0:post_stim - 1) / 512 * 1000;  % in ms
    


    figure();
    plot(((smoothdata(pz_pre_target))),'Color','b','LineWidth',2);
    hold on
    plot(((smoothdata(pz_post_target))),'Color','r','LineWidth',2);
    hold on
    plot(((smoothdata(pz_pre_target+mean(stdpre,3)))),'Color','b','LineWidth',0.5);
    hold on
    plot(((smoothdata(pz_pre_target-mean(stdpre,3)))),'Color','b','LineWidth',0.5);
    hold on
    plot(((smoothdata(pz_post_target+mean(stdpost,3)))),'Color','r','LineWidth',0.5);
    hold on
    plot(((smoothdata(pz_post_target-mean(stdpost,3)))),'Color','r','LineWidth',0.5);
    grid on
    xlim([0 900])
    xlabel('Time (ms)');ylabel("Amplitude (uV)")
    legend('pre tACS target','post tACS target'); %,'pre relax non-target','post relax non-target');
    title('P300 Analysis with standard error')


end


%legend('pre tACS target','post tACS target'); %,'pre tACS non-target','post tACS non-target');


%% ALL SUBJECTS ANALYSIS 
    cpp_all_pre{1,n} = cpp_pre;
    cpp_all_post{1,n} = cpp_post;
    p300_all_pre{1,n} = pz_pre_target;
    p300_all_post{1,n} = pz_post_target;


    %%
cpp_pre=[];cpp_post=[];p300_pre=[];p300_post=[];
for i=1:1:n
    cpp_pre = cat(2,cpp_pre,cpp_all_pre{1,i});
    cpp_post = cat(2,cpp_post,cpp_all_post{1,i});
    % p300_pre = cat(2,p300_pre,p300_all_pre{1,i});
    % p300_post = cat(2,p300_post,p300_all_post{1,i});

end
p300_pre_avg = mean(cpp_pre,2);p300_post_avg = mean(cpp_post,2);

stdpre = std(cpp_pre)/(size(cpp_pre,2)*16);
stdpost = std(cpp_post)/(size(cpp_post,2)*16);
%time_vector = (0:post_stim - 1) / 512 * 1000;  % in ms
  figure();
    plot(((smoothdata(p300_pre_avg))),'Color','b','LineWidth',2);
    hold on
    plot(((smoothdata(p300_post_avg)))+0.2,'Color','r','LineWidth',2);
    hold on
    plot(((smoothdata(p300_pre_avg+mean(stdpre,2)))),'Color','b','LineWidth',0.5);
    hold on
    plot(((smoothdata(p300_pre_avg-mean(stdpre,2)))),'Color','b','LineWidth',0.5);
    hold on
    plot(((smoothdata(p300_post_avg+mean(stdpost,2))))+0.2,'Color','r','LineWidth',0.5);
    hold on
    plot(((smoothdata(p300_post_avg-mean(stdpost,2))))+0.2,'Color','r','LineWidth',0.5);
    grid on
    xlim([0 600])
    xlabel('Time (ms)');ylabel("Amplitude (uV)")
    legend('pre relax target','post relax target'); %,'pre relax non-target','post relax non-target');
    title('AVG P300 with standard error for all subject')



function eeg_clean = remove_eog_artifacts(eeg_data, eog_data)
% Detrend and normalize
eeg_data = detrend(eeg_data);
eog_data = detrend(eog_data(:,1:2));

% Normalize
eeg_z = zscore(eeg_data);
eog_z = zscore(eog_data);

% Regularization parameter (to avoid instability)
lambda = 1e-5;

% Ridge regression
C = eog_z' * eog_z + lambda * eye(size(eog_z,2));  % [2 x 2]
B = C \ (eog_z' * eeg_z);  % [2 x EEG_channels]

eog_proj = eog_z * B;  % [samples x EEG_channels]
eeg_clean_z = eeg_z - eog_proj;

% Restore original scaling
eeg_clean = eeg_clean_z .* std(eeg_data) + mean(eeg_data);
end
