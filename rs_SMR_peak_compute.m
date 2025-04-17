%To compute the SMR predictor, measure the difference between mu band power
%and the 1/f noise floor at c3 and C4 during rest.
%Compute the maximum difference between the PSD curve and a fit of the 1/f
%noise, and then average the values at C3 and C4
%close all
%clear
clc
% addpath(genpath('../PL_enhance_neuromod/code/matlab_code/functions/'))
% %PSD from 2-35Hz
% currentpath = pwd;
% sub = '033';
% session = 'post';
% path = ([currentpath '/f1_data/' session '/Subject_' sub '_tACSMI/' 'Subject_' sub '_tACSMI_Session_001/']);
% dirInfo = dir(path);
% rest_file = contains({dirInfo.name}, 'rest');
% eeg = sload([path '/' dirInfo(find(rest_file==1)).name '/' 'Subject' '_' sub '*.gdf']);
% %eeg = sload([currentpath '/' dirInfo(find(rest_file==1)).name '/' 'Subject' '_101*.gdf']);
% %using only 22 channels
load('ch32Locations.mat');
channel_names = {'C3', 'C4', 'FC5', 'FC1', 'CP5', 'CP1', 'Cz', 'FC2', 'FC6', 'CP2', 'CP6'};
channel_names = {'C3', 'C4'};
neighbors = {{9,10,20,21}; {11,12,22,23}; {4,5,15}; {5,6,15,16}; {15,24,25}; {15,16,25,26}; {10,11,21,22}; {6,7,16,17}; {7,8,17}; {16,17,26,27}; {17, 27,28} };
%%
eeg = eeg_t_ecpre{1,1};
figure;
for i = 1:length(channel_names)
    L = 512;
    noverlap = 256;
    f = 1:30;
    fs = 512;
    channels_id = find(ismember({ch32Locations.labels}, channel_names{i}));
    %spatial filter (The use of small laplacian)
    eeg(:,channels_id) = eeg(:,channels_id) - mean(eeg(:,neighbors(i,:)), 2);
    signal = eeg(15*fs:end-15*fs, channels_id);
    [pxx,freq] = pwelch(signal,hamming(L),noverlap, f,fs);
    subplot(3,4,i)
    plot(freq,10*log10(pxx)); hold on;
    % f = 2:0.5:35;
    %higher frequency resolution and lower variance
%      [pxx,freq] = pmtm(signal,[], f,fs);
%      figure;
%      plot(freq,10*log10(pxx)); 
    
    %Optimisation for 9 parameters
    modelFun =  @(p,x) p(8)/(p(1)* sqrt(2*pi)) .* exp(-1/2.*((x-p(2))./p(1)).^2) + p(9)/(p(3)* sqrt(2*pi)) .* exp(-1/2.*((x-p(4))./p(3)).^2) + p(5) + p(6)./(x.^p(7));
    %all coefficients are positive
    startingVals = [1 1 1 1 1 1 1 1 1];
    for k = 1:10 %this actually helps
        %nlModel = fitnlm(freq',pxx',modelFun,startingVals);
        nlModel = lsqcurvefit(modelFun,startingVals, freq',10*log10(pxx)');
        %fit it with actual data to see the accuracy of fit.
        %line(freq,10*log10(predict(nlModel,freq')),'Color','r');
        line(freq,modelFun(nlModel,freq'),'Color','r');
        startingVals = nlModel;
    end
    
    %Fsumsquares = @(p)sum((modelFun(p,freq) - pxx).^2);
    %[xunc,ressquared,eflag,outputu] = fminunc(Fsumsquares,startingVals,'Algorithm','quasi-newton');
    
    psd_estimate = modelFun(nlModel,freq');
    noise_estimate_raw = nlModel(5) + nlModel(6)./(freq'.^nlModel(7));
    noise_estimate = 10*log10(noise_estimate_raw);
    noise_estimate(find(noise_estimate_raw < 0)) = 0;
    signal_psd_estimate = psd_estimate - noise_estimate; %decibel normalisation
    peak_freq = freq(find(signal_psd_estimate == max(signal_psd_estimate)));
    plot(freq,signal_psd_estimate); hold on;
    title (['Analysis at ' channel_names{i}]);
end

%% RUN THIS
eeg = eeg_t_ecpre{1,1};

%close all
figure; %with small laplacian 
eeg_Lap = eeg(:,1:32);
eeg_pre = eeg_Lap;
for i = 1:length(channel_names)
    fs = 512;
    channels_id = find(ismember({ch32Locations.labels}, channel_names{i}));
    %spatial filter (The use of small laplacian)
    eeg_Lap(:,channels_id) = eeg_pre(:,channels_id) - mean(eeg_pre(:,cell2mat(neighbors{i})),2);
    signal = eeg_Lap(10*fs:end-5*fs, channels_id);
    subplot(3,4,i)
    [pxx,freq] = pwelch(signal,512,256,[1:30],512);
    %[pxx,freq] = pmtm(signal, [], [1:30], 512);
    plot(freq,10*log10(pxx)); hold on;
    ylim([-30 20]);
    
    %manually select points to estimate noise floor
%     [x_noise,y_noise] = getpts;
%     modelFun = @(p,x) p(1) + p(2)./(x.^p(3));
%     startingVals = [1 1 1];
%     for k = 1:1 %this actually helps
%         nlModel = lsqcurvefit(modelFun,startingVals, x_noise',y_noise'); 
%         startingVals = nlModel;
%     end
%     noise_floor_estimate = modelFun(nlModel,x_noise');
%     line(x_noise,noise_floor_estimate,'Color','b', 'DisplayName', 'noise floor initial estimate');
%     
%     modelFun =  @(p,x) p(3)/(p(1)* sqrt(2*pi)) .* exp(-1/2.*((x-p(2))./p(1)).^2) + p(4) + p(5)./(x.^p(6));
    % + p(7)/(p(8)* sqrt(2*pi)) .* exp(-1/2.*((x-p(9))./p(8)).^2);
    %Use the estimated noise floor coefficients as starting values in the
    %estimation
%     [x_alpha,y_alpha] = getpts;
%     [x_beta,y_beta] = getpts;
%     startingVals = [1 x_alpha 1 nlModel(1) nlModel(2) nlModel(3)];
%     %startingVals = [1 x_alpha 1 nlModel(1) nlModel(2) nlModel(3) 1 x_beta 1];
%     for k = 1:1 %this actually helps
%         options = optimoptions('lsqcurvefit','OptimalityTolerance', 1e-16, 'FunctionTolerance', 1e-16, 'MaxFunctionEvaluations',9e+08, 'MaxIterations', 4e+05, 'StepTolerance', 1e-20);
%         nlModel = lsqcurvefit(modelFun,startingVals, freq',10*log10(pxx)',[],[],options);
%         
%         startingVals = nlModel;
%     end
%     line(freq,modelFun(nlModel,freq'),'Color','r');
    %Fsumsquares = @(p)sum((modelFun(p,freq) - pxx).^2);
    %[xunc,ressquared,eflag,outputu] = fminunc(Fsumsquares,startingVals,'Algorithm','quasi-newton');
    
%     psd_estimate = modelFun(nlModel,freq');
%     %noise_estimate_raw = nlModel(5) + nlModel(6)./(freq'.^nlModel(7));
%     noise_estimate_raw = nlModel(4) + nlModel(5)./(freq'.^nlModel(6));
%     signal_psd_estimate = psd_estimate - noise_estimate_raw; 
%     line(freq,noise_estimate_raw,'Color','b', 'DisplayName', 'noise floor');
%     peak_freq = freq(find(signal_psd_estimate == max(signal_psd_estimate)));
    %plot(freq,signal_psd_estimate); hold on;
    title (['Analysis at ' channel_names{i} ' ;Estimated peak freq: ']);
end


figure; % with CAR
eeg_CAR = eeg(:,1:32);
eeg_pre = eeg_CAR;
for i = 2:length(channel_names)
    fs = 512;
    channels_id = find(ismember({ch32Locations.labels}, channel_names{i}));
    %spatial filter (The use of small laplacian)
    %eeg(:,channels_id) = eeg_pre(:,channels_id) - mean(eeg_pre(:,cell2mat(neighbors{i})),2);
    %spatial filter (CAR)
    channels_all = 1:32;
    channels_all(channels_id) = []; 
    eeg_CAR(:,channels_id) = eeg_pre(:,channels_id) - mean(eeg_pre(:,channels_all), 2);
    signal = eeg_CAR(15*fs:end-15*fs, channels_id);
    subplot(3,4,i)

    [pxx,freq] = pwelch(signal,512,256,[1:20],512);
    plot(freq,10*log10(pxx)); hold on;
    %ylim([-30 20]);
   %manually select points to estimate noise floor
    [x_noise,y_noise] = getpts;
    modelFun = @(p,x) p(1) + p(2)./(x.^p(3));
    startingVals = [1 1 1];
    for k = 1:1 %this actually helps
        nlModel = lsqcurvefit(modelFun,startingVals, x_noise',y_noise'); 
        startingVals = nlModel;
    end
    noise_floor_estimate = modelFun(nlModel,x_noise');
    line(x_noise,noise_floor_estimate,'Color','b', 'DisplayName', 'noise floor initial estimate');
    
    modelFun =  @(p,x) p(3)/(p(1)* sqrt(2*pi)) .* exp(-1/2.*((x-p(2))./p(1)).^2) + p(4) + p(5)./(x.^p(6));
    % + p(7)/(p(8)* sqrt(2*pi)) .* exp(-1/2.*((x-p(9))./p(8)).^2);
    %Use the estimated noise floor coefficients as starting values in the
    %estimation
    [x_alpha,y_alpha] = getpts;
    [x_beta,y_beta] = getpts;
    startingVals = [1 x_alpha 1 nlModel(1) nlModel(2) nlModel(3)];
    %startingVals = [1 x_alpha 1 nlModel(1) nlModel(2) nlModel(3) 1 x_beta 1];
    for k = 1:1 %this actually helps
        options = optimoptions('lsqcurvefit','OptimalityTolerance', 1e-16, 'FunctionTolerance', 1e-16, 'MaxFunctionEvaluations',9e+08, 'MaxIterations', 4e+05, 'StepTolerance', 1e-25);
        nlModel = lsqcurvefit(modelFun,startingVals, freq',10*log10(pxx)',[],[],options);
        
        startingVals = nlModel;
    end
    line(freq,modelFun(nlModel,freq'),'Color','r', 'DisplayName', 'estimate'); hold on
    
    %Fsumsquares = @(p)sum((modelFun(p,freq) - pxx).^2);
    %[xunc,ressquared,eflag,outputu] = fminunc(Fsumsquares,startingVals,'Algorithm','quasi-newton');
    
    psd_estimate = modelFun(nlModel,freq');
    noise_estimate_raw = nlModel(4) + nlModel(5)./(freq'.^nlModel(6));
    %noise_estimate_raw = nlModel(5) + nlModel(6)./(freq'.^nlModel(7));
    line(freq,noise_estimate_raw,'Color','b', 'DisplayName', 'noise floor');
    signal_psd_estimate = psd_estimate - noise_estimate_raw; 
    peak_freq = freq(find(signal_psd_estimate == max(signal_psd_estimate)));
    %plot(freq,signal_psd_estimate); hold on;
    title (['Analysis at ' channel_names{i} ' ;Estimated peak freq: ' num2str(peak_freq)]);
end
