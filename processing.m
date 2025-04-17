function [BLP_power, PSD_norm,A,f ] = processing(eeg,ch,fs,fl,fh,rec)

fprintf("PROCESSING THE DATA\n");
% get the broadband signal

if strcmp(rec,'Eyes Closed pre')
    t = 'Eyes Closed pre-Relaxation'; % Eyes closed data
elseif strcmp(rec, 'Relax') 
    t = 'Relaxation'; % Relaxation data
elseif strcmp(rec, 'Eyes Closed post') 
    t = 'Eyes Closed post-Relaxation'; % Relaxation data
elseif strcmp(rec, 'Eyes Closed post tACS') 
    t = 'Eyes Closed post-tACS'; % tACS data
elseif strcmp(rec, 'Eyes Closed pre tACS') 
    t = 'Eyes Closed pre-tACS'; % tACS data
end

ds = 6;     % downsample to 40 Hz

%notch filter to remove line noise
Wo = 60/(fs/2);  
BW = Wo/35;
[bn,an] = iirnotch(Wo,BW); 
fprintf("Notch Applied\n");
filtnotch = filter(bn,an,eeg);

% Remove M1,M2
fprintf("Removing M1,M2\n");

remove_M1_M2 = filtnotch;
remove_M1_M2(:,19) = [];
remove_M1_M2(:,13) = [];

%car
avg_channel = mean(remove_M1_M2,2);
car_channels = [];

for i=1:1:size(filtnotch,2)
        car_channels(:,i) = filtnotch(512:end,i)-avg_channel(512:end,:);
end
fprintf("CAR applied\n");
low = 0.1;
high = 30;
nyquist = 0.5 * fs;


[b,a] = butter(4,[low high]/nyquist,'bandpass');
bp = filtfilt(b,a,car_channels);
bp=bp(1:ds:end,:);      % downsample, size(X)=[28,16]

bp(:,19)=zeros(size(bp,1),1);bp(:,13)=zeros(size(bp,1),1);


fprintf(sprintf("BP applied in %s - %s Hz : ",fl,fh));

fprintf("Exrtacting the power and plotting . . .\n ");

[A,f] = pwelch(bp,8*fs,4*fs,0.1:0.1:128,fs);
%figure();
%loglog(f,A)
BLP_power = sum(A(f>=fl & f<= fh,:),1);
figure();
topoplot(BLP_power,ch,'maplimits','maxmin','whitebk','on');
colorbar
if fh == 12.5
title(sprintf("Band limited power in the alpha band during '%s'",t));
elseif fh == 30.5
title(sprintf("Band limited power in the beta band during '%s'",t));
elseif fh == 8.5
title(sprintf("Band limited power in the theta band during '%s'",t));
end

% get the normalized power to see the contribution of every channel to the
% % overall power

AA = A./sum(A(:,:),1);
PSD_norm = sum(AA(f>=fl & f<=fh,:),1);

figure();
topoplot(PSD_norm,ch,'maplimits','maxmin','whitebk','on');
colorbar
if fh == 12.5
title(sprintf("Normalized power in the alpha band during '%s'",t));
elseif fh == 30.5
title(sprintf("Normalized power in the beta band during '%s'",t));
elseif fh == 8.5
title(sprintf("Normalized power in the theta band during '%s'",t));
end


end
