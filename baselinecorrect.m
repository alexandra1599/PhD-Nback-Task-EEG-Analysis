function signal = baselinecorrect(eeg,fs)

   dbaseline = mean(eeg(1:0.1*fs,:,:),1);
   signal = eeg - dbaseline;

end