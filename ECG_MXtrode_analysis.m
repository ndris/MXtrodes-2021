%% ECG data analysis recorded from MXtrode and Natus 
intan_read_RHS2000_DR % NOTE: use 2nd MXene data file, use Redo Natus data file (skin prep) 
%% rename relevant variables
amplifier_data_Natus=amplifier_data(1,:);
t_Natus=t;
fs = frequency_parameters.amplifier_sample_rate; % samples/s

%% keep only relevant variables
keep amplifier_data_MX t_MX fs amplifier_data_Natus t_Natus
%% plot raw data - MXene
figure() 
% plot(t_MX,amplifier_data_MX,'k')
plot(t_MX-t_MX(1)-0.6,filt_MX2(1,:)./1e3,'k') %filtered data 
title('MXene ECG')
xlabel('Time (s)') 
ylabel('Signal (mV)')
ylim([-0.6 2.7])
xlim([0 10.5])

%% plot raw data - Natus
figure() 
% plot(t_MX,amplifier_data_MX,'k')
plot(t_Natus-t_Natus(1),filt_Natus2(1,:)./1e3,'k') %filtered data 
title('Natus ECG')
xlabel('Time (s)') 
ylabel('Signal (mV)')
ylim([-0.6 2.7])
xlim([0 10.5])

%% filtering 
filt_MX=[];
filt_Natus=[];
d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',58,'HalfPowerFrequency2',62, ...
               'DesignMethod','butter','SampleRate',fs);
filt_MX= filtfilt(d,amplifier_data_MX);
filt_Natus= filtfilt(d,amplifier_data_Natus);

filt_MX2=[];
filt_Natus2=[];
[b,a]=butter(2,[0.1,150]/(fs/2),'bandpass');
filt_MX2=filtfilt(b,a,filt_MX);
filt_Natus2=filtfilt(b,a,filt_Natus);

%% detect peaks 
[MX_val,MX_ind]=findpeaks(filt_MX2,'minpeakheight',1500); 
[N_val,N_ind]=findpeaks(filt_Natus2,'minpeakheight',1500); 

%% align peaks, plot each waveform, and overlay average 
segs_MX=[];segs_N=[];
start=0.3*fs;
stop=0.5*fs;
for i=1:length(MX_ind)-1
    segs_MX(i,:)=filt_MX2(MX_ind(i)-start:MX_ind(i)+stop);
end
for i=2:length(N_ind)-1
    segs_N(i,:)=filt_Natus2(N_ind(i)-start:N_ind(i)+stop);
end
%%
figure() % overlaid waveforms
hold on 
for i=1:length(MX_ind)-1
plot((0:1/fs:(start+stop)/fs),segs_MX(i,:)./1000,'color',[0.75 0.75 0.75])
plot((0:1/fs:(start+stop)/fs),(mean(segs_MX,1))./1000,'k','linewidth',1.5)
end
ylim([-1 3])
title('ECG waveforms - MXene') 
xlabel('Time (s)')
ylabel('Signal (mV)')

figure() % overlaid waveforms
hold on 
for i=1:length(N_ind)-1
plot((0:1/fs:(start+stop)/fs),segs_N(i,:)./1000,'color',[0.75 0.75 0.75])
plot((0:1/fs:(start+stop)/fs),(mean(segs_N,1))./1000,'k','linewidth',1.5)
end 
ylim([-1 3])
title('ECG waveforms - Natus') 
xlabel('Time (s)')
ylabel('Signal (mV)')
