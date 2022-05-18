signal = 'wnoise';
fs = 9600;
freqLimits = [1 4800];
filtMode = 'filtfilt';
filtOrder = 10;
bpo = 3;
dispProgress = true;

% Load Signal
s = load(signal);
x = s.(signal);
t = (0:length(x)-1)*1/fs;

% Bandpass Filtering
AnsiFilterBank = ansiBankFilterDesign(fs,bpo,'FrequencyLimits',freqLimits,...
    'FilterOrder',filtOrder);
[xb0,fc] = ansiBankFilter(AnsiFilterBank,x,'MetricsOutput',0,'ZeroPadding',false,...
    'FilterMode',filtMode,'DisplayProgress',dispProgress); 
[xb1,~] = ansiBankFilter(AnsiFilterBank,x,'MetricsOutput',1,'ZeroPadding',false,...
    'FilterMode',filtMode,'DisplayProgress',dispProgress); 
[xb0_zp,~] = ansiBankFilter(AnsiFilterBank,x,'MetricsOutput',0,'ZeroPadding',true,...
    'FilterMode',filtMode,'DisplayProgress',dispProgress); 
[xb1_zp,~] = ansiBankFilter(AnsiFilterBank,x,'MetricsOutput',1,'ZeroPadding',true,...
    'FilterMode',filtMode,'DisplayProgress',dispProgress); 
Xb1 = 20*log10(xb1(1,:));
Xb1_zp = 20*log10(xb1_zp(1,:));

% Plotting
figure
hold on
plot(fc,Xb1,'b')
plot(fc,Xb1_zp,'m')
xlabel('Central Frequency [Hz]')
ylabel('Band Level [dBV]')
title(sprintf('IIR Bandpass Filter Bank \\rm(%s)',signal))
set(gca,'XScale','log')
box on
axis([fc(1) fc(end) -180 0])
legend({'Without Zero Padding','With Zero Padding'},'Location','SouthEast')
fname = sprintf ('ANSI BPF (%s)',signal);
savefig(fname)
print(fname,'-dpng','-r200')
