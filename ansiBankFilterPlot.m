%  hfig = ANSIBANKFILTERPLOT(AnsiFilterBank)
%
%  DESCRIPTION
%  Plots the decimators (red) and bank of bandpass filters (magenta) defined 
%  in the input structure ANSIFILTERBANK (see ANSIBANKFILTERDESIGN), along 
%  with the masks of the three classes defined in the ANSI S1.11 standard. 
%  The diamond markers on the decimator curves highlight the central frequency 
%  of their parent octaves and half-power points. The function returns the
%  handle to the generated figure.
% 
%  INPUT VARIABLES
%  - AnsiFilterBank: structure containing the filtering objects and information
%    of a filter bank generated with ANSIBANKFILTERDESIGN.
%   
%  OUTPUT VARIABLES
%  - hfig: figure handle
%
%  FUNCTION CALL
%  hfig = ansiBankFilterPlot(AnsiFilterBank)
%
%  FUNCTION DEPENDENCIES
%  - ansiFilterMask
%
%  TOOLBOX DEPENDENCIES
%  - MATLAB (Core)
%
%  EXAMPLE
%  % 1) Configuration Data
%  fs = 44100;
%  bpo = 3;
%  freqLims = [20 10e3];
%  filtOrder = 10;
%  dispProgress = true;
%
%  % 2) Filter Bank Design
%  AnsiFilterBank = ANSIBANKFILTERDESIGN(fs,bpo,'FrequencyLimits',freqLims,...
%      'FilterOrder',filtOrder,'DisplayProgress',dispProgress)
%
%  % 3) Plot Filter Bank
%  hfig = ANSIBANKFILTERPLOT(AnsiFilterBank)
%
%  REFERENCES
%  - ANSI (2004), "ANSI S1.11: Specification for Octave, Half-Octave and
%    Third Octave Band Filter Sets", American National Standards Institute
%  - BSI (2014), "Electroacoustics - Octave-band and fractional-octave-band
%    filters. Part 1: Specifications". European Standard BS EN 61260-1:2014
%    published by the British Standards Institution (BSI).
%
%  See also ANSIBANKFILTERDESIGN, ANSIFILTERMASK, ANSIBANKFILTERCOMPLIANCE, 
%  ANSIBANKFILTER

%  VERSION 2.0
%  Date: 06 Aug 2021
%  Author: Guillermo Jimenez Arranz
%  - Renamed functions to inclulde the word 'bank'. This is done to avoid
%    conflict with the ANSI filtering toolbox with an individual filter.
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  30 May 2020

function hfig = ansiBankFilterPlot(AnsiFilterBank)

% Error Control (check filter structure)
genfieldsValid = {'sampleRate','freqLimits','freqLimitMode','octaveRatio', ...
    'decimatorBank','bandpassBank'};
decfieldsValid = {'filterType','filterResponse','filterOrder',...
    'filterObject','halfPowerFreqn','decimationFactor','sampleRate',...
    'halfPowerFreq','parentOctaveIndex','parentOctaveCentralFreq',...
    'parentOctaveHalfPowerFreq1','parentOctaveHalfPowerFreq2'};
bpafieldsValid = {'filterType','filterResponse','bandsPerOctave',...
    'filterOrder','filterClass','filterObject','centralFreqn',...
    'halfPowerFreqn1','halfPowerFreqn2','filterObjectIndex',...
    'decimationFactor','sampleRate','nominalFreq','centralFreq',...
    'halfPowerFreq1','halfPowerFreq2','parentOctaveIndex',...
    'parentOctaveCentralFreq','parentOctaveHalfPowerFreq1',...
    'parentOctaveHalfPowerFreq2'};

genfields = fieldnames(AnsiFilterBank)';
if isequal(genfieldsValid,genfields)
    decfields = fieldnames(AnsiFilterBank.decimatorBank)';
    bpafields = fieldnames(AnsiFilterBank.bandpassBank)';
    if ~isequal(decfieldsValid,decfields) || ~isequal(bpafieldsValid,bpafields)
        error('The input argument ANSIFILTERBANK is not a valid filter structure')
    end
else
    error('The input argument ANSIFILTERBANK is not a valid filter structure')
end

% Load Parameters
decimBank = AnsiFilterBank.decimatorBank;
bpassBank = AnsiFilterBank.bandpassBank;
G = AnsiFilterBank.octaveRatio;
freqLims = AnsiFilterBank.freqLimits;
fs = AnsiFilterBank.sampleRate;
bpo = bpassBank.bandsPerOctave;
filtClass = bpassBank.filterClass;
filtOrder = bpassBank.filterOrder;
decimFilter = decimBank.filterObject;
decimFactor = decimBank.decimationFactor;
fdn1 = decimBank.halfPowerFreqn;
fbc = bpassBank.centralFreq;
fbr = bpassBank.sampleRate;
ibpassFilter = bpassBank.filterObjectIndex;
ioctParent = bpassBank.parentOctaveIndex;

% Compliance with ANSI Standard
if any(filtClass == -1)
    ansiComplStr = 'Non-Compliant with ANSI Standard'; 
else
    ansiComplStr = sprintf('Class %d',max(filtClass));
end

% Plot Filter Bank with ANSI Masks
hfig = [];
nOctBands = length(decimFactor); % number of parent octave bands
if nOctBands > 0
    hfig = figure;
    hold on
    nPointsInBand = 25; % plotting points per fractional octave band
    nFraOctBands = length(fbc);
    fmin = freqLims(1)*2^-3; % lower frequency limit for calculations
    fmax = fs/2; % upper frequency limit for calculations
    nPoints = round(nPointsInBand * log(fmax/fmin)/log(G^(1/bpo))); % total number of plotting points for each curve
    fAxis = fmin * (fmax/fmin).^((0:nPoints-1)/(nPoints-1)); % frequency axis
    nDecimators = sum(~isnan(decimBank.halfPowerFreq)); % number of decimation filters
    decimHandleIndex = nan(1,nDecimators);
    mark1HandleIndex = nan(1,nDecimators);
    mark2HandleIndex = nan(1,nDecimators);
    plotCounter = 0;
    iDecim = 0;
    for m = 1:nOctBands
        % Plot Decimation Filter for Current Octave Band
        D = decimFactor(m);    
        hDecimMark = [];
        if D >= 2
            fr = 2*fs/D; % sampling rate (before downsampling by 2)
            fAxis = fAxis(fAxis <= fr/2); % limit frequency axis to Nyquist frequency
            decimResp = 20*log10(abs(freqz(decimFilter,fAxis*2/fr,2))); % frequency response of decimator

            plot(fAxis,decimResp,'r','LineWidth',1); % plot decimation filter before decimation
            [~,iMark1] = min(abs(fAxis*2/fr - fdn1*2^-1.5)); % axis index for central frequency of parent octave band
            [~,iMark2] = min(abs(fAxis*2/fr - fdn1)); % axis index of half-power frequency of decimator
            hDecimMark = plot(fdn1*2^-1.5*fr/2,decimResp(iMark1),'r',...
                'LineWidth',1,'Marker','diamond','MarkerFaceColor','w',...
                'MarkerEdgeColor','r','MarkerSize',5);  % plot mark 1 (central frequency of octave)
            plot(fdn1*fr/2,decimResp(iMark2),'r','LineWidth',1,'Marker',...
                'diamond','MarkerFaceColor','w','MarkerEdgeColor','r',...
                'MarkerSize',5); % plot mark 2 (half-power frequency of decimator)
            plotCounter = plotCounter + 3; 
            iDecim = iDecim + 1;
            decimHandleIndex(iDecim) = plotCounter-2; % handle indices for decimator curve
            mark1HandleIndex(iDecim) = plotCounter-1; % handle indices for mark 1 (central freq of octave)
            mark2HandleIndex(iDecim) = plotCounter; % handle indices for mark 2 (hasl-power freq of decimator)
        end

        % Plot Bandpass Filters for Current Octave Band
        fr = fs/D; % sampling rate (after downsampling by 2)
        fAxis = fAxis(fAxis <= fr/2); % limit frequency axis to Nyquist frequency
        iFraInOct = find(ioctParent == m); % index of fractional octave bands within current octave band (1...nFraBands)
        nFraInOct = length(iFraInOct); % number of fractional octave bands within current octave band (1...bpo)
        for n = 1:nFraInOct
            k = iFraInOct(n); % absolute index of fractional octave band (1...nFraBands)
            bpassFilter = bpassBank.filterObject(ibpassFilter(k));
            fbnc = 2*fbc(k)/fr; % central frequency of fractional octave band

            [fmn,amin0,amax0] = ansiFilterMask(fbnc,bpo,0,0); % class 0 filter mask for current band
            [~,amin1,amax1] = ansiFilterMask(fbnc,bpo,1,0); % class 1 filter mask for current band
            [~,amin2,amax2] = ansiFilterMask(fbnc,bpo,2,0); % class 2 filter mask for current band
            bpassResp = 20*log10(abs(freqz(bpassFilter,fAxis*2/fr,2))); % current bandpass filter response

            hClass0 = plot(fmn * fr/2,-amin0,'g'); % plot class 0 mask (lower limit)
            plot(fmn * fbr(k)/2,-amax0,'g') % plot class 0 mask (upper limit)
            hClass1 = plot(fmn * fr/2,-amin1,'b'); % plot class 1 mask (lower limit)
            plot(fmn * fbr(k)/2,-amax1,'b') % plot class 1 mask (upper limit)
            hClass2 = plot(fmn * fr/2,-amin2,'k'); % plot class 2 mask (lower limit)
            plot(fmn * fbr(k)/2,-amax2,'k') % plot class 2 mask (upper limit)
            hBpass = plot(fAxis,bpassResp,'m','LineWidth',1.5); % plot bandpass filter response

            plotCounter = plotCounter + 7;
        end    
    end

    % Plot Settings
    set(gca,'XScale','log')
    title({'ANSI S1.11 Compliance of Filter Bank';...
        sprintf('\\rm\\fontsize{10}%d x 1/%d Oct Bands, IIR Order %d, %s',...
        nFraOctBands,bpo,filtOrder,ansiComplStr)})
    warning('off','MATLAB:legend:IgnoringExtraEntries')
    legend([hClass0 hClass1 hClass2 hBpass hDecimMark ],{'Class 0','Class 1',...
        'Class 2','Bandpass Filter','Decimator (fc,2.8fc)'},'Location','East')
    axis([fbc(end)*2^(-1/bpo) fs/2 -80 10])
    xlabel('Frequency [Hz]')
    ylabel('Filter Gain [dB]')
    box on
    grid on
    set(gcf,'Units','Normalized','OuterPosition',[0.05 0.07 0.9 0.90])

    % Place Decimator Curves on Top
    hplo = get(gca,'Children');
    nLines = length(hplo);
    decimHandleIndex = nLines - decimHandleIndex + 1;
    mark1HandleIndex = nLines - mark1HandleIndex + 1;
    mark2HandleIndex = nLines - mark2HandleIndex + 1;
    uistack([hplo(mark2HandleIndex); hplo(mark1HandleIndex); ...
        hplo(decimHandleIndex)],'top');
else
    warning('Input ANSIFILTERBANK is an empty filter bank structure')
end
