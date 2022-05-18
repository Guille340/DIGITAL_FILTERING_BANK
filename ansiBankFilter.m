%  [xb,fc] = ansiBankFilter(AnsiFilterBank,x,varargin)
%
%  DESCRIPTION
%  Filters the input signal X with the filter bank ANSIFILTERBANK generated 
%  with ANSIBANKFILTERDESIGN. The function returns the filtered data XB and 
%  the central frequencies of the fractional octave bands FC. 
%
%  Property METRICSOUTPUT = FALSE assigns the bandpass-filtered signals to 
%  XB. The function can also return the band metrics (RMS, Exposure, Peak, 
%  Peak-to-Peak) calculated with the ANSI filter bank (METRICSOUTPUT = TRUE). 
%
%  The filtering function can be selected with FILTERMODE property (see 
%  FILTER and MYFILTFILT). 
%
%  An appropriate amount of zero-padding will also be applied to X when 
%  property ZEROPADDING = TRUE to achieve accurate filtering at low frequencies
%  on short-duration signals. 
%
%  In order to return a filtered signal XB with identical duration as X when 
%  zero-padding is applied, set property DATAWRAP = TRUE for the filtered 
%  signal to be wrapped using MATLAB's function DATAWRAP.
%
%  The progress of the filtering process can also be displayed by setting
%  property DISPLAYPROGRESS = TRUE.
% 
%  INPUT VARIABLES
%  - AnsiFilterBank: structure containing the filtering objects and 
%    information of a filter bank generated with ANSIBANKFILTERDESIGN.
%  - x: input signal. The sampling frequency must be identical to the one
%    used for the design of ANSIFILTERBANK. A difference in the sampling
%    rate will result in a shift in band central frequencies and the ANSI
%    standard not being met.
%  - PROPERTIES (varargin): the user can set up to five function properties. 
%    The properties must be given in (Name,Value) pairs, separated by comma 
%    with the property name first. The available property names and the 
%    description of their values are given below.
%    ¬ 'MetricsOutput': logical or numeric value [0,1] indicating whether the 
%      output variable XB shall return the waveforms or metrics (RMS, exposure,
%      peak, peak-to-peak) for each fractional octave band defined in 
%      ANSIFILTERBANK. If this property is omitted, the default METRICSOUTPUT 
%      = FALSE is used.
%      # FALSE (or 0): XB is a matrix of dimensions [LENGTH(X) LENGTH(FC)] 
%        that contains the bandpass-filtered signals, one per column. The 
%        length of XB is that of the zero-padded signal if ZEROPADDED = 
%        TRUE and DATAWRAP = FALSE. Otherwise, the length of XB will be the 
%        same as X. Slower than METRICSOUTPUT = TRUE, but gives the filtered 
%        waveforms for further processing. DEFAULT option.
%      # TRUE (or 1): XB is a matrix of dimensions [4, LENGTH(FC)] containing 
%        the signal metrics (RMS, exposure, peak, and peak-to-peak), one per 
%        band. The band metrics are calculated using the filters in the bank.
%    ¬ 'ZeroPadding': logical or [0,1] numeric value indicating whether zero-
%       padding will be applied to the input signal X. Zero-padding is
%       necessary in those cases where the signal is short enough that lowpass
%       or bandpass filtering will cause a delay that will displace the 
%       filtered signal beyond its original window length. This will result 
%       in a noticeable energy error. Depending on the signal duration and 
%       top cutoff frequency of the filter, the error can be as large as 
%       several tens of decibels. One-sided zero padding is applied when 
%       FILTERMODE = 'filter', and two-sided zero padding when using 
%       FILTERMODE = 'filtfilt'. Any filtering approach will suffer from this 
%       "time leakage" effect. If this property is omitted, ZEROPADDING = TRUE 
%       (default).
%     ¬ 'FilterMode': character string specifying the function to be used for
%       filtering. If this property is omitted, FILTERMODE = 'filter'.
%       # 'filter': fastest option (DEFAULT).
%       # 'filtfilt': zero-phase filtering option for IIR filters. Slower
%         than 'filter'. Use if the phase of the filtered signal is important 
%         and want its envelope to follow the envelope of the original signal 
%         (for example, for peak amplitude calculations). Note that custom 
%         function MYFILTFILT is used instead of MATLAB's FILTFILT due to large 
%         amplitude errors; these errors are caused by FILTFILT trying to 
%         minimise transient effects at the edges of the signal to be filtered.
%     ¬ 'DataWrap': logical or [0,1] numeric value indicating that wrapping
%        will be applied to the bandpass filtered signals in XB when choosing
%        the option 'METRICSOUTPUT' = FALSE. The function uses MATLAB's 
%        DATAWRAP to reduce the zero-padded signal ('ZEROPADDING' = TRUE) to 
%        the length of the original input signal X. The wrapping will result 
%        in a certain signal energy error due to amplitude cancellations, in 
%        particular at the lowest frequencies. Those minor errors are 
%        compensated by applying a correction factor to the wrapped signals.
%        If this property is omitted, DATAWRAP = TRUE (default). This property 
%        will be ignored when 'METRICSOUTPUT' = TRUE.
%     ¬ 'DisplayProgress': TRUE for the progress of the filter bank design to 
%        be displayed in the command window. If this property is omitted, the 
%        default DISPPROGRESS = FALSE is used.
%   
%  OUTPUT VARIABLES
%  - xb: matrix of filtered results with same units as X. For METRICSOUTPUT = 
%    FALSE, XB is a matrix of dimensions [LENGTH(X) LENGTH(FC)] containing 
%    the bandpass filtered signals, one per column. For METRICSOUTPUT = TRUE, 
%    XB is a matrix of dimensions [4, LENGTH(FC)] containing the RMS, exposure, 
%    peak, and peak-to-peak amplitudes for each fractional octave band in 
%    ANSIFILTERBANK.
%  - fc: central frequencies of the fractional octave bands in ANSIFILTERBANK
%    [Hz].
%
%  FUNCTION CALL
%  [xb,fc] = ANSIBANKFILTER(AnsiFilterBank,x)
%  [xb,fc] = ANSIBANKFILTER(...,PROPERTYNAME,PROPERTYVALUE)
%
%  If no properties are specified, their default values will be used. The list 
%  of default values is:
%
%   METRICSOUTPUT = FALSE;
%   ZEROPADDING = TRUE;
%   FILTERMODE = 'filter';
%   DATAWRAP = TRUE;
%   DISPROGRESS = FALSE; 
%
%  FUNCTION DEPENDENCIES
%  - ansiBankFilterDesign
%  - myfiltfilt
%
%  TOOLBOX DEPENDENCIES
%  - MATLAB (Core)
%
%  % EXAMPLE 1: "Calculate TOB Spectrum with FFT and IIR Filtering"
%  % 1) Configuration Data
%  fs = 44100;
%  bpo = 3;
%  freqLims = [20 10e3];
%  filtOrder = 10;
%  dispProgress = true;
%  T = 1; % signal duration [s]
%  x = 2*rand(1,T*fs) - 1; % signal
%  filtMode = 'filtfilt';
%
%  % 2) Filter Bank Design
%  AnsiFilterBank = ansiBankFilterDesign(fs,bpo,'FrequencyLimits',freqLims,...
%      'FilterOrder',filtOrder,'DisplayProgress',dispProgress)
%
%  % 3) Filter Signal with IIR and FFT Methods
%  [xb_iir,fc] = ansiBankFilter(AnsiFilterBank,x,'MetricsOutput',1,...
%      'FilterMode',filtMode);
%
%  % 4) Plot Band Level Spectrum
%  Xb_iir = 20*log10(xb_iir(1,:));
%  figure
%  hold on
%  plot(fc,Xb_iir,'b')
%  xlabel('Central Frequency [Hz]')
%  ylabel('Band Level [dBV]')
%  set(gca,'XScale','log')
%  box on
%
%  % EXAMPLE 2: "Calculate TOB Spectrum from IIR Banpass Filtered Signal"
%  % 1) Configuration Data 
%  [Identical to EXAMPLE 1]
%
%  % 2) Filter Bank Design 
%  [Identical to EXAMPLE 1]
%
%  % 3) Filter Signal with IIR Method (Waveform)
%  [xb_wav1,fc] = ansiBankFilter(AnsiFilterBank,x,'MetricsOutput',0,...
%      'FilterMode',filtMode,'DataWrap',false);
%  [xb_wav2,fc] = ansiBankFilter(AnsiFilterBank,x,'MetricsOutput',0,...
%      'FilterMode',filtMode,'DataWrap',true);
%
%  % 4) Calculate TOB RMS Amplitude
%  xb_iir1 = rms(xb_wav1)' * sqrt(size(xb_wav1,1)/length(x));
%  xb_iir2 = rms(xb_wav2)';  
%
%  % 5) Plot Band Level Spectrum
%  Xb_iir1 = 20*log10(xb_iir1);
%  Xb_iir2 = 20*log10(xb_iir2);
%  figure
%  hold on
%  plot(fc,Xb_iir1,'b')
%  plot(fc,Xb_iir2,'r')
%  xlabel('Central Frequency [Hz]')
%  ylabel('Band Level [dBV]')
%  legend({'Zero Pad','Data Wrap'},'Location','SouthEast')
%  set(gca,'XScale','log')
%  box on
%
%  REFERENCES
%  - ANSI (2004), "ANSI S1.11: Specification for Octave, Half-Octave and
%    Third Octave Band Filter Sets", American National Standards Institute
%  - BSI (2014), "Electroacoustics - Octave-band and fractional-octave-band
%    filters. Part 1: Specifications". European Standard BS EN 61260-1:2014
%    published by the British Standards Institution (BSI).
%  - Jiménez-Arranz, G. (2020). Filtering and Zero-Padding (Notes). Notes
%    prepared on 5 Jun 2020.
%
%  See also ANSIBANKFILTERDESIGN, ANSIFILTERMASK, ANSIBANKFILTERCOMPLIANCE, 
%  ANSIBANKFILTERPLOT, MYFILTFILT

%  VERSION 3.0
%  Date: 06 Aug 2021
%  Author: Guillermo Jimenez Arranz
%  - Renamed functions to inclulde the word 'bank'. This is done to avoid
%    conflict with the ANSI filtering toolbox with an individual filter.
%  - Removed the FFT filtering functionality. This will be implemented
%    separately for both bank and single filtering.
%
%  VERSION 2.0
%  Date: 16 May 2021
%  Author: Guillermo Jimenez Arranz
%  - Replaced TIC and TOC with NOW function due to potential conflict with
%    external TIC calls.
%  - Corrected the length of zero-padding (line 425).
%  - Added 'DataWrap' input option to reduce the zero-padded signal to the
%    length of the original input signal X.
%  - Included SEL, peak, and peak-to-peak metrics in XB when 'MetricsOutput'
%    1 or 2 are selected.
%  - Replaced FILTFILT function with MYFILTFILT due to problems associated
%    with the former that result in large transients at the start and end
%    of the signal.
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  10 Jun 2020

function [xb,fc] = ansiBankFilter(AnsiFilterBank,x,varargin)

% Check Number of Input Arguments
narginchk(2,12)
nVarargin = nargin - 2;
if rem(nVarargin,2)
    error('Property and value input arguments must come in pairs')
end

% Initialise Default Parameters
metricsOutput = 0; 
zeroPad = true;
filtMode = 'filter';
dataWrap = true;
dispProgress = false;

% Retrieve Input Variables
for m = 1:2:nVarargin
    inputProperty = lower(varargin{m}); % case insensitive
    inputProperties = lower({'MetricsOutput','ZeroPadding','FilterMode',...
        'DataWrap','DisplayProgress'});
    if ~ismember(inputProperty,inputProperties)
        error('Invalid input property')
    else
        switch inputProperty
            case 'metricsoutput'
                metricsOutput = varargin{m+1};
            case 'zeropadding'
                zeroPad = varargin{m+1};
            case 'filtermode'
                filtMode = varargin{m+1};
            case 'datawrap'
                dataWrap = varargin{m+1};
            case 'displayprogress'
                dispProgress = varargin{m+1};
        end
    end
end

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

% Error Control (check input signal)
if ~isreal(x) || ~isvector(x)
    error('The input signal must be a real numeric vector')
end

% Error Control (PROPERTY: 'MetricsOutput')
if all(metricsOutput ~= [0 1])
    metricsOutput = 0;
    warning(['Non-valid value for METRICSOUTPUT property. '...
        'METRICSOUTPUT = 0 will be used'])
end

% Error Control (PROPERTY: 'ZeroPadding')
if ~any(zeroPad == [0 1])
    zeroPad = true;
    warning(['ZEROPADDING property must be [0,1] or logical. '...
        'ZEROPADDING = TRUE will be used'])
end

% Error Control (PROPERTY: 'FilterMode')
switch filtMode
    case 'filter'
        filtFun = @filter;
    case 'filtfilt'
        filtFun = @myfiltfilt;
    otherwise
        filtFun = @filter;
end

% Error Control (PROPERTY: 'DataWrap')
if ~islogical(dataWrap) && all(dataWrap ~= [0 1])
    dataWrap = true;
    warning(['DATAWRAP property must be [0,1] or logical. '...
        'DATAWRAP = TRUE will be used'])
end

% Error Control (PROPERTY: 'DisplayProgress')
if ~any(dispProgress == [0 1])
    dispProgress = true;
    warning(['DISPLAYPROGRESS property must be [0,1] or logical. '...
        'DISPLAYPROGRESS = TRUE will be used'])
end

% Load Parameters
fs = AnsiFilterBank.sampleRate; % sampling rate of filter bank and input signal
decimBank = AnsiFilterBank.decimatorBank; % substructure of decimation filters
bpassBank = AnsiFilterBank.bandpassBank; % substructure of bandpass filters
bpo = bpassBank.bandsPerOctave; % band factor or no. bands per octave
fc = bpassBank.centralFreq; % central frequency of bandpass filters
bpassFilterBank = bpassBank.filterObject; % unique bandpass filter objects
ibpassFilter = bpassBank.filterObjectIndex; % indices of unique bandpass filters
ioctParent = bpassBank.parentOctaveIndex; % indices of parent octaves for each bandpass filter
decimFilter = decimBank.filterObject; % unique decimation filter object
decimFactor = decimBank.decimationFactor; % decimation factors for each parent octave

% Filter Input Signal
x = x(:); % convert x into column vector
L = length(x); % length of input signal [samples]
Lz = L; % intialise signal length to no zero-padding [samples]

% Zero-Padding
if zeroPad 
    nPeriods = 11.5*bpo;
    Lz = nPeriods * fs/min(fc) + 1230; % length of zero-padded signal (single-sided padding)
    Lz = 2^ceil(log2(Lz)); % round to highest power of 2
    Lz = max(Lz,L);
    nZeros = Lz - L; % number of zeros for single-sided zero padding
    if strcmp(filtMode,'filter')
        x = [x; zeros(nZeros,1)];
    else % filtfilt
        Lz = L + 2*nZeros; % length of zero-padded signal (double-sided padding)
        x = [zeros(nZeros,1); x; zeros(nZeros,1)]; % double-sided zero-padding
    end
end

% Apply Decimators and Bandpass Filters
nFraBands = length(fc); % number of fractional octave bands
nOctBands = length(decimFactor); % number of parect octave bands
if metricsOutput
    xb = nan(1,nFraBands);
else
    if dataWrap, xb = nan(L,nFraBands); 
    else, xb = nan(Lz,nFraBands); 
    end
end
D_prev = 1;
for m = 1:nOctBands
    if dispProgress, fprintf('Parent octave band %d/%d ',m,nOctBands); end

    % Decimate
    D = decimFactor(m);
    nDecimSteps = log2(D/D_prev);
    for n = 1:nDecimSteps % decimate X down to the current parent octave
        x = filtFun(decimFilter,x);
        x = downsample(x,2);
    end

    % Bandpass Filter
    iFraInOct = find(ioctParent == m); % index of fractional octave bands within current octave band (1...nFraBands)
    nFraInOct = length(iFraInOct); % number of fractional octave bands within current octave band (1...bpo)
    for n = 1:nFraInOct
        k = iFraInOct(n); % absolute index of current fractional octave band
        bpassFilter = bpassFilterBank(ibpassFilter(k));
        switch metricsOutput 
            case 0 % Calculate Filtered Waveform at Original Sampling Rate
                xbd = filtFun(bpassFilter,x); % bandpass-filtered signal
                nInterpSteps = log2(D); % number of interpolation steps
                xbi = xbd; % filtered signal to be interpolated to fs

                % Reconstruct Signals to Original Zero-Padding Length
                for p = 1:nInterpSteps
                    Li = ceil(Lz/2^(nInterpSteps-p)); % target size of current xbi vector
                    xbi = 2 * upsample(xbi,2); % upsample by 2
                    % Trim or Pad signal to Right Length
                    if strcmp(filtMode,'filter')
                        nZeros = max(Li - length(xbi),0);
                        nTrimSamples = max(length(xbi) - Li,0);
                        startSample = 1;
                        endSample = length(xbi) - nTrimSamples;
                        xbi = [xbi(startSample:endSample); zeros(nZeros,1)]; % interpolated signal trimmed to right size
                    else % filtfilt
                        nZeros = max(Li - length(xbi),0);
                        nZerosLeft = round(nZeros/2);
                        nZerosRight = nZeros - nZerosLeft;
                        nTrimSamples = max(length(xbi) - Li,0);
                        nTrimSamplesLeft = round(nTrimSamples/2);
                        nTrimSamplesRight = nTrimSamples - nTrimSamplesLeft;
                        startSample = 1 + nTrimSamplesLeft;
                        endSample = length(xbi) - nTrimSamplesRight;
                        xbi = [zeros(nZerosLeft,1); ...
                            xbi(startSample:endSample); ...
                            zeros(nZerosRight,1)]; % interpolated signal trimmed/padded to right size                            
                    end
                    % Decimate
                    xbi = filtFun(decimFilter,xbi); % decimated signal (anti-aliasing filtering)                      
                end 
                % Data Wrapping
                if dataWrap && zeroPad
                    xbw = datawrap(xbi,L); % wrap over original length
                    if strcmp(filtMode,'filtfilt')
                        nZeros = round((Lz - L)/2); % zero-padding samples (one-sided) for original signal
                        Lc = rem(nZeros,L); % number of samples for circular shift
                        xbw = circshift(xbw,-Lc);
                    end
                    weight = sqrt(sum(xbi.^2)/sum(xbw.^2)); % correction factor to compensate for any energy cancellation from wrapping
                    xbi = xbw*weight;
                end
                xb(:,k) = xbi; % matrix of bandpass-filtered signals
            case 1 % Calculate Band Metrics of Input Signal x (with Filter Bank)
                weight = sqrt(Lz/L); % correction factor to compensate the effect of zero-padding on the RMS
                xbd = filtFun(bpassFilter,x); % bandpass-filtered signal
                xb(1,k) = rms(xbd)*weight; % vector root-mean-square values
                xb(2,k) = rms(xbd)^2*Lz/fs; % vector exposure values (= sum(xbd.^2)*D/fs)
                xb(3,k) = max(abs(xbd)); % vector peak values
                xb(4,k) = max(xbd) - min(xbd); % vector peak-to-peak values   
        end
    end
    D_prev = D; % set previous decimation step to current one
    if dispProgress, fprintf('[%s]\n',datestr(now,'dd/mm/yyyy HH:MM:SS')); end
end
    
xb = fliplr(xb); % bandpass-filtered results (ascending order)
fc = fliplr(fc); % central frequencies of bandpass filters (ascending order)

