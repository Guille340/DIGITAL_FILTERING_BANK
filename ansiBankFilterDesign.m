%  AnsiFilterBank = ansiBankFilterDesign(fs,bpo,varargin)
%
%  DESCRIPTION
%  Designs a bank of digital IIR bandpass filters following ANSI S1.11 
%  standard. The function returns a structure ANSIFILTERBANK containing the 
%  filtering objects and information about the decimator and bandpass filters 
%  in the bank. ANSIFILTERBANK is used as input of ANSIBANKFILTER to extract 
%  the bandpass filtered signals or their metrics (RMS, Exposure, Peak, 
%  Peak-to-Peak) from a given broadband signal. 
% 
%  The filter bank and decimation filters can be plotted along with the ANSI 
%  masks with ANSIFILTERPLOT, allowing the user to inspect the response of 
%  the filters and their compliance with any ANSI class (0,1,2). For details 
%  about the content of the structure check the "OUTPUT VARIABLES" section.
% 
%  INPUT VARIABLES
%  - fs: sampling rate on which the filter bank is based. Note that FS needs 
%    to be identical to the sampling rate of signal to be filtered for the 
%    bandpass central frequencies to meet the standard (see ANSIBANKFILTER).
%  - bpo: bandwidth factor. Positive integer that specifies the number of 
%    bands in an octave (i.e. BPO = 1, octave; BPO = 2, half-octave; BPO = 
%    3, third-octave, etc). BPO determines the nominal bandwidth of the band-
%    pass filter. BPO is reciprocal of the bandwidth designator (see "Terms and
%    Definitions" in ANSI S1.11 or BS EN 61260-1:2014).
%  - PROPERTIES (varargin): the user can set up to four function properties. 
%    The properties must be given in (Name,Value) pairs, separated by comma 
%    with the property name first. The available property names and the 
%    description of their values are given below:
%    ¬ 'FrequencyLimits': two-element vector indicating the bottom and top
%      frequency limits [Hz]. The top frequency must be lower than the
%      Nyquist frequency (FREQLIMS(2) < FS/2). If this property is
%      omitted, the default FREQLIMS = [20 20e3] is used.
%    ¬ 'FilterOrder': order of the decimator and bandpass filters. Both
%      types of filters use an IIR Unconstrained Butterworth response. The 
%      order must be even and >= 2. Use FILTORDER >= 4 for reliable bandpass 
%      filtering (but bear in mind that only FILTORDER >= 8 meets the ANSI
%      standard). If this property is omitted, the default FILTORDER = 8 is
%      used. If this property is omitted, but the 'FilterClass' property is
%      provided, the filter order will be calculated based on the target
%      class (see description of 'FilterClass' property below). If both 
%      properties 'FilterOrder' and 'FilterClass', the first will have
%      priority over the second.
%    ¬ 'FilterClass': class of the fractional octave bandpass filter according 
%      to ANSI standard. There are three filter classes: 0, 1 and 2, the first 
%      being the most restrictive. FILTCLASS may adopt any of these four 
%      values: 0 (Class-0), 1 (Class-1), 2 (Class-2) or -1 (no class 
%      restrictions). When using the FILTCLASS property, ANSIBANKFILTERDESIGN
%      will find the minimum common order that will result in a class equal
%      to or lower than FILTCLASS for all bandpass filters in the bank.
%      Therefore, and unlike the 'FilterOrder' property, 'FilterClass' will 
%      require some additional calculations; these may take up to 50 seconds 
%      depending on the number of bandpass filters to be designed. If this
%      property is omitted, the default FILTCLASS = -1 is used and the default
%      FILTORDER = 8 is applied. If this property is omitted, but the 
%      'FilterOrder' property is provided, the provided filter order will
%      be used (see description of 'FilterOrder' property above).
%    ¬ 'DisplayProgress': TRUE for the progress of the filter bank design to be
%      displayed in the command window. If this property is omitted, the 
%      default DISPPROGRESS = FALSE is used.
%   
%  OUTPUT VARIABLES
%  - AnsiFilterBank: structure containing information and filtering objects for 
%    the decimators and bandpass filters in the bank. ANSIFILTERBANK comprises
%    the following fields and subfields:
%
%    ~~~ Level 1 Fields (AnsiFilterBank.<Level1Field>) ~~~
%    ¬ 'sampleRate': sampling rate used to design the filter bank [Hz] (see
%      input variable FS). Note that the signal to be filtered must have the 
%      same sampling rate if we want the central frequencies of the bandpass 
%      filters in the bank to meet the ANSI S1.11 (2004) and BS EN 61260-1 
%      (2014) standards.
%    ¬ 'freqLimits': two-element vector indicating the bottom and top frequency 
%      limits of the filter bank [Hz] (see 'FrequencyLimits' property).
%    ¬ 'freqLimitMode': reference point of first and last fractional band used 
%      to determine the bands falling within the specified frequency limits.
%      The reference can be the central frequency (fc), the bandedge frequency 
%      (fc * 2^(1/(2*bpo))) or the bandstop frequency (fc * 2^(1/bpo)). The 
%      function uses the 'bandedge' option, meaning that any band whose
%      bandegdge frequencies are within the specified frequency limits will
%      be included in the filter bank.
%    ¬ 'octaveRatio': number of times a frequency is higher than other for it
%      to be considered an octave of the first. The ANSI S1.11 (2004) and BS EN 
%      61260-1 (2014) standards give two options: base 10 (10^(3/10)) or
%      base 2 (2^1). A base-2 octave ratio is used for consistency with a
%      decimation factor multiple of 2.
%    ¬ decimatorBank: substructure containing the filtering object and specific
%      information for the unique decimator used in the filter bank design.
%    ¬ bandpassBank: substructure containing the filtering objects and specific
%      information for the fractional octave bandpass filters in the bank.
%
%    ~~~ Level 2 Fields (decimatorBank.<Level2Field>) ~~~
%    ¬ 'filterType': 'lowpass' for decimation filter.
%    ¬ 'filterResponse': 'IIR' or infinite impulse response (no FIR filters)
%    ¬ 'filterOrder': order of the decimator. The filter order is given 
%      specifically as an input of ANSIBANKFILTERDESIGN or calculated from the
%      target filter class. The same order is used for all filters in the
%      bank, including the decimator.
%    ¬ 'filterObject': digitalFilter class for an unconstrained, lowpass, IIR
%      Butterworth filter. This is the anti-aliasing filter to be applied
%      before downsampling. The filter is designed for its bandedge frequency
%      to be 1.5 times the central frequency of the octave band to be filtered.
%      Decimation by 2 is applied when the central frequency of the octave band
%      to be filtered is at least 2^2.5 times lower than the current sampling 
%      rate (i.e. sampling rate before downsampling by 2). This means that the 
%      normalised frequencies of the filter (and the filter itself) are the 
%      same for all decimation stages, therefore only one decimation filter is 
%      needed. A signal can be filtered with this object using the functions 
%      FILTER or FILTFILT (see ANSIBANKFILTER).
%    ¬ 'halfPowerFreqn': normalised half-power (-3 dB) frequency for the
%      decimation filter. The absolute frequency is calculated by
%      multiplying HALFPOPWERFREQN * FR/2, where FR is the sampling rate
%      for a given decimation stage (FR = 2*FS/D, before downsampling by 2).
%    ¬ 'decimationFactor': vector of factors by which the original sampling 
%      frequency FS and the number of samples in the original signal are 
%      reduced at each decimation stage (after downsampling by 2). As many 
%      elements as octave bands containing at least one fractional octave band.
%    ¬ 'sampleRate': sampling rate at each decimation stage before downsampling 
%      by 2 (FR = 2*FS/D). As many elements as octave bands containing at least
%      one fractional octave band.
%    ¬ 'halfPowerFreq': absolute half-power (-3 dB) frequency for the
%      decimation filter at each decimation stage. As many elements as octave
%      bands containing at least one fractional octave band. HALFPOWERFREQ
%      is NaN for the parent octave bands that do not have a decimation filter.
%    ¬ parentOctaveIndex: vector of indices of parent octave bands. As many
%      elements as octave bands with at least one fractional octave band.
%    ¬ parentOctaveCentralFreq: central frequency of parent octave bands [Hz]. 
%      As many elements as octave bands with at least one fractional band.
%    ¬ parentOctaveHalfPowerFreq1: low-edge frequency of parent octave bands 
%      [Hz]. As many elements as octave bands with at least one fractional band. 
%    ¬ parentOctaveHalfPowerFreq2: high-edge frequency of parent octave bands 
%      [Hz]. As many elements as octave bands with at least one fractional band. 
%    
%    ~~~ Level 2 Fields (bandpassBank.<Level2Field>) ~~~
%    ¬ filterType: 'bandpass' for the fractional octave filters.
%    ¬ filterResponse: 'IIR' or infinite impulse response (no FIR filters)
%    ¬ bandsPerOctave: bandwidth factor or number of fractional bands in an 
%      octave (see input variable BPO).
%    ¬ filterOrder: order of the unique bandpass filters. The filter order is 
%      given specifically as an input of ANSIBANKFILTERDESIGN or calculated 
%      from the target filter class.
%    ¬ filterClass: class of each of unique fractional octave bandpass filters.
%      Whilst the filter order is fixed for all bandpass filters, the class 
%      may vary.
%    ¬ filterObject: vector of unconstrained, bandpass, IIR Butterworth filters
%      of digitalFilter class. Each of these filters corresponds to the
%      unique filters in the bank, from highest to lowest frequency. BPO unique
%      filters that are reused after each decimation stage (note that the 
%      number of unique filters can be as high as 2*BPO, depending on how close
%      FMAX is to FS/2). A signal can be filtered with each of these filter 
%      objects using the functions FILTER or FILTFILT (see ANSIBANKFILTER).
%    ¬ centralFreqn: vector of normalised central frequencies for the unique
%      bandpass filters, in descending order. The absolute frequency is 
%      calculated by multiplying CENTRALFREQN * FR/2, where FR is the sampling 
%      rate for a given decimation stage (FR = FS/D, after downsampling by 2).
%    ¬ halfPowerFreqn1: vector of normalised low-edge half-power (-3 dB) 
%      frequencies of the unique bandpass filters, in descending order. The 
%      absolute frequency is calculated by multiplying HALFPOWERFREQN1 * FR/2, 
%      where FR is the sampling rate for a given decimation stage (FR = FS/D, 
%      after downsampling by 2).
%    ¬ halfPowerFreqn2: vector of normalised high-edge half-power (-3 dB) 
%      frequencies of the unique bandpass filters, in descending order. The 
%      absolute frequency is calculated by multiplying HALFPOWERFREQN2 * FR/2, 
%      where FR is the sampling rate for a given decimation stage (FR = FS/D, 
%      after downsampling by 2).
%    ¬ filterObjectIndex: vector with the indices of the unique filter object
%      (see 'filterObject') that each bandpass filter in the bank must use. 
%      As many elements as bandpass filters are in the bank.
%    ¬ decimationFactor: vector of factors by which the original sampling 
%      frequency FS and the number of samples in the original signal are 
%      reduced at each decimation stage (after downsampling by 2). As many 
%      elements as bandpass filters are in the bank.
%    ¬ sampleRate: sampling rate at each decimation stage after downsampling 
%      by 2 (FR = FS/D). As many elements as bandpass filters are in the bank.
%    ¬ centralFreq: vector of absolute central frequencies for the unique
%      bandpass filters, in descending order. As many elements as bandpass 
%      filters are in the bank.
%    ¬ halfPowerFreq1: vector of absolute low-edge half-power (-3 dB) 
%      frequencies for the unique bandpass filters, in descending order. As 
%      many elements as bandpass filters are in the bank.
%    ¬ halfPowerFreq2: vector of absolute high-edge half-power (-3 dB) 
%      frequencies for the unique bandpass filters, in descending order. As 
%      many elements as bandpass filters are in the bank.
%    ¬ parentOctaveIndex: vector of indices of parent octave bands. As many
%      elements as fractional octave bands.
%    ¬ parentOctaveCentralFreq: central frequency of parent octave bands [Hz]. 
%      As many elements as octave bands with at least one fractional band.
%    ¬ parentOctaveHalfPowerFreq1: low-edge frequency of parent octave bands 
%      [Hz]. As many elements as octave bands with at least one fractional band. 
%    ¬ parentOctaveHalfPowerFreq2: low-edge frequency of parent octave bands 
%      [Hz]. As many elements as octave bands with at least one fractional band. 
%
%  FUNCTION CALL
%  AnsiFilterBank = ANSIBANKFILTERDESIGN(fs,bpo)
%  AnsiFilterBank = ANSIBANKFILTERDESIGN(fs,bpo,PROPERTYNAME,PROPERTYVALUE,...)
%
%  If no properties are specified, their default values will be used. The list 
%  of default values is:
%
%   FREQUENCYLIMITS = [20 20e3];
%   FILTERORDER = 8;
%   FILTERCLASS = -1;
%   DISPLAYPROGRESS = FALSE; 
%
%  FUNCTION DEPENDENCIES
%  - fractionalOctaveBands
%  - ansiBankFilterCompliance
%  - ansiFilterMask
%
%  TOOLBOX DEPENDENCIES
%  - MATLAB (Core)
%
%  EXAMPLE
%  1) Configuration Data
%  fs = 44100;
%  bpo = 3;
%  freqLims = [20 10e3];
%  filtOrder = 10;
%  dispProgress = true;
%
%  2) Filter Bank Design
%  AnsiFilterBank = ANSIBANKFILTERDESIGN(fs,bpo,'FrequencyLimits',freqLims,...
%      'FilterOrder',filtOrder,'DisplayProgress',dispProgress)
%
%  REFERENCES
%  - ANSI (2004), "ANSI S1.11: Specification for Octave, Half-Octave and
%    Third Octave Band Filter Sets", American National Standards Institute
%  - BSI (2014), "Electroacoustics - Octave-band and fractional-octave-band
%    filters. Part 1: Specifications". European Standard BS EN 61260-1:2014
%    published by the British Standards Institution (BSI).
%
%  See also ANSIBANKFILTER, ANSIBANKFILTERPLOT, ANSIBANKFILTERMASK, 
%  ANSIBANKFILTERCOMPLIANCE

%  VERSION 2.0
%  Date: 06 Aug 2021
%  Author: Guillermo Jimenez Arranz
%  - Renamed functions to inclulde the word 'bank'. This is done to avoid
%    conflict with the ANSI filtering toolbox with an individual filter.
%
%  VERSION 1.1
%  Date: 16 May 2021
%  Author: Guillermo Jimenez Arranz
%  - Replaced TIC and TOC with NOW function due to potential conflict with
%    external TIC calls.
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  29 May 2020

function AnsiFilterBank = ansiBankFilterDesign(fs,bpo,varargin)

% Check Number of Input Arguments
narginchk(2,10)
nVarargin = nargin - 2;
if rem(nVarargin,2)
    error('Property and value input arguments must come in pairs')
end

% Initialise Default Parameters
freqLims = [20 20e3]; % frequency limits [Hz]
filtOrder = 8; % filter order
filtClass = -1; % filter class (-1,0,1,2)
dispProgress = false; % display progress flag
fmode = 'bandedge'; % filtering limits mode
base10 = 0; % fractional octave band factor (base10 = 0 -> base 2)
G = 2; % octave ratio (base 2)

% Retrieve Input Variables
for m = 1:2:nVarargin
    inputProperty = lower(varargin{m}); % case insensitive
    inputProperties = lower({'FrequencyLimits','FilterOrder',...
        'FilterClass','DisplayProgress'});
    if ~ismember(inputProperty,inputProperties)
        error('Invalid input property')
    else
        switch inputProperty
            case 'frequencylimits'
                freqLims = varargin{m+1};
            case 'filterorder'
                filtOrder = varargin{m+1};
            case 'filterclass'
                filtClass = varargin{m+1};
            case 'displayprogress'
                dispProgress = varargin{m+1};
        end
    end
end

% Error Control for Input Arguments
if numel(fs) > 1 || ~isreal(fs) || fs < 0
    error('FS must be a real positive number')
end
if numel(bpo) > 1 || ~isreal(bpo) || bpo < 1
    error('BPO must be a real number equal to 1 or higher')
end
if rem(bpo,1)
    bpo = round(bpo);
    warning('BPO must be a non-decimal number. B = %d will be used',bpo)
end
if numel(freqLims) == 2 && isreal(freqLims)
    fmin = min(freqLims); % minimum frequency of analysis [Hz]
    fmax = max(freqLims); % maximum frequency of analysis [Hz]
    if any(freqLims < 0)
        error('FREQUENCYLIMITS must be positive')
    end
    if fmax > fs/2
        fmax = fs/2;
        warning(['The top frequency limit must be lower than or equal to '...
            'half of the sample rate FS (FREQUENCYLIMITS(2) <= FS/2). '...
            'FREQUENCYLIMITS(2) = FS/2 will be used'])
    end
else
    error('FREQUENCYLIMITS must be a 2-element real vector')
end
if numel(filtOrder) > 1 || ~isreal(filtOrder)
    error('FILTORDER must be a real number')
end

if filtOrder < 2 || rem(filtOrder,2) || rem(filtOrder,1)
    filtOrder = round(filtOrder/2)*2;
    if filtOrder < 2, filtOrder = 8; end
    warning(['FILTORDER must be an even number greater than or equal to 2. '...
        'FILTORDER = %d will be used'],filtOrder)
end
inputProperties = varargin(1:2:nVarargin);
if all(ismember({'FilterOrder','FilterClass'},inputProperties))
    filtClass = -1;
    warning(['FILTERORDER has priority over FILTERCLASS. The latter will '...
        'be ignored'])
elseif ismember('FilterClass',inputProperties)
    if length(filtClass) > 1 || ~ismember(filtClass,[-1 0 1 2])
        filtClass = 2;
        warning(['FILTERCLASS must be 0, 1, 2 or -1. The default filter '...
            'class will be used (FILTERCLASS = %d)'],filtClass)
    end
end

% Central and Edge frequencies of Fractional Octave Bands
[f,fc_nom] = fractionalOctaveBands('BandsPerOctave',bpo,'FrequencyLimits',...
    [fmin fmax],'LimitMode',fmode,'Base10',base10); 
f = flipud(f); % sort frequencies in descending order
f1 = f(:,1)'; % low-edge frequencies for fractional octave bands
fc = f(:,2)'; % central frequencies for fractional octave bands
f2 = f(:,3)'; % high-edge frequencies for fractional octave bands
nFraBands = length(fc); % number of fractional octave bands

% Initialise General Parameters
fc_oct = [];
f1_oct = [];
f2_oct = [];

% Initialise Paramters of Decimation Filter
decimFiltType = 'lowpass'; % type of filter
decimFiltResponse = 'IIR'; % impulse response of filter
decimFiltOrder = 0; % filter order
decimFilter = []; % filter object
decimFn1 = []; % normalised half-power frequency
decimDecimFactor = []; % decimation factor (after downsampling by 2)
decimFd = fs; % sampling rate (before downsampling by 2)
decimF1 = []; % half-power frequencies
decimOctIndex = []; % indices of parent octave bands

% Initialise Parameters of Bandpass Filter
bpassFiltType = 'bandpass'; % type of filter
bpassFiltResponse = 'IIR'; % impulse response of filter
bpassFiltOrder = 0; % filter order
bpassFiltClass = -1; % filter class (-1,0,1,2)
bpassFilter = []; % filter object
bpassFnc = []; % normalised central frequency
bpassFn1 = []; % normalised low-edge half-power frequency
bpassFn2 = []; % normalised high-edge half-power frequency
bpassFilterIndex = []; % indices of unique bandpass filters
bpassDecimFactor = []; % decimation factors (after downsampling by 2)
bpassFd = fs; % sampling rate (after downsampling by 2)
bpassFc = []; % central frequencies
bpassF1 = []; % low-edge half-power frequencies
bpassF2 = []; % high-edge half-power frequencies
bpassOctIndex = []; % indices of parent octave bands

if nFraBands > 0
    % Central and Edge frequencies of Octave Bands
    [f_oct,~] = fractionalOctaveBands('BandsPerOctave',1,'FrequencyLimits',[fmin/2, fmax*2], ...
        'LimitMode',fmode,'Base10',base10); % central and edge frequencies of octave bands
    f_oct = flipud(f_oct); % sort frequencies in descending order
    f1_oct = f_oct(:,1)'; % low-edge frequencies for parent octave bands
    fc_oct = f_oct(:,2)'; % central frequencies for parent octave bands
    f2_oct = f_oct(:,3)'; % high-edge frequencies for parent octave bands

    % Calculate Index of Parent Octave Band for Each Fractional Octave Band
    bpassOctIndex = zeros(nFraBands,1);
    for ifra = 1:nFraBands
        bpassOctIndex(ifra) = find(fc(ifra) > f1_oct & fc(ifra) < f2_oct); % index of parent octave bands
    end
    bpassOctIndexUnique = unique(bpassOctIndex);
    fc_oct = fc_oct(bpassOctIndexUnique); % central frequency of octave bands containing at least one fractional octave band [Hz]
    f1_oct = f1_oct(bpassOctIndexUnique); % central frequency of octave bands containing at least one fractional octave band [Hz]
    f2_oct = f2_oct(bpassOctIndexUnique); % central frequency of octave bands containing at least one fractional octave band [Hz]
    bpassOctIndex = bpassOctIndex - bpassOctIndex(1) + 1; % start from 1 for correspondance with fc_oct
    nOctBands = length(fc_oct); % number of octave bands containing at least one fractional octave band

    % Calculate Number of Unique Bandpass Filters 
    bpassFnc = nan(1,nFraBands); % normalised central frequencies of fractional bandpass filters
    ioct = 1; % index of parent octave band (initialise)
    nUniqueFilters = 0; % number of unique fractional bandpass filters
    isFilterUnique = true; 
    while isFilterUnique && ioct <= nOctBands
        fd1 = fc_oct(ioct) * G^1.5; % half-power frequency of current decimator
        frMin = min([fs, 2*fd1]); % minimum sampling rate to avoid aliasing
        D = 2^floor(log2(fs/frMin)); % decimation factor
        fr = fs/D; % sampling rate (after downsampling)

        iFraInOct = find(bpassOctIndex == ioct); % index of fractional octave bands within current octave band (1...nFraBands)
        nFraInOct = length(iFraInOct); % number of fractional octave bands within current octave band (1...bpo)
        for n = 1:nFraInOct
            ifra = iFraInOct(n); % absolute index of fractional octave band within current octave band (1...nFraBands)
            fbnc = 2 * fc(ifra)/fr; % normalised central frequency of current fractional octave band
            isFilterUnique = ~ismember(round(fbnc,14),round(bpassFnc,14)); % unique filter = unique normalised central freq.
            if isFilterUnique, nUniqueFilters = nUniqueFilters + 1; end % unique filter count
            bpassFnc(ifra) = fbnc; % normalised central frequencies of fractional bandpass filters
        end
        ioct = ioct + 1; % index of octave band
    end
    clear bpassFnc 

    % Find Common Filter Order to Meet Class Specification for ALL Bandpass Filters
    if filtClass ~= -1
        filtOrderFilterBank = nan(1,nFraBands); % vector of filter orders for all fractional octave bandpass filters
        bpassFnc = nan(1,nFraBands); 
        for ioct = 1:nOctBands
            fd1 = fc_oct(ioct) * G^1.5; % half-power frequency of current decimator
            frMin = min([fs, 2*fd1]); % minimum sampling rate to avoid aliasing
            D = 2^floor(log2(fs/frMin)); % decimation factor for current octave (after downsampling by 2)
            fr = fs/D; % decimated sampling rate (after downsampling by 2)

            iFraInOct = find(bpassOctIndex == ioct); % absolute index of fractional octave bands within current octave band (1...nFraBands)
            nFraInOct = length(iFraInOct); % number of fractional octave bands within current octave band (1...bpo)
            for n = 1:nFraInOct
                ifra = iFraInOct(n); % absolute index of fractional octave band within current octave band (1...nFraBands)
                fbnc = 2 * fc(ifra)/fr; % normalised low-edge frequency for current fractional octave band
                fbn1 = 2 * f1(ifra)/fr; % normalised central frequency for current fractional octave band
                fbn2 = 2 * f2(ifra)/fr; % normalised high-edge frequency for current fractional octave band
                isFilterUnique = ifra <= nUniqueFilters; % unique bandpass filters are consecutive
                if isFilterUnique % design bandpass filters
                    if dispProgress, fprintf(['Calculating common filter order '...
                            '%d/%d '],ifra,nUniqueFilters); end
                    filtClassForBand = 3; % class for current fractional octave bandpass filter
                    filtOrderForBand = 2; % order for current fractional octave bandpass filter
                    while filtClassForBand == -1 || filtClassForBand > filtClass
                        filtBandpass = designfilt('bandpassiir','FilterOrder',filtOrderForBand,...
                            'HalfPowerFrequency1',fbn1,'HalfPowerFrequency2',fbn2); % desing bandpass filter
                        filtClassForBand = ansiBankFilterCompliance(fbnc,bpo,filtBandpass); % calculate filter class
                        filtOrderForBand = filtOrderForBand + 2; % test even filter orders from 2
                    end
                    filtOrderFilterBank(ifra) = filtOrderForBand; % vector of order for filters in bank
                    if dispProgress, fprintf('[%s]\n',datestr(now,'dd/mm/yyyy HH:MM:SS')); end
                end
                bpassFnc(ifra) = fbnc;
            end
        end
        filtOrder = max(filtOrderFilterBank); % common filter order to meet class condition in all bandpass filters
    end
    
    % Design Decimator (Anti-Alias Filter)
    if dispProgress, fprintf('Designing decimator '), end
    decimOctIndex = 1:nOctBands;
    ioct = find(fc_oct <= fs/2 * 2^-2.5,1,'first'); % index of first octave with decimation filter
    if ~isempty(ioct) % if decimation filter exists
        fd1 = fc_oct(ioct)*2^1.5; % half-power frequency of current decimator
        frMin = 2*fd1; % minimum sampling rate to avoid aliasing
        D = 2^floor(log2(fs/frMin)); % decimation factor (after downsampling by 2)
        fr = 2*fs/D; % sampling rate (resampled) (before downsampling by 2)
        decimFn1 = 2*fd1/fr; % normalised half-power frequency (before downsampling by 2)
        decimFiltOrder = filtOrder;
        decimFilter = designfilt('lowpassiir','FilterOrder',decimFiltOrder,...
                    'HalfPowerFrequency',decimFn1); % design decimator (anti-alias filter)
    end
    if dispProgress, fprintf('[%s]\n',datestr(now,'dd/mm/yyyy HH:MM:SS')), end

    % Design Bandpass Filters
    bpassFiltClass = nan(1,nUniqueFilters); % class of designed bandpass filters
    bpassFnc = nan(1,nUniqueFilters); % normalised central frequency of designed bandpass filters
    bpassFn1 = nan(1,nUniqueFilters); % normlalised low-edge frequency of designed bandpass filters
    bpassFn2= nan(1,nUniqueFilters); % normlalised high-edge frequency of designed bandpass filters
    bpassFilterIndex = nan(1,nFraBands); % index of designed filter for each fractional octave band
    bpassFd = nan(1,nFraBands); % decimated sampling frequency for each fractional octave band 
    bpassDecimFactor = nan(1,nFraBands); % decimation factor for each fractional octave band
    bpassFc_nom = nan(1,nFraBands); % nominal central frequency for each fractional octave band
    bpassFc = nan(1,nFraBands); % central frequency for each fractional octave band
    bpassF1 = nan(1,nFraBands); % low-edge frequency for each fractional octave band
    bpassF2 = nan(1,nFraBands); % high-edge frequency for each fractional octave band
    decimFd = nan(1,nOctBands); % decimated sampling frequency for each parent octave band
    decimDecimFactor = nan(1,nOctBands); % decimation factor for each octave band
    decimF1 = nan(1,nOctBands); % edge frequency of decimator for each octave band
    bpassFiltOrder = filtOrder;
    clear bpassFilter % clear variable to avoid class conflict
    for m = 1:nOctBands
        fd1 = fc_oct(m) * G^1.5; % half-power frequency of current decimator
        frMin = min([fs, 2*fd1]); % minimum sampling rate to avoid aliasing
        D = 2^floor(log2(fs/frMin)); % decimation factor for current octave (after downsampling)
        fr = fs/D; % decimated sampling rate (after downsampling)
        if D == 1, fd1 = NaN; end % no decimation filter for unity decimation factor

        % Parameters of Decimation Filterbank
        decimF1(m) = fd1; % absolute half-power frequency of decimator
        decimDecimFactor(m) = D; % decimation factor (for downsampling after LPF)
        decimFd(m) = min([fs, 2*fr]); % sampling rate (before downsampling by 2)

        % Design Bandpass Filters
        iFraInOct = find(bpassOctIndex == m); % index of fractional octave bands that are within current octave band
        nFraInOct = length(iFraInOct); % number of fractional octave bands within current octave band
        for n = 1:nFraInOct
            ifra = iFraInOct(n); % absolute index of fractional octave band within current octave band (1...nFraBands)
            fbnc = 2 * fc(ifra)/fr; % normalised central frequency of current bandpass filter
            fbn1 = 2 * f1(ifra)/fr; % normalised low-edge frequency of current bandpass filter
            fbn2 = 2 * f2(ifra)/fr; % normalised high-edge frequency of current bandpass filter
            isFilterUnique = ifra <= nUniqueFilters; % unique bandpass filters are consecutive
            if isFilterUnique % design bandpass filters
                if dispProgress, fprintf(['Designing bandpass filter '...
                        '%d/%d '],ifra,nUniqueFilters); end
                filtBandpass = designfilt('bandpassiir','FilterOrder',...
                    bpassFiltOrder,'HalfPowerFrequency1',fbn1,...
                    'HalfPowerFrequency2',fbn2); 
                filtClassForBand = ansiBankFilterCompliance(fbnc,bpo,filtBandpass);
                bpassFilter(ifra) = filtBandpass;  %#ok<AGROW>

                % Store Bandpass Filter Specs in ¦filterStructure¦
                bpassFiltClass(ifra) = filtClassForBand; 
                bpassFnc(ifra) = fbnc; % vector of normalised central frequences for bandpass filters
                bpassFn1(ifra) = fbn1; % vector of normalised low-edge frequences for bandpass filters
                bpassFn2(ifra) = fbn2; % vector of normalised high-edge frequences for bandpass filters

                if dispProgress, fprintf('[%s]\n',datestr(now,'dd/mm/yyyy HH:MM:SS')); end
            end

            bpassFilterIndex(ifra) = find(round(fbnc,14) == ...
                round(bpassFnc,14),1,'first');
            bpassDecimFactor(ifra) = D;
            bpassFc_nom(ifra) = fc_nom(ifra); % vector of nominal central frequencies for bandpass filters
            bpassFc(ifra) = fc(ifra); % vector of central frequences for bandpass filters
            bpassF1(ifra) = f1(ifra); % vector of low-edge frequences for bandpass filters
            bpassF2(ifra) = f2(ifra); % vector of high-edge frequences for bandpass filters
            bpassFd(ifra) = fr; % vector of sampling rates (after downsampling by 2)
        end    
    end
end

% Store General Parameters in ¦filterStructure¦
AnsiFilterBank.sampleRate = fs; % original sampling frequency [Hz]
AnsiFilterBank.freqLimits = [fmin fmax]; % frequency range [fmin fmax]
AnsiFilterBank.freqLimitMode = fmode; % calculation method for fractional octave bands (fmode = 'bandedge')
AnsiFilterBank.octaveRatio = G; % octave ratio

% Store Decimator Filterbank in ¦filterStructure¦
AnsiFilterBank.decimatorBank.filterType = decimFiltType;
AnsiFilterBank.decimatorBank.filterResponse = decimFiltResponse;
AnsiFilterBank.decimatorBank.filterOrder = decimFiltOrder;
AnsiFilterBank.decimatorBank.filterObject = decimFilter;
AnsiFilterBank.decimatorBank.halfPowerFreqn = decimFn1;
AnsiFilterBank.decimatorBank.decimationFactor = decimDecimFactor;
AnsiFilterBank.decimatorBank.sampleRate = decimFd;
AnsiFilterBank.decimatorBank.halfPowerFreq = decimF1;
AnsiFilterBank.decimatorBank.parentOctaveIndex = decimOctIndex;
AnsiFilterBank.decimatorBank.parentOctaveCentralFreq = fc_oct;
AnsiFilterBank.decimatorBank.parentOctaveHalfPowerFreq1 = f1_oct;
AnsiFilterBank.decimatorBank.parentOctaveHalfPowerFreq2 = f2_oct;

% Store Bandpass Filterbank in ¦filterStructure¦
AnsiFilterBank.bandpassBank.filterType = bpassFiltType;
AnsiFilterBank.bandpassBank.filterResponse = bpassFiltResponse;
AnsiFilterBank.bandpassBank.bandsPerOctave = bpo;
AnsiFilterBank.bandpassBank.filterOrder = bpassFiltOrder;
AnsiFilterBank.bandpassBank.filterClass = bpassFiltClass;
AnsiFilterBank.bandpassBank.filterObject = bpassFilter;
AnsiFilterBank.bandpassBank.centralFreqn = bpassFnc;
AnsiFilterBank.bandpassBank.halfPowerFreqn1 = bpassFn1;
AnsiFilterBank.bandpassBank.halfPowerFreqn2 = bpassFn2;
AnsiFilterBank.bandpassBank.filterObjectIndex = bpassFilterIndex;
AnsiFilterBank.bandpassBank.decimationFactor = bpassDecimFactor;
AnsiFilterBank.bandpassBank.sampleRate = bpassFd;
AnsiFilterBank.bandpassBank.nominalFreq = bpassFc_nom;
AnsiFilterBank.bandpassBank.centralFreq = bpassFc;
AnsiFilterBank.bandpassBank.halfPowerFreq1 = bpassF1;
AnsiFilterBank.bandpassBank.halfPowerFreq2 = bpassF2;
AnsiFilterBank.bandpassBank.parentOctaveIndex = bpassOctIndex';
AnsiFilterBank.bandpassBank.parentOctaveCentralFreq = fc_oct;
AnsiFilterBank.bandpassBank.parentOctaveHalfPowerFreq1 = f1_oct;
AnsiFilterBank.bandpassBank.parentOctaveHalfPowerFreq2 = f2_oct;
