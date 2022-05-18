%  [fb,amin,amax] = ANSIFILTERMASK(fc,bpo,filtClass,base10)
%
%  DESCRIPTION
%  Calculates the maximum and minimum attenuation limits for an ANSI S1.11 
%  bandpass filter of bandwidth 1/BPO octaves, central frequency FC, and class
%  FILTCLASS.
% 
%  INPUT VARIABLES
%  - fc: central frequency of bandpass filter. The frequency should meet the 
%    ANSI S1.11 standard (see FRACTIONALOCTAVEBANDS), but the attenuation mask 
%    will be calculated regardless of which value is selected for FC.
%  - bpo: bandwidth factor. Positive integer that specifies the number of 
%    bands in an octave (i.e. BPO = 1, octave; BPO = 2, half-octave; BPO = 3, 
%    third-octave, etc). BPO determines the nominal bandwidth of the band-pass
%    filter. BPO is reciprocal of the bandwidth designator (see "Terms and
%    Definitions" in ANSI S1.11 or BS EN 61260-1:2014).
%  - filtClass: class of the fractional octave bandpass filter according to 
%    ANSI standard. There are three filter classes: 0, 1 and 2, the first being
%    the most restrictive. FILTCLASS may adopt any of these four values: 0 
%    (Class-0), 1 (Class-1), 2 (Class-2) or -1 (no class restrictions).
%  - base10: true for a base-10 octave ratio. The octave ratio is the number of
%    times a frequency is higher than other for it to be considered an octave 
%    of the first. The ANSI S1.11 (2004) and BS EN 61260-1 (2014) standards 
%    give two options: base-10 and base 2 (2^1).
%    ¬ FALSE: base 2 octave ratio
%    ¬ TRUE: base 10 octave ratio
%   
%  OUTPUT VARIABLES
%  - fb: breakpoint frequencies of the attenuation mask for bandpass filter
%    of central frequency FC [Hz]
%  - amin: relative attenuation at FB for minimum limit mask (positive) [dB]
%  - amax: relative attenuation at FB for maximum limit mask (positive) [dB]
%
%  FUNCTION CALL
%  [fb,amin,amax] = ANSIFILTERMASK()
%  [fb,amin,amax] = ANSIFILTERMASK(fc)
%  [fb,amin,amax] = ANSIFILTERMASK(fc,bpo)
%  [fb,amin,amax] = ANSIFILTERMASK(fc,bpo,filtClass)
%  [fb,amin,amax] = ANSIFILTERMASK(fc,bpo,filtClass,base10)
%
%  FC, B, FILTCLASS, or BASE10 can be left empty ([]). Whenever any of these
%  input arguments are left empty or omitted from the call, their default 
%  values will be used. The list of default values is:
%
%  FC = 1000
%  B = 1
%  FILTCLASS = 2
%  BASE10 = FALSE
%
%  FUNCTION DEPENDENCIES
%  - None
%
%  TOOLBOX DEPENDENCIES
%  - MATLAB (Core)
%
%  EXAMPLE
%  % 1) Configuration Data
%  fs = 1000;
%  fc = 250;
%  bpo = 3;
%  filtClass = 0;
%  filtOrder = 10;
%  nPoints = 1000;
%
%  % 2) Class-0 ANSI S1.11 Mask and Bandpass Butterworth Filter Order 10
%  fcn = 2*fc/fs;
%  fn1 = fcn * 2^(-1/(2*bpo)); % normalised bottom half-power frequency
%  fn2 = fcn * 2^(1/(2*bpo)); % normalised top half-power frequency
%  [fbn,amin,amax] = ansiFilterMask(fcn,bpo,filtClass);
%  bandpassFilter = designfilt('bandpassiir','FilterOrder',filtOrder,...
%  'HalfPowerFrequency1',fn1,'HalfPowerFrequency2',fn2);
%  fbni = (0:nPoints-1)*(fbn(end)-fbn(1))/(nPoints-1) + fbn(1);
%  fbni = fbni(fbni <= 1);
%  hb = 20*log10(abs(freqz(bandpassFilter,fbni,2)));
%
%  % 3) Plot ANSI Mask and Filter
%  figure
%  hold on
%  hmask = plot(fbn,-amin,'b');
%  plot(fbn,-amax,'b');
%  hfilt = plot(fbni,hb,'m','LineWidth',1.5); 
%  axis([0 1 -80 10])
%  xlabel('Frequency [*2/fs]')
%  ylabel('Gain [dB]')
%  title({'Class-0 ANSI S1.11 Mask for a Third-Octave Bandpass Filter';...
%        '\rm\fontsize{10}Central Frequency = 250 Hz, Filter Order = 10'})
%  legend([hmask,hfilt],{'Class-0 ANSI S1.11','Bandpass Filter'})
%  box on
%
%  REFERENCES
%  - ANSI (2004), "ANSI S1.11: Specification for Octave, Half-Octave and
%    Third Octave Band Filter Sets", American National Standards Institute
%
%  See also ANSIBANKFILTERCOMPLIANCE, ANSIBANKFILTERDESIGN, ANSIBANKFILTER, 
%  ANSIBANKFILTERPLOT 

%  VERSION 2.0
%  Date: 06 Aug 2021
%  Author: Guillermo Jimenez Arranz
%  - Renamed functions to inclulde the word 'bank'. This is done to avoid
%    conflict with the ANSI filtering toolbox with an individual filter.
%
%  VERSION 1.1
%  Date: 22 May 2020
%  Author: Guillermo Jimenez Arranz
%  - Added error control for input arguments
%  - Additional input argument BASE10 to allow the user to specify the base
%    for the octave ratio (before, BASE10 = 1 was the default).
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  19 May 2019

function [fb,amin,amax] = ansiFilterMask(fc,bpo,filtClass,base10)

% Error Control (input argument check)
if nargin < 1 || isempty(fc)
    fc = 1000; % ref frequency of 1 kHz as band central frequency
end
if nargin < 2 || isempty(bpo)
    bpo = 1; % one octave bandwidth
end
if nargin < 3 || isempty(filtClass)
    filtClass = 2; % Class-2 ANSI mask
end
if nargin < 4 || isempty(base10)
    base10 = 0; % ovatave ratio with base 2
end

% Octave Frequency Ratio G
if base10
    G = 10^(3/10); % base-10 octave ratio
else
    G = 2; % base-2 octave ratio
end 

% Frequency Breakpoints
fn = [G^(1/8), G^(1/4), G^(3/8), G^(1/2) - eps, G^(1/2), G^(1/2) + eps,...
    G, G^2, G^3, G^4, G^5];
fnbHi = 1 + ((G^(1/(2*bpo)) - 1)/(G^(1/2) - 1)) * (fn - 1); % normalised breakpoint frequencies above central
fnbLo = 1./(1 + ((G^(1/(2*bpo)) - 1)/(G^(1/2) - 1)) * (fn - 1)); % normalised breakpoint frequencies below central
fnb = [fliplr(fnbLo), 1, fnbHi]; % normalised breakpoint frequencies
fb = fnb * fc; % absolute breakpoint frequencies

% Relative Attenuation for Minimum and Maximum Limits
switch filtClass
    case 0
        amax = [-0.15, -0.15, -0.15, -0.15, 2.30, 2.30, 18.0, 42.5, 62.0, 75.0, 75.0];
        amin = [ 0.20,  0.40,  1.10,  4.50, 4.50,  200,  Inf,  Inf,  Inf,  Inf,  Inf];
        amax = [fliplr(amax), -0.15, amax]; 
        amin = [fliplr(amin),  0.15, amin]; 
    case 1
        amax = [-0.30, -0.30, -0.30, -0.30, 2.00, 2.00, 17.5, 42.0, 61.0, 70.0, 70.0];
        amin = [ 0.40,  0.60,  1.30,  5.00, 5.00,  200,  Inf,  Inf,  Inf,  Inf,  Inf];
        amax = [fliplr(amax), -0.30, amax]; 
        amin = [fliplr(amin),  0.30, amin]; 
    case 2
        amax = [-0.50, -0.50, -0.50, -0.50, 1.60, 1.60, 16.5, 41.0, 55.0, 60.0, 60.0];
        amin = [ 0.60,  0.80,  1.60,  5.50, 5.50,  200,  Inf,  Inf,  Inf,  Inf,  Inf];
        amax = [fliplr(amax), -0.50, amax]; 
        amin = [fliplr(amin),  0.50, amin]; 
end

