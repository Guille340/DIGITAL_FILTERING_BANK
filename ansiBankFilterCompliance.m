%  ansiClass = ANSIBANKFILTERCOMPLIANCE(fn,bpo,varargin)
%
%  DESCRIPTION
%  Returns the ANSI S1.11 class that the input bandpass filter is compliant 
%  with (0,1,2). The function returns -1 if the input filter does not meet the
%  class 2 specification, the least restrictive.
%
%  INPUT VARIABLES
%  - fn: normalised central frequency of bandpass filter (FN = 2*FC/FS) 
%    [Hz]. FN is a decimal value between 0 (= 0 Hz) and 1 (= FS/2). The 
%    code does not verify if FN is an ANSI S1.11 midband frequency. FN 
%    must be set to the normalised central frequency used to design the 
%    input bandpass filter for the compliance test to work.
%  - bpo: bandwidth factor or number of bands per octave (i.e. b = 1, 
%    octave; b = 2, half-octave; b = 3, third-octave, etc). Reciprocal 
%    of the bandwidth designator (see "Terms and Definitions" in BS EN 
%    61260-1:2014 for definition).
%  - b,a (varargin{1},varargin{2}): coefficients of the input filter. For
%    a FIR filter, A = 1.
%  - filtobj (varargin{1}): filter object of class type 'dfilt.df2sos' 
%    (see fdesign.octave) or 'digitalFilter' (see designfilt).
%
%  OUTPUT VARIABLES
%  - ansiClass: integer indicating the ANSI S1.11 class that the input 
%    bandpass filter is compliant with. Four possible values:
%    ¬  0: class 0 (the most demanding filter response).
%    ¬  1: class 1
%    ¬  2: class 2 (the least demanding filter response).
%    ¬ -1: filter not compliant with any class in ANSI S1.11.
%
%  INTERNALLY CALLED FUNCTIONS
%  - ansiFilterMask
%
%  FUNCTION CALLS
%  - ansiClass = ANSIBANKFILTERCOMPLIANCE(fn,bpo,b,a)
%    ¬ for FIR or IIR filters with coefficients (b,a)
%  - ansiClass = ANSIBANKFILTERCOMPLIANCE(fn,bpo,filtobj)
%    ¬ for filters designed with fdesign.octave (class 'dfilt.df2sos')
%      or designfilt ('digitalFilter').
%
%  See also ANSIFILTERMASK, ANSIBANKFILTERDESIGN, ANSIBANKFILTERPLOT, 
%  ANSIBANKFILTER

%  VERSION 1.1
%  Date: 30 May 2020 
%  Author: Guillermo Jimenez Arranz
%  - Updated help
%  - Normalised frequencies now use FS = 2 rather than FS = 2*PI
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  20 May 2019

function ansiClass = ansiBankFilterCompliance(fn,bpo,varargin)

switch nargin
    case {0 1 2}
        error('Not enough input arguments')
    case 3
        filtobj = varargin{1}; % filter object
        
        % Error Control
        strclass = class(filtobj); % class type of filter object
        if ~any(strcmp(strclass,{'dfilt.df2sos','digitalFilter'}))
            error('The first input argument is not a valid filter class')
        end 
        
        % Class Compliance Analysis
        ansiClass = -1; % initialize filter compliance vector
        for filtClass = 0:2
            if ansiClass == -1 % if filter does not comply with previous class ...
                [fb,amin,amax] = ansiFilterMask(fn,bpo,filtClass,1);
                ival = find(fb <= 1); % index of breakpoints below Nyquist frequency
                H = 20*log10(abs(freqz(filtobj,fb,2))); % magnitude response of bandpass filter [dB]
                Hmin = -amin; % magnitude response of mask's upper limit for ANSI S1.11 class filtClass
                Hmax = -amax; % magnitude response of mask's lower limit for ANSI S1.11 class filtClass
                if all((H(ival) >= Hmin(ival)) & (H(ival) <= Hmax(ival))) % if filter complies with current class
                    ansiClass = filtClass;
                end
            end
        end
        
    case 4
        b = varargin{1}; % numerator filter coefficients
        a = varargin{2}; % denominator filter coefficients
        
        % Error Control
        if ~(isnumeric(a) && isnumeric(b))
            error('First and second input arguments must be numeric vectors')
        end
        
        % Class Compliance Analysis
        ansiClass = -1; % initialize filter compliance vector
        for filtClass = 0:2
            if ansiClass == -1 % if filter does not comply with previous class ...
                [fb,amin,amax] = ansiFilterMask(fn,bpo,filtClass,1);
                ival = find(fb <= 1); % index of breakpoints below Nyquist frequency
                H = 20*log10(abs(freqz(b,a,fb,2))); % magnitude response of bandpass filter [dB]
                Hmin = -amin; % magnitude response of mask's upper limit for ANSI S1.11 class filtClass
                Hmax = -amax; % magnitude response of mask's lower limit for ANSI S1.11 class filtClass
                if all((H(ival) >= Hmin(ival)) & (H(ival) <= Hmax(ival))) % if filter complies with current class
                    ansiClass = filtClass;
                end
            end
        end
           
    otherwise
        error('Invalid number of input arguments')
end
