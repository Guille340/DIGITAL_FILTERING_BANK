Digital Filtering (Bank, ANSI)
===============================

MATLAB code for designing and applying a bank of bandpass filters.

The filter bank meets the ANSI S1.11 standard. 

The filtering function can either return the set of bandpass filtered signals 
or the metrics per band (SPLrms, SEL, SPLpk, SPLpp). 

Zero-padding is applied when the duration of the signal is comparable to or 
smaller than the period of the lowest frequency to be filtered, to avoid issues 
associated with delay of filters. 

Zero-phase filtering can be applied to keep the phase information of the signal 
intact. This can be useful for estimating the peak amplitude of filtered bands.

The further details on the code, useful articles, and ANSI standard check the
"Docs" folder.

To use this toolbox with your code simply execute the command ADDPATH(FPATH), 
where FPATH is the absolute path of the toolbox.

[Guillermo Jiménez Arranz, 16 Jun 2021]





