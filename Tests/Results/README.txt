RESULTS
==================

- Zero-Padding: zero-padding was applied to the signal before filtering.
  If the word is not included in the brackets of the file name, zero-padding
  was then not applied.
- blow: the TOB spectrum is calculated over a piling blow.
- wnoise: the TOB spectrum is calculated over white Gaussian noise.
- Two-Sided: the filtering (using myfiltfilt) is applied to the signal
  zero-padded as both sides. This is the method used in the ANSI Toolbox.
  The reason for padding with zeroes at both ends of the signal is that
  myfiltfilt (which is a modified version of filtfilt) filters the signal
  first in the forward and then in the backward direction, resulting in
  energy leakage over time towards the two ends of the signal.
- One-Sided: the filtering (using myfiltfilt) is applied to the signal
  zero-padded to the right only. This example is only for testing. The
  comparison between the Two-Sided and One-Sided figures shows that 
  zero-padding at both ends of the signal when myfiltfilt (or filtfilt) 
  is used is necessary to achieve accurate results.