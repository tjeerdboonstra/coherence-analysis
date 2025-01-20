# Time-frequency coherence

Code for analysing time-frequency coherence based on Fourier decomposition

tfcohf: Function to compute time-frequency coherence between a pair of continuous signals. The cross and auto spectra are smoothed to estimate time-frequency coherence. See https://doi.org/10.1186/1687-6180-2013-73 for further details

multi_spec: Function to compute coherence in multivariate data. Data are segmented with a sliding window and spectral estimates are aggregated across time.

multi_er_spec: Function to compute event-related coherence in multivariate data. Data are segmented relative to temporal events and spectral estimates are aggregated across these trials (events).

multi_er_spec_segm: Function to compute event-related coherence in multivariate data that is already segmented.
