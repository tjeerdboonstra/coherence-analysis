# time-freq

Code for analysing time-frequency coherence based on Fourier decomposition

tfcohf: Function to compute time-frequency coherence between two continuous signal. The cross and auto spectra are smoothed to estimate time-frequency coherence. See https://doi.org/10.1186/1687-6180-2013-73 for further details

multi_er_spec: Function to compute event-related coherence in multivariate data. Data are segment relative to temporal events and spectral estimates are aggregated across these trials.
