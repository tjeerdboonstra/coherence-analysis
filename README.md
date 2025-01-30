# Time-frequency analysis

Code for performing time-frequency analysis based on Fourier decomposition. Spectral measures include power, coherence, complex-valued coherency and inter-trial coherence.

tfcohf: Function to compute time-frequency coherence between a pair of continuous signals. The cross and auto spectra are smoothed to estimate time-frequency coherence. See https://doi.org/10.1186/1687-6180-2013-73 for further details

multi_spec: Function to perform spectral analysis in multivariate data. Data are segmented with a sliding window and spectral estimates are aggregated across time.

multi_er_spec: Function to perform event-related time-frequency analysis in multivariate data. Data are segmented relative to temporal events and spectral estimates are aggregated across these trials (events). Spectral measures included event-related coherence and inter-trial coherence. 

multi_er_spec_segm: Function to perform time-frequency analysis in multivariate data that is already segmented.
