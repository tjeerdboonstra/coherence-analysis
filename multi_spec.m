function results = multi_spec(x,nfft,window, fs, alpha)

% MULTI_SPEC Multivariate spectral decomposition
% Estimates several spectral measures using Fourier decomposition 
% of data matrix x, including power spectral density and coherence 
% between signals. 
%
% All spectral measures are computed using Welch's periodogram method, i.e. 
% by segmenting the signals using a moving window. Signals are dived into 
% overlapping segments with 50% overlap, each of which is detrended and windowed 
% by the WINDOW parameter then zero padded to length NFFT. The spectra X, Y, 
% and XY are estimated for each segment. Spectral coefficients are then 
% aggregated across trials (events) for the estimation of spectral measures. 
%
% ARGUMENTS:
%           x           --  signal matrix [N samples, M channels]
%           nfft        --  length of fft, zero-padded if window has less
%                           than n points
%           window      --  length of window used for spectral
%                           decomposition in samples
%           fs          --  sample frequeny
%           alpha       --  significance threshold for the confidence
%                           interval
%
%
% OUTPUTS:                  Results have the following fields
%           Px          --  spectral power
%                           [F frequencies, M channels]
%           cxx         --  complex-valued coherency [F, K] 
%                           K = M*(M-1)/2 channel combinations
%           Cxx         --  magnitude-squared coherence [F, K] 
%           freq        --  frequency vector [F,1]
%           CI          --  confidence interval for Cxx
%           combi       --  [M x M] matrix with channel combinations in Cxx
%
%
% T.W. Boonstra          16-June-2022
% University of Maastricht, The Netherlands
%
% See also FFT


% determine frequency axes
if rem(nfft,2)    % nfft odd
    select = 1:(nfft+1)/2;
else
    select = 1:nfft/2+1;   % include DC AND Nyquist
end
freq = (select - 1)'*fs/nfft;

overlap = floor(window/2); % 50% overlap between windows

% initiate variables
NF = length(freq);  % number of frequencies
NS = window;        % length of window used to segment data
NW = floor((size(x,1)-NS)/overlap)+1; % number of segments
Nx = size(x,2);     % number of channels

results.Px = zeros(NF,Nx);
results.cxx = zeros(NF,Nx*(Nx-1)/2);
results.Cxx = zeros(NF,Nx*(Nx-1)/2);

% segment data
segment = [1:NS]';
segments = zeros(NS,NW,Nx);
for t = 1:NW
    segments(:,t,:) = x(segment + (t-1)*overlap,:);
end

% confidence interval
conf_int = 1-alpha^(1/(NW-1)); %see https://doi.org/10.1016/S0165-0270(96)02214-5
conf_int = (11/9) * conf_int; % correction for using overlapping windows, see https://doi.org/10.1109/TAU.1967.1161901

% compute spectral measures
window = repmat(hanning(window),1,NW);

% spectral decomposition of x
Xx = zeros(NF,NW,Nx);
for c = 1:Nx
    X = fft(detrend(segments(:,:,c),'constant').*window,nfft);
    Xx(:,:,c) = X(select,:);
end

% compute auto spectra
XX = abs(Xx).^2;

% and average across segments to compute spectral power density
results.Px = squeeze(mean(XX,2)); % PSD of x

% compute coherence between all signal pairs
counter = 1;
for c1 = 1:Nx-1
    for c2 = c1+1:Nx
        XY = Xx(:,:,c1) .* conj(Xx(:,:,c2)); % cross-spectra
        XY = XY./sqrt(XX(:,:,c1).*XX(:,:,c2)); % complex coherency
        
        results.cxx(:,counter) = mean(XY,2); % complex-valued coherency
        results.Cxx(:,counter) = abs(mean(XY,2)).^2; % magnitude-squared coherence
        counter = counter + 1;
    end
end

results.freq = freq;
results.CI = conf_int;
results.combi = squareform([1:Nx*(Nx-1)/2]);