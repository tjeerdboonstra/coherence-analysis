function results = multi_er_spec_segm(x,nfft,event,time,window,fs,alpha)
% MULTI_ER_SPEC_SEGM Multivariate event-related spectral decomposition
% Estimates several time-frequency measures using Fourier decomposition
% of segmented data (into trials), including power-spectral density, 
% inter-trial coherence and coherence between signals. Time-frequency 
% measures are computed across these trials.
%
% All measures are derived from the cross- and autospectra estimated using
% Welch's periodogram method. Signals are divided into overlapping segments,
% each of which is detrended and windowed by the WINDOW parameter then
% zero padded to length NFFT. TIME defines the time points of these segments
% relative to the EVENT. The spectra X, Y, and XY are estimated for each
% segment. Spectral coefficients are then aggregated across trials for the
% estimation of time-frequency measures.
%
% ARGUMENTS:
%           x           --  segmented data [N samples, M channels, P trials]
%           nfft        --  length of fft, zero-padded if window has less
%                           than n points
%           event       --  sample defined as t0
%           time        --  time points relative to event
%           window      --  length of window used for spectral
%                           decomposition
%           fs          --  sample frequency
%           alpha       --  significance threshold for the confidence
%                           interval
%
%
% OUTPUTS:                  Results have the following fields
%           Px          --  event-related power
%                           [F frequencies, T timepoints, M channels]
%           Phx         --  mean phase
%                           [F frequencies, T timepoints, M channels]
%           ITC         --  inter-trial coherence
%                           [F frequencies, T timepoints, M channels]
%           Cxx         --  event-related coherence
%                           [F frequencies, T timepoints, N]
%           cxx         --  event-related complex-valued coherency
%                           [F frequencies, T timepoints, N]
%                           N = M*(M-1)/2 channel combinations
%           freq        --  frequency vector
%           time        --  time vector
%           CI          --  confidence interval for ITC and Cxx
%           combi       --  [M x M] matrix with channel combinations in Cxx
%
%
% T.W. Boonstra          10-July-2019
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

% define segment
if ~rem(window,2)
    window = window-1;
end
segment = [-fix(window/2):fix(window/2)]';

% convert time to samples
time = round(time*fs);

% check if the time points fit within the segmented trials
if event-time(1)-segment(1)<0 || event+time(end)+segment(end)>size(x,1)
    event-time(1)-segment(1)
    event+time(end)+segment(end)
    size(x,1)
    error('Requested time points fall outside segmented data.')
end

% initiate variables
NF = length(freq);
NT = length(time);
NE = size(x,3);
NCH = size(x,2);

results.Ax = zeros(NF,NT,NCH);
results.Px = zeros(NF,NT,NCH);
results.Phx = zeros(NF,NT,NCH);
results.ITC = zeros(NF,NT,NCH);
results.Cxx = zeros(NF,NT,NCH*(NCH-1)/2);
results.cxx = zeros(NF,NT,NCH*(NCH-1)/2);

% confidence interval
conf_int = 1-alpha^(1/(NE-1));

% compute event-related spectral measures
window = repmat(hanning(window),1,NCH);
for t = 1:NT
    % extract segment relative to event
    ii = event + time(t) + segment;

    % spectral decomposition of x
    Xx = zeros(NF,NE,NCH);
    for n = 1:NE
        X = fft(detrend(x(ii,:,n),'constant').*window,nfft);
        Xx(:,n,:) = X(select,:);
    end

    % compute auto spectra
    XX = abs(Xx).^2;

    % and average across segments to compute event-related power
    results.Ax(:,t,:) = mean(abs(Xx),2); % average amplitude spectrum
    results.Px(:,t,:) = mean(XX,2); % PSD

    % compute phase
    Ph = angle(Xx);

    % mean phase and coherence across trials
    for c = 1:NCH
        [m,~] = circstat(Ph(:,:,c)');
        results.Phx(:,t,c) = m;  % mean phase
        results.ITC(:,t,c) = abs(mean(Xx(:,:,c),2)./mean(sqrt(XX(:,:,c)),2)).^2; % inter-trial coherence
    end

    % compute cross spectra and coherence
    counter = 1;
    for c1 = 1:NCH-1
        for c2 = c1+1:NCH
            XY = Xx(:,:,c1) .* conj(Xx(:,:,c2)); % cross-spectra
            XY = XY./sqrt(XX(:,:,c1).*XX(:,:,c2)); % complex coherency

            results.Cxx(:,t,counter) = abs(mean(XY,2)).^2; % magnitude-squared coherence
            results.cxx(:,t,counter) = mean(XY,2); % complex-valued coherency
            counter = counter+1;
        end
    end

end

results.freq = freq;
results.time = time/fs;
results.CI = conf_int;
results.combi = squareform(1:counter-1);
end


function [m,u]=circstat(a,b)
%CIRCSTAT(a,b) compute mean and variance of angular variables

if nargin==1
    sa=sin(a);
    ca=cos(a);
else
    sa=a;
    ca=b;
    for k=1:size(sa,2)
        rr=sqrt(sa(:,k).^2+ca(:,k).^2);
        sa(:,k)=sa(:,k)./rr;
        ca(:,k)=ca(:,k)./rr;
    end
end
sa=mean(sa);
ca=mean(ca);

m=mod(atan2(sa,ca),2*pi);

if nargout==2
    u=m;
    for k=1:length(u)
        u(k)=1-sqrt(sa(:,k)^2+ca(:,k)^2);
    end
end
end