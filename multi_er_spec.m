function results = multi_er_spec(data,nfft,events,time,window,fs,alpha)
% MULTI_ER_SPEC Multivariate event-related spectral decomposition
% Estimates several time-frequency measures using Fourier decomposition 
% of data matrix x, including power-spectral density, inter-trial coherence
% and coherence between signals. 
%
% Time-frequency measures are computed relative to temporal events and
% aggregated across these trials. All measures are derived from the cross- 
% and autospectra estimated using Welch's periodogram method. Signals are 
% devided into overlapping segments, each of which is detrended and windowed 
% by the WINDOW parameter then zero padded to length NFFT. 
% 
% TIME defines the time points of these segments relative to the events. If 
% events are time points (e.g. stimulus onsets), EVENTS should be a vector
% and the time points in TIME are in samples. If events are intervals, EVENTS 
% should have two columns (start and end sample) and the time points in
% TIME are relative (0 corresponds to start sample and 1 to end sample).
%
% and XY are estimated for each segment. Spectral coefficients are then 
% aggregated across trials (events) for the estimation of time-frequency
% measures. 
%
% ARGUMENTS:
%           data        --  signal matrix of continuous data [N samples,
%                           M channels]
%           nfft        --  length of fft, zero-padded if window has less
%                           than n points
%           events      --  the events or markers in the continuous data 
%                           used to segment the data (in samples); the 
%                           spectral measures are computed (averaged) 
%                           across these segments
%                           [P, 1] -- P events are time points
%                           [P, 2] -- P events are intervals (start, end)
%           time        --  time points relative to event at which spectral
%                           measures are computed. If events is a vector,
%                           time is in samples. If events has two columns,
%                           time is relative (0 corresponds to start sample
%                           and 1 to end sample)
%           window      --  length of window used for spectral
%                           decomposition
%           fs          --  sample frequeny
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
%           cxx         --  event-related coherency (complex-valued)
%                           [F frequencies, T timepoints, N] 
%                           N = M*(M-1)/2 channel combinations          
%           Cxx         --  event-related coherence
%                           [F frequencies, T timepoints, N] 
%                           N = M*(M-1)/2 channel combinations
%           freq        --  frequency vector
%           time        --  time vector
%           CI          --  confidence interval for ITC and Cxx
%           num_events  --  number of events used for averaging spectral
%                           measures
%           combi       --  [M x M] matrix with channel combinations in Cxx
%
%
% T.W. Boonstra          27-March-2019
%                        06-January-2025 (updated for event intervals)
% University of Maastricht, The Netherlands
%
% See also FFT

% determine type of events
if isvector(events) % events are time points
    time_matrix = repmat(time(:),1,length(events)); % use equal time points for all events
    event_type = 'time points';
elseif size(events,2) == 2 % events are interval with start and end
    time_matrix = nan(length(time),size(events,1));
    for e = 1:size(events,1)
        time_matrix(:,e) = round(time*diff(events(e,:))); % normalize time points to the duration of each interval
    end
    events = events(:,1); % time points are relative to start sample
    event_type = 'interval';
else
    error('Events have the wrong dimensions.\n')
end

% determine frequency axes
if rem(nfft,2)    % nfft odd
    select = 1:(nfft+1)/2;
else
    select = 1:nfft/2+1;   % include DC AND Nyquist
end
freq = (select - 1)'*fs/nfft;

events = events(:)';
window = round(window); % ensure length of window is integer

if ~rem(window,2)
    window = window-1;
end
segment = [-fix(window/2):fix(window/2)]';

% check events
if min(events(:))+min(time_matrix(:))+segment(1) < 0
    error('Not enough data samples before the first event.\n')
elseif max(events(:))+max(time_matrix(:))+segment(end) > size(data,1)
    error('Not enough data samples after the last event.\n')
end

% initiate variables
NF = length(freq);
NT = size(time_matrix,1);
NE = length(events);
NS = window;
Nx = size(data,2);

results.Px = zeros(NF,NT,Nx);
results.Phx = zeros(NF,NT,Nx);
results.ITC = zeros(NF,NT,Nx);
results.cxx = zeros(NF,NT,Nx*(Nx-1)/2);
results.Cxx = zeros(NF,NT,Nx*(Nx-1)/2);

% confidence interval
conf_int = 1-alpha^(1/(NE-1)); % https://doi.org/10.1016/S0165-0270(96)02214-5

% compute event-related spectral measures seperately for each time point 
window = repmat(hanning(window),1,NE); % [NS length window, NE number of events]
for t = 1:NT 
    % extract segments relative to events
    segments = repmat(segment,1,NE) + repmat(events + time_matrix(t,:),NS,1); % [NS length window, NE number of events]
    
    % spectral decomposition of x
    Xx = zeros(NF,NE,Nx);
    for c = 1:Nx
        signal = data(:,c); % select one channel
        X = fft(detrend(signal(segments),'constant').*window,nfft); % segment, detrend, window, and perform fft
        Xx(:,:,c) = X(select,:); % put complex-valued frequency representation of segmented data of channel c
    end
    
    % compute auto spectra
    XX = abs(Xx).^2;
    
    % and average across segments to compute event-related power
    results.Px(:,t,:) = mean(XX,2); % PSD of x
    
    % compute phase
    Ph = angle(Xx);
    
    % mean phase and coherence across trials
    for c = 1:Nx
        [m,~] = circstat(Ph(:,:,c)');
        results.Phx(:,t,c) = m;  % mean phase       
        results.ITC(:,t,c) = abs(mean(Xx(:,:,c),2)./mean(sqrt(XX(:,:,c)),2)).^2; % inter-trial coherence
    end
    
    % compute cross spectra and coherence
    counter = 1;
    for c1 = 1:Nx-1
        for c2 = c1+1:Nx
            XY = Xx(:,:,c1) .* conj(Xx(:,:,c2)); % cross-spectra
            XY = XY./sqrt(XX(:,:,c1).*XX(:,:,c2)); 
            
            results.cxx(:,t,counter) = mean(XY,2); % complex-valued coherency 
            results.Cxx(:,t,counter) = abs(results.cxx(:,t,counter)).^2; % magnitude-squared coherence
            counter = counter+1;
        end
    end
end

results.freq = freq;
if strcmp(event_type,'time points')
    results.time = time/fs; % convert to seconds
elseif strcmp(event_type,'interval')
    results.time = time*100; % convert to percentage
else
    error('Something went wrong with event type.\n')
end
results.CI = conf_int;
results.num_events = NE;
results.combi = squareform(1:counter-1);
end


function [m,u]=circstat(a,b)
%CIRCSTAT(a,b) compute mean and variance of angular variables
%
%   [m,v]=CIRCSTAT(x) x is an angle in radians!!!
%   [m,v]=CIRCSTAT(x,y) x,y are read/imaginary or sin/cos components
%
%                                   (c) 8/2000 A. Daffertshofer
%
%   See also MEAN, STD.

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