%% Three-component FK analysis: 3C FK processing
% 
% This script takes the output that was computed by FK3C_Fourier and
% performs 3-component array processing on it. More information is found in
% the parameter definition cell below and in particular inside the 3C
% processor compFK3C().
% 
% 
% 
% 
% Nima Riahi, UC San Diego, nriahi@ucsd.edu
% First version out of local biotope: 2014-JUL-15
% 
% Early version, please don't distribute yet.
% 
% 


clear;

% Location where Fourier data is stored
inpdir = './SynthDat/';


%% Define 3C FK analysis parameters

% Number of time windows over which block-averaging of the SDM is performed
procpars.Nblock = 11;

% Estimate a 3C wavenumber spectrum every ntsubsample windows (this makes
% sense because the Block averaging will introduce correlations between
% adjacent time windows because they are based on strongly overlapping
% data). It is reasonable to subsample at about half the block-averaging
% duration.
procpars.ntsubsample = 6;

% Frequency range within which all bins will be processed [Hz]
procpars.frange = [.95 1.05];

% Wavenumber grid over which to compute the FK spectra [/km]
procpars.kgrid = linspace(0,1/1000,50)';

% Use polar wavenumber grid (rather than cartesian)
procpars.dopolar = true;

% Methods that are supported:
% 'DS': conventional delay-and-sum beamforming
% 'MU': MUSIC beamforming (inverse of orthogonal distance to signal subspace)
% 'HR': Capon high-resolution beamformer (estimated number of sources must be
% provided, default is 1)
procpars.method = 'DS';


% There's a lot of tuning an parameterization that can be done in the
% actual beamformer function, compFK3C().


%% Batch process all frequencies in requested frequency range

% Get first MAT file to find available frequency range
l = dir([inpdir '*.mat']);
load([inpdir l(1).name]);
f = DFT.h.f;
% Crop to requested frequency range
mask = f >= procpars.frange(1) & f<=procpars.frange(2);
f = f(mask);

for i = 1:length(f)
    
    % Load DFT data for given frequency
    f0 = f(i);
    fname = sprintf('%s%05.2f.mat',inpdir,f0);
    load(fname);
    
    % Construct time vector (start time plus signal duration converted to
    % MATLAB days)    
    t = DFT.h.t0 + (0:size(DFT.data,2)-1)*DFT.h.dt /24/3600;
    
    % Compute a series of SDMs for each time window, including
    % block-averaging (effectively a moving average)
    [f0,S] = compSDM(DFT.h.f0, reshape(DFT.data,1,size(DFT.data,2),size(DFT.data,3)) , procpars.Nblock);
    
    
    %% Subample in time
    % It is reasonable to subsample at about half the block-averaging
    % duration
    ind = 1:procpars.ntsubsample:length(t);
    t = t(ind);
    S = S(:,:,:,ind);
    
    
    %% Compute 3C wavenumber spectra
    [P,Q,kr,kth,EVals,polstates] = compFK3C(DFT.h.coords,S,...
        procpars.kgrid,...
        procpars.method,...
        procpars.dopolar);
    
    
    %% Save the results for this frequency bin somewhere and move on...
    % Use extrema2 to find local peaks in the response
    
    [];
    
end


%% Plot a part of the result

% Time window 
itime = 3;

% This plots the wavenumber response of the 3C array processing. What is
% plotted is the maximum response over all polarization states for any
% given wave vector (i.e. wavenumber and propagation azimuth). To see
% exactly what was stored in P(:,:,:) go to compFK3C(), cell 'Compute
% response'
P1 = P(:,:,itime);
subplot(211);
imagesc(kr*1000,(180/pi)*kth,P1); axis xy;
set(gca,'ytick',-180:45:180);
xlabel('wavenumber [/km]');
ylabel('propagation azimuth [deg]');
title('Array output','fontsize',14);


% This plot gives the polarization state index that led to the maximum
% response plotted above
Q1 = Q(:,:,itime);
subplot(212);
title('Polarization state ID');
imagesc(kr*1000,(180/pi)*kth,Q1); axis xy;
set(gca,'ytick',-180:45:180);
xlabel('wavenumber [/km]');
ylabel('propagation azimuth [deg]');
title('Polarization ID','fontsize',14);



% Find the strongest peak in the spectrum
[~,idx] = extrema2(P1);
idx = idx(1); % Only pick strongest maximum
[i,j] = ind2sub(size(P1),idx);

% Index of polarization state that is associated with the peak in the
% wavenumber spectrum:
Q1(i,j)

% Parameters of that polarization state ellipse. To understand the numbers
% consider an ellipse with major axis a and minor axis b.
% Entry 1: azimuth of major axis a relative to propagation direction
%           (cartesian, radians). In the current version this is always 0.
% Entry 2: Angle between major axis and its projection to the horizontal
%           plane ('dip')
% Entry 3: Ellipticity measure: 0 linear along major axis, 1 circular, 2
%           linear along minor axis
% Entry 4: Rotation of the ellipse along the major axis. Takes essentially
%           three values: 0 for retrograde Rayleigh waves, pi for prograde Rayleigh
%           waves, and pi/2 in combination with ellipticity=2 for SH waves
polstates(Q1(i,j),:)



    

