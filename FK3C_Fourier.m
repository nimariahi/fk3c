%% Three-component FK analysis: Fourier computation
% 
% This script loads a trace gather of 3*n traces and computes Fourier
% coefficients for all traces. The coefficients of all sensors for all time 
% are stored per frequency bin.
% 
% The storage format of the Fourier coefficients is understood by the next
% script that computes the 3C FK spectra
% 
% 
% Nima Riahi, UC San Diego, nriahi@ucsd.edu
% First version out of local biotope: 2014-JUL-15
% 
% Early version, please don't distribute yet.
% 

clear;


%% Data locations

% Fourier coefficient data is stored here
outdir = './SynthDat/';

% Create directory
if ~exist(outdir,'dir')
    mkdir(outdir);
end


%% Load data

% DAT = FK3C_dataload();
load('./SynthDat.mat');

% 'DAT' is a structure that contains a trace gather from all 3*n receivers
% of the n element 3C array. It has the following contents:
% 
% DAT.h             : Header
% DAT.h.coords      : (n, 2) orthogonal projection (if possible) of coordinates of receivers
% DAT.h.stations    : (n, p) custom table containing metadata about the
%                       receivers. This field is ignored in processing but
%                       may be useful for analysis and troubleshooting
% DAT.h.t0          : (1, 1) Start time of trace gather in MATLAB days
% DAT.h.dt          : (1, 1) Sampling period of trace gather
% 
% 
% DAT.procpars      : (struct) Processing parameters that went into the data, e.g.
%                       bandpass filter, etc. (same applies as with the
%                       field 'h.stations')
% 
% 
% DAT.data          : (L, 3*n) Actual trace gather. Each column contains
%                       the time recording of a receivers. There are 3*n
%                       receivers for an n element 3C array. THE ORDERING
%                       IS KEY!
%                       - First n MUST be the East components,
%                       - Next n MUST be the North components with the
%                       same order as the East
%                       - Last n MUST be the vertical components with same
%                       ordering as East
%                       
%                       It is CRUCIAL that the field DAT.h.coords ALSO has
%                       the identical ordering as the columns in DAT.data
% 

%  While we're at it, here's what the structure of the variable that this
%  script will write to disk
% 
% DFT.h
% DFT.h.coords      : Same as DAT
% DFT.h.stations    : Same as DAT
% DFT.h.t0          : As above but adjusted for the larger time step
% DFT.h.dt          : As above, but now the step between windows
% DFT.h.f           : (n, 1) frequency vector from FFT
% DFT.h.f0          : (1, 1) specific frequency bin the data of which is
%                       stored in this DFT variable.
% 
% DFT.procpars      : Processing paramters set in this script
% 
% DFT.data          : (1, L, 3*n) 3D matrix containing spectral coefficients
%                       for all times and all receivers. The first
%                       dimension is a rudiment from ealier versions that
%                       also contained all frequencies. This was abandoned
%                       for faster processing.
% 
% 



%% Fourier parameters

% Window and step size are defined in samples and have to be converted from
% an appropriate physical time
procpars.nwin = 512; % Sample size of windows (use power of two to produce less confusion about normalization)
procpars.nstep = 512/2; % Sample size of time step

% Tukey window is used as preprocessing to compute spectral coefficients
procpars.tukeyfrac = 0.2; 

% Frequency range that will be stored
procpars.frange = [0 3];


%% Loop through all 3C data channels and compute Fourier coefficients
sps = 1/DAT.h.dt;
Nloc = size(DAT.h.coords,1);
for i = 1:Nloc
    

    % Compute short-time FFT with tukey window (incl. power correction).
    % Set imaginary parts to 'real' NaN (there are imaginary NaN's which
    % imagesc cannot handle).
    win = tukeywin(procpars.nwin,procpars.tukeyfrac); % for the 4096 cache
    [DFTE,f_,~,~] = spectrogram(DAT.data(:,i)       ,win,procpars.nwin-procpars.nstep,procpars.nwin,sps);
    [DFTN,f_,~,~] = spectrogram(DAT.data(:,i+Nloc)  ,win,procpars.nwin-procpars.nstep,procpars.nwin,sps);
    [DFTZ,f_,~,~] = spectrogram(DAT.data(:,i+2*Nloc),win,procpars.nwin-procpars.nstep,procpars.nwin,sps);
    
    % Initialize spectrogram containers for each component
    if ~exist('DFTZS','var')
        
        % Create time and frequency vectors
%         t = (0:size(DFTZ,2)-1)*DAT.h.dt / 24/3600 + DAT.h.t0;
        
        fmask = f_>=procpars.frange(1) & f_<=procpars.frange(2);
        f = f_(fmask);
        
        % Normalization factor, relating FFT output to PSD, see MATLAB
        % documentation for 'spectrogram'
        procpars.pnorm = 1/ ( (1/DAT.h.dt)*sum( win.^2 ) );
                
        % Actual initialization
        DFTES = zeros(sum(fmask),size(DFTZ,2),size(DAT.h.coords,1));
        DFTNS = zeros(sum(fmask),size(DFTZ,2),size(DAT.h.coords,1));
        DFTZS = zeros(sum(fmask),size(DFTZ,2),size(DAT.h.coords,1));
    end
    
    % Store DFT data into spectrogram containers
    DFTES(:,:,i) = DFTE(fmask,:);
    DFTNS(:,:,i) = DFTN(fmask,:);
    DFTZS(:,:,i) = DFTZ(fmask,:);
    
end

% Some clean up

% Reshape it such that a 3-dim matrix results with
% 1- freq / 2- time / 3- data sources K*E,K*N,K*Z
DFTS = reshape([DFTES(:);DFTNS(:);DFTZS(:)], size(DFTES,1), size(DFTES,2) , size(DFTES,3)*3 );
clear DFTES DFTNS DFTZS;
clear DFTE DFTN DFTZ;


%% Save DFT data to disk 
% Note that the 3C array processor operates on 1 or a few frequency bins
% only

% Create directory (incl. subdirectories) if necessary
if ~exist(outdir,'dir')
    mkdir(outdir);
end

% Store spectra of all 3C channels in the array per frequency bin
% All data and metadata are stored in a structure that is supposed to
% contain all relevant information of where the data came from and what
% processing it received to this point.
for i = 1:size(DFTS,1)
    
    % Extract spectral data for a given frequency
    DFT.data = DFTS(i,:,:);
    
    % Store Fourier processing parameters in a header
    DFT.procpars = procpars;
    
    % Add frequency of this spectral 'layer' to header
    DFT.h.f0 = f(i);
    
    % Add header information (start time, dt, coordinates, etc.)
    DFT.h.t0 = DAT.h.t0 + procpars.nwin*DAT.h.dt/2/24/3600; % Shift start time by half the window length (/24/3600 to convert sec to MATLAB days)
    DFT.h.coords = DAT.h.coords;
    DFT.h.stations = DAT.h.stations;
    DFT.h.dt = procpars.nstep * DAT.h.dt;
    DFT.h.f = f; % Also add the full frequency vector
    
    % The file name contains the frequency. Depending on frequency a naming
    % convention of mHz or dHz may procide more tractable file names
    fname = sprintf('%s%05.2f.mat',outdir,DFT.h.f0);
    save(fname,'DFT');
    
    % Clean up
    clear DFT;
    
end

