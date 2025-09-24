%% =========================================================================
%  StatFilter Imaging (Static Matched Filtering Method, 512x512, No Crop)
%
%  Description
%  -----------
%  This script implements **StatFilter**, denoting the static matched
%  filtering method. StatFilter represents the earliest paradigm in
%  array-SAR imaging, where the image is directly reconstructed via a
%  fixed analytical operator, without exploiting any prior knowledge of
%  the scene. All parameters remain fixed once defined (i.e., no adaptive
%  updates across intermediate states).
%
%  What this script does
%  ---------------------
%  1) Reads raw array-SAR echoes from 'output_rawdata.mat';
%  2) Applies simple array-plane motion correction (row-wise shift);
%  3) Performs frequency-domain static matched filtering per range bin (ID);
%  4) Aggregates results via Maximum Intensity Projection (MIP) over IDs;
%  5) Visualizes a 512x512 magnitude image in dB scale (no cropping).
%
%  Notes
%  -----
%  - Change `idRange` to focus different beat-frequency bins.
%  - Ensure the raw data struct matches the expected fields below.
%  - Requires MATLAB R2016b+ (for local functions in scripts).
% =========================================================================

clear; close all; clc;

%% ----------------------------- I/O Settings ------------------------------
% Raw data file produced by your ADC pipeline. The file must contain a struct:
%   adcRawData.data  % cell array where each cell is [numRx x numSamples]
rawMatFile = 'output_rawdata';
assert(exist([rawMatFile '.mat'],'file')==2, ...
    'Cannot find %s.mat in current folder.', rawMatFile);

load(rawMatFile, 'adcRawData');
assert(isfield(adcRawData,'data') && ~isempty(adcRawData.data), ...
    'adcRawData.data is missing or empty.');

%% ----------------------------- Data Geometry -----------------------------
numFrames   = numel(adcRawData.data);  % # frames
numTx       = 1;                       % # transmitting antennas
numRx       = 4;                       % # receiving  antennas
numSamples  = 256;                     % # time samples per chirp

% Repack raw echoes (frame-major) into [numTx*numRx*numFrames] x [numSamples]
rawEcho = zeros(numFrames*numTx*numRx, numSamples);
for ii = 1:numTx*numFrames
    rawEcho((ii-1)*numRx+1:ii*numRx, :) = squeeze(adcRawData.data{ii});
end

%% ----------------- Array-Plane Motion Correction (Static) ----------------
% Array sampling grid: Nx (horizontal), Nz (vertical)
Nx = 407;                 % horizontal sampling points (x)
Nz = 200;                 % vertical   sampling points (y)
numMimo = numTx*numRx;    % MIMO-equivalent channels (here: 4)

% Select 1T1R data (first RX of each MIMO group) into Sr
Sr = rawEcho(1:numMimo:end, :);     % [Nx*Nz] x [numSamples]
Echo = zeros(Nx*Nz, numSamples);    % corrected 1T1R echo buffer

% Static, linear motion correction model: progressive row shift
numErrPts = 43;                      % number of "movement error" points
for iz = 1:Nz
    kk = floor(iz/Nz * numErrPts);   % row-wise shift amount
    Echo((iz-1)*Nx+1:iz*Nx, :) = Sr((iz-1)*Nx+1+kk : iz*Nx+kk, :);
end

% Reshape to [Nx, Nz, numSamples]
Echo = reshape(Echo, [Nx, Nz, numSamples]);

%% ----------------------------- System Params -----------------------------
% Static (fixed) acquisition/processing parameters (no adaptive updates)
dxMm          = 1;          % sampling step along x (mm)
dyMm          = 2;          % sampling step along y (mm)
nFftSpace     = 512;        % spatial FFT size -> final image is 512x512

cLight        = physconst('lightspeed');  % 2.99792458e8 m/s
f0Hz          = (77 + 1.8) * 1e9;         % start (or ref) frequency (Hz)
fsHz          = 5e6;                      % ADC sampling rate (Hz)
tsSec         = 1/fsHz;                   % ADC sampling period (s)
slopeHzPerSec = 70.295e12;                % FMCW chirp slope (Hz/s)
tInstSec      = 6.2516e-10;               % instrument delay for range calib (s)

k0 = 2*pi*f0Hz/cLight; % reference wavenumber (rad/m)

%% -------------------------- Range FFT (Inspection) -----------------------
% Fixed time-FFT length and axis. Used by StatFilter to index range bins.
nFftTime = numSamples;                   % time-FFT length
SrFft = fft(Echo, nFftTime, 3);          % [Nx, Nz, nFftTime]

% (Optional) quick look at range spectra (static preview)
figure('Name','Range FFT (Quick Look)');
imagesc(abs(reshape(SrFft, [], nFftTime)));
xlabel('Range-FFT Bin'); ylabel('Spatial Index'); title('Range FFT Magnitude');
colormap('turbo'); colorbar;

%% ------------- StatFilter: Static Matched Filtering (per ID) -------------
% Choose a *fixed* set of range bins (IDs). No adaptive ID selection.
idRange = 12:20;                         % beat-frequency bins to focus
numIds  = numel(idRange);

% Precompute kx, ky grids (in rad/m). Convert dx, dy (mm) to meters.
% These spectral supports are fixed for the entire run (static operator).
wx = 2*pi / (dxMm*1e-3);
kx = linspace(-wx/2, wx/2, nFftSpace);   % [1 x 512]
wy = 2*pi / (dyMm*1e-3);
ky = linspace(-wy/2, wy/2, nFftSpace).'; % [512 x 1]

% Propagating wavenumber support: |k_t| <= 2k0; evanescent part set to zero.
% kProp is the static, precomputed kernel backbone used by matched filtering.
kProp = single(sqrt( max((2*k0)^2 - (kx.^2 + ky.^2), 0) ));  % [512 x 512]

% Allocate MIP stack over IDs (StatFilter results aggregated by max)
mipStack = zeros(nFftSpace, nFftSpace, numIds, 'single');

for it = 1:numIds
    thisId  = idRange(it);

    % Pick this range-bin slice and transpose to [Nz x Nx] (x as columns).
    % This is a *fixed* selection; no state-dependent reweighting.
    sarData = squeeze(SrFft(:,:,thisId)).';   % [Nz, Nx]

    % Serpentine line correction: flipping even rows to align scan direction.
    % Again, this is a static pre-process (not learned/adaptive).
    for iz = 2:2:Nz
        sarData(iz,:) = fliplr(sarData(iz,:));
    end

    % Beat-bin -> range (meters); used to build the static phase kernel.
    Rm = cLight/2 * ( thisId/(slopeHzPerSec*tsSec*nFftTime) - tInstSec );

    % Build the *static* 2-D matched filtering kernel in k-domain.
    phaseKernel = buildPhaseKernel(kProp, kx, ky, Rm, k0);

    % Size alignment (zero-pad) to match operator dimensions. Static padding.
    [sarDataPad, phaseKernelPad] = padToMatch(sarData, phaseKernel);

    % Frequency-domain multiplication (static matched filtering operator).
    sarDataF = fft2(sarDataPad, nFftSpace, nFftSpace);
    imgCplx  = ifft2( sarDataF .* phaseKernelPad );  % complex image 512x512

    % Magnitude of StatFilter output (per ID). No adaptive rescaling.
    mipStack(:,:,it) = abs(imgCplx);
end

% Maximum-Intensity Projection (StatFilter results aggregated over IDs).
mipImg = max(mipStack, [], 3);   % 512x512

%% ------------------------------ Visualization ----------------------------
% Fixed (static) dB visualization settings.
dbCfg.clipRange    = 40;
dbCfg.displayRange = 40;

dbImg = toDbImage(mipImg, dbCfg.clipRange, dbCfg.displayRange);

figure('Name','StatFilter Image (MIP, 512x512)');
imagesc(dbImg.data, [-dbImg.range, 0]);
axis image off;
colormap('turbo'); colorbar;
title(sprintf('StatFilter (Static Matched Filtering) — MIP of ID %d–%d, 512×512', ...
      idRange(1), idRange(end)));

%% ============================== Local Functions ==========================
function phaseKernel = buildPhaseKernel(kProp, kx, ky, Rm, k0)
%BUILDPHASEKERNEL Build the static 2-D matched filtering kernel in k-domain.
%   StatFilter uses a fixed analytical operator. This function constructs
%   a centered (fftshift-ed) kernel with evanescent components removed.
%
%   Inputs:
%     kProp : [Ny x Nx]  sqrt((2k0)^2 - kx^2 - ky^2), nonnegative (static)
%     kx    : [1 x Nx]   spectral grid along x (rad/m)
%     ky    : [Ny x 1]   spectral grid along y (rad/m)
%     Rm    : scalar     range (m) corresponding to the chosen range bin
%     k0    : scalar     reference wavenumber (rad/m)
%
%   Output:
%     phaseKernel : [Ny x Nx], centered and masked for |k_t|>2k0
%
%   Note:
%     The multiplication by kProp is a standard dispersion-compensation
%     choice in ω–k domain processing; it remains fixed across the run.

    % Evanescent region mask (|k_t| > 2k0) -> zero out (static support)
    [KX,KY] = meshgrid(kx, ky);
    outCone = (KX.^2 + KY.^2) > (2*k0)^2;

    % Base phase term with static propagation kernel
    phase0 = exp(-1i * Rm .* kProp);
    phase0(outCone) = 0;

    % Fixed amplitude weighting (no adaptive modulation)
    phaseKernel = kProp .* phase0;

    % Center in frequency so that ifft2 returns a spatially centered image
    phaseKernel = fftshift(fftshift(phaseKernel,1),2);
end

function [Aout, Bout] = padToMatch(A, B)
%PADTOMATCH Zero-pad two matrices to the same size (centered).
%   This is a static alignment step; no data-dependent resizing/warping.
%
%   Inputs:
%     A, B : 2-D arrays (real/complex). The smaller one is padded to match.
%   Outputs:
%     Aout, Bout : size-matched arrays.

    [ay, ax] = size(A);
    [by, bx] = size(B);
    Aout = A; Bout = B;

    % Match X dimension
    if bx > ax
        pre  = floor((bx - ax)/2);
        post = ceil( (bx - ax)/2);
        Aout = padarray(Aout, [0 pre], 0, 'pre');
        Aout = padarray(Aout, [0 post], 0, 'post');
    elseif ax > bx
        pre  = floor((ax - bx)/2);
        post = ceil( (ax - bx)/2);
        Bout = padarray(Bout, [0 pre], 0, 'pre');
        Bout = padarray(Bout, [0 post], 0, 'post');
    end

    % Match Y dimension
    [ay, ax] = size(Aout);
    [by, bx] = size(Bout);
    if by > ay
        pre  = floor((by - ay)/2);
        post = ceil( (by - ay)/2);
        Aout = padarray(Aout, [pre 0], 0, 'pre');
        Aout = padarray(Aout, [post 0], 0, 'post');
    elseif ay > by
        pre  = floor((ay - by)/2);
        post = ceil( (ay - by)/2);
        Bout = padarray(Bout, [pre 0], 0, 'pre');
        Bout = padarray(Bout, [post 0], 0, 'post');
    end
end

function out = toDbImage(X, clipRange, displayRange)
%TODBIMAGE Convert a magnitude image to normalized dB scale (static view).
%   out = toDbImage(X, clipRange)
%   out = toDbImage(X, clipRange, displayRange)
%
%   - X: real/complex 2-D image (magnitude will be used).
%   - clipRange (dB): values below -clipRange are set to -displayRange.
%   - displayRange (dB): colorbar span [-displayRange, 0].
%
%   Returns a struct:
%       out.data  : dB image normalized to its global max (= 0 dB peak)
%       out.range : displayRange (for convenience)

    if nargin < 3
        displayRange = clipRange;
    end

    Xmag = abs(X);
    xmax = max(Xmag(:));
    if xmax == 0, xmax = eps; end

    XdB = 20*log10( max(Xmag, eps) / xmax );
    XdB(XdB < -clipRange) = -displayRange;

    out = struct('data', XdB, 'range', displayRange);
end
