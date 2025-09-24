%% =========================================================================
%  2-D RMA Imaging with Multi-Range MIP (512x512, No Crop)

%  Description
%  -----------
%  This script performs 2-D Range Migration Algorithm (RMA) imaging for a
%  planar scan (Nx x Nz) using single-TX/multi-RX MIMO-equivalent data.
%  It computes range FFT, applies per-range-bin (ID) RMA focusing, stacks
%  the magnitude images over a set of range bins, and visualizes the
%  Maximum Intensity Projection (MIP) along the range-bin dimension.
%
%  Key features
%  -----------
%  - Reads raw data layout used in 'output_rawdata' (Matlab .mat).
%  - Performs simple motion correction in array plane (row shift).
%  - RMA focusing with frequency-domain phase compensation.
%  - No cropping: final image is exactly 512x512.
%  - dB visualization with configurable clipping/display ranges.
%
%  Notes
%  -----
%  - Change `idRange` to sweep different beat-frequency bins.
%  - Ensure your raw data struct matches fields used below.
%  - Requires MATLAB R2016b+ (for local functions in scripts).
% =========================================================================

clear; close all; clc;

%% ----------------------------- I/O Settings ------------------------------
rawMatFile = 'output_rawdata';   % .mat file produced by your ADC pipeline
assert(exist([rawMatFile '.mat'],'file')==2, ...
    'Cannot find %s.mat in current folder.', rawMatFile);

load(rawMatFile, 'adcRawData');  % expected: adcRawData.data is a cell array
assert(isfield(adcRawData,'data') && ~isempty(adcRawData.data), ...
    'adcRawData.data is missing or empty.');

%% ----------------------------- Data Geometry -----------------------------
numFrames   = numel(adcRawData.data);  % number of frames
numTx       = 1;                       % # of transmitting antennas
numRx       = 4;                       % # of receiving  antennas
numSamples  = 256;                     % # of time samples per chirp

% Repack raw echoes (frame-major) into [numTx*numRx*m] x [numSamples]
rawEcho = zeros(numFrames*numTx*numRx, numSamples);
for ii = 1:numTx*numFrames
    rawEcho((ii-1)*numRx+1:ii*numRx, :) = squeeze(adcRawData.data{ii});
end

%% ------------------- Motion Correction in Array Plane --------------------
% Array sampling: Nx (horizontal), Nz (vertical)
Nx = 407;                 % horizontal sampling points (x)
Nz = 200;                 % vertical   sampling points (y)
numMimo = numTx*numRx;    % MIMO-equivalent channels (here: 4)

% Select 1T1R data (first RX of each MIMO group) into Sr
Sr = rawEcho(1:numMimo:end, :);     % [Nx*Nz] x [numSamples]
Echo = zeros(Nx*Nz, numSamples);    % corrected 1T1R echo buffer

% Simple linear motion correction: progressive row shift
numErrPts = 43;                      % number of "movement error" points
for iz = 1:Nz
    kk = floor(iz/Nz * numErrPts);   % row-wise shift amount
    Echo((iz-1)*Nx+1:iz*Nx, :) = Sr((iz-1)*Nx+1+kk : iz*Nx+kk, :);
end

% Reshape to [Nx, Nz, numSamples]
Echo = reshape(Echo, [Nx, Nz, numSamples]);

%% ----------------------------- System Params -----------------------------
dxMm          = 1;          % sampling step along x (mm)
dyMm          = 2;          % sampling step along y (mm)
nFftSpace     = 512;        % spatial FFT size -> final image is 512x512

cLight        = physconst('lightspeed');  % 2.9979e8 m/s
f0Hz          = (77 + 1.8) * 1e9;         % start (or ref) frequency
fsHz          = 5e6;                      % ADC sampling rate
tsSec         = 1/fsHz;                   % ADC sampling period
slopeHzPerSec = 70.295e12;                % FMCW chirp slope (Hz/s)
tInstSec      = 6.2516e-10;               % instrument delay for range calib

k0 = 2*pi*f0Hz/cLight; % reference wavenumber

%% -------------------------- Range FFT (Inspection) -----------------------
nFftTime = numSamples;                   % time-FFT length
SrFft = fft(Echo, nFftTime, 3);          % [Nx, Nz, nFftTime]

% (Optional) quick look at range spectra
figure('Name','Range FFT (Quick Look)'); 
imagesc(abs(reshape(SrFft, [], nFftTime)));
xlabel('Range-FFT Bin'); ylabel('Spatial Index'); title('Range FFT Magnitude');
colormap('turbo'); colorbar;

%% -------------------- RMA Imaging over Multiple Range IDs ----------------
idRange = 12:20;                         % beat-frequency bins to focus
numIds  = numel(idRange);

% Precompute kx, ky grids (in rad/m). Convert dx, dy (mm) to meters.
wx = 2*pi / (dxMm*1e-3);
kx = linspace(-wx/2, wx/2, nFftSpace);   % [1 x 512]
wy = 2*pi / (dyMm*1e-3);
ky = linspace(-wy/2, wy/2, nFftSpace).'; % [512 x 1]

% Propagation wavenumber support (|k_t| <= 2k0); evanescent set to zero
% kProp is 2D matrix over ky,kx
kProp = single(sqrt( max((2*k0)^2 - (kx.^2 + ky.^2), 0) ));  % [512 x 512]

% Allocate MIP stack
mipStack = zeros(nFftSpace, nFftSpace, numIds, 'single');

for it = 1:numIds
    thisId   = idRange(it);
    % pick this range-bin slice, transpose to [Nz x Nx] so x is columns
    sarData  = squeeze(SrFft(:,:,thisId)).';   % [Nz, Nx]

    % serpentine (snake) correction: flip even rows to align scan direction
    for iz = 2:2:Nz
        sarData(iz,:) = fliplr(sarData(iz,:));
    end

    % Beat-bin -> range (meters), then phase kernel for RMA
    Rm = cLight/2 * ( thisId/(slopeHzPerSec*tsSec*nFftTime) - tInstSec );

    % Build the 2-D phase compensation kernel (centered in k-domain)
    phaseKernel = buildPhaseKernel(kProp, kx, ky, Rm, k0);

    % Size alignment (zero-pad either data or kernel to match)
    [sarDataPad, phaseKernelPad] = padToMatch(sarData, phaseKernel);

    % Frequency-domain multiplication (RMA focusing)
    sarDataF = fft2(sarDataPad, nFftSpace, nFftSpace);
    imgCplx  = ifft2( sarDataF .* phaseKernelPad );  % complex image 512x512

    mipStack(:,:,it) = abs(imgCplx);                 % magnitude image
end

% Maximum-Intensity Projection (along ID dimension)
mipImg = max(mipStack, [], 3);   % 512x512

%% ------------------------------ Visualization ----------------------------
% dB image with clipping/display range = 40 dB
dbCfg.clipRange    = 40;
dbCfg.displayRange = 40;

dbImg = toDbImage(mipImg, dbCfg.clipRange, dbCfg.displayRange);

figure('Name','SAR Image - 2D RMA (MIP, 512x512)');
imagesc(dbImg.data, [-dbImg.range, 0]);
axis image off;
colormap('turbo'); colorbar;
title(sprintf('SAR Image - 2D RMA (MIP of ID %d–%d, 512×512)', idRange(1), idRange(end)));

%% ============================== Local Functions ==========================
function phaseKernel = buildPhaseKernel(kProp, kx, ky, Rm, k0)
%BUILDPHASEKERNEL Build the centered 2-D RMA phase compensation kernel.
%   kProp : [Ny x Nx] propagating wavenumber, sqrt((2k0)^2 - kx^2 - ky^2)
%   kx    : [1 x Nx] spectral grid along x (rad/m)
%   ky    : [Ny x 1] spectral grid along y (rad/m)
%   Rm    : range (m) corresponding to a beat-bin
%   k0    : reference wavenumber at f0 (rad/m)
%
% Returns:
%   phaseKernel : [Ny x Nx], centered (fftshift) and masked for |k_t|>2k0

    % Evanescent region mask (|k_t| > 2k0) -> zero out
    [KX,KY] = meshgrid(kx, ky);
    outCone = (KX.^2 + KY.^2) > (2*k0)^2;

    % Base phase term and amplitude weighting (standard ω-domain RMA)
    phase0 = exp(-1i * Rm .* kProp);
    phase0(outCone) = 0;

    % Common choice: multiply by kProp to compensate dispersion (optional)
    phaseKernel = kProp .* phase0;

    % Center in frequency (so that ifft2 returns spatially centered image)
    phaseKernel = fftshift(fftshift(phaseKernel,1),2);
end

function [Aout, Bout] = padToMatch(A, B)
%PADTOMATCH Zero-pad 2 matrices to the same size (centering the content).
%   A, B : input 2-D arrays (real/complex). The smaller one will be padded
%          around to match the larger one's size (both dims).
% Returns:
%   Aout, Bout : size-matched arrays.

    [ay, ax] = size(A);
    [by, bx] = size(B);
    Aout = A; Bout = B;

    % Match X
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

    % Match Y
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
%TODBIMAGE Convert a real/complex image to normalized dB scale.
%   out = toDbImage(X, clipRange)
%   out = toDbImage(X, clipRange, displayRange)
%
%   - X: real or complex 2-D image.
%   - clipRange (dB): values below -clipRange will be set to -displayRange.
%   - displayRange (dB): colorbar span [-displayRange, 0].
%
%   Returns a struct:
%       out.data  : dB image of X normalized to its global max (= 0 dB peak)
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
