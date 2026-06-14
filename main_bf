% main_bf.m
% Rebuild of the beamforming project - Etapa A
%
% Scope of this stage (validation gate before anything else is added):
%   TX (QPSK + pulse shaping) -> SIMO channel (direct path) ->
%   independent AWGN per branch -> matched filter + integer timing ->
%   per-branch demodulation -> BER per branch vs SISO theory.
%
% Pass criterion: each branch BER must follow the theoretical QPSK
% SISO curve within Monte Carlo uncertainty. No beamforming yet.
%
% Next stages (not in this file yet):
%   B - fixed beamformer w0 = a(theta0)/M  (+6 dB check)
%   C - spatial LMS
%   D - multipath channel
%   E - sps reduction + fractional timing

clc; clear; close all;
rng(2026);

%% ------------------------------------------------------------------
%% Configuration (single place for ALL parameters)
%% ------------------------------------------------------------------

% --- Physical environment
env.c     = 1500;            % sound speed [m/s]
env.fc    = 26e3;            % carrier frequency [Hz]
env.fsADC = 156250;          % ADC sampling frequency [Hz]
env.decim = 8;               % decimation factor
env.fs    = env.fsADC / env.decim;   % baseband sampling rate [Hz]
env.lambda = env.c / env.fc; % wavelength [m]

% --- Array (fixed project geometry: 4 hydrophones, linear, 25 cm)
array.M = 4;
array.d = 0.25;                              % spacing [m]
array.positions = (0:array.M-1) * array.d;   % sensor positions [m]
array.thetaDeg = 20;                         % nominal AoA [deg]
array.thetaRad = deg2rad(array.thetaDeg);

% --- Modulation and frame
p.M = 4;                         % QPSK
p.k = log2(p.M);                 % bits per symbol
p.sps = 8;                       % samples per symbol
p.Rs = env.fs / p.sps;           % symbol rate [sym/s]
p.phaseOffset = pi/4;

p.trainingLengthSymbols = 1024;
p.preambleLengthBits = p.trainingLengthSymbols * p.k;
p.payloadLengthBits = 4000;

% --- Simulation control
p.EbNo_dB = 0:7;
p.numFrames = 20;
p.useBestSamplingPhase = true;   % integer timing: best of sps phases

% --- Pulse shaping ('rect' reproduces legacy results; 'rrc' for later)
p.pulseShape = 'rect';

switch p.pulseShape
    case 'rect'
        txFilter = ones(p.sps,1) / sqrt(p.sps);
    case 'rrc'
        rolloff = 0.25;
        span = 8;                % filter span in symbols
        txFilter = rcosdesign(rolloff, span, p.sps, 'sqrt').';
    otherwise
        error('Unknown pulse shape.');
end
rxFilter = txFilter;             % matched filter

% --- Channel: Etapa A uses ONLY the direct path
ch.Npaths = 1;
ch.delaysSec = 0;
ch.gains = 1;
ch.phases = 0;
ch.thetaDeg = array.thetaDeg;
ch.thetaRad = array.thetaRad;

%% ------------------------------------------------------------------
%% Preamble and theoretical reference
%% ------------------------------------------------------------------

preambleBits = randi([0 1], p.preambleLengthBits, 1);
preSym = qpskBitsToSymbols(preambleBits, p);

preLenSym = length(preSym);
numFrameSym = preLenSym + p.payloadLengthBits / p.k;

berTheorySISO = berawgn(p.EbNo_dB, 'psk', p.M, 'nondiff');

fprintf('--- Etapa A | direct channel | pulse = %s | sps = %d ---\n', ...
    p.pulseShape, p.sps);

%% ------------------------------------------------------------------
%% Monte Carlo sweep: BER per branch vs Eb/N0
%% ------------------------------------------------------------------

nSNR = length(p.EbNo_dB);
berBranches = zeros(nSNR, array.M);

for iSNR = 1:nSNR

    EbNo_dB = p.EbNo_dB(iSNR);

    % Eb/N0 -> SNR per SAMPLE:
    % +10log10(k) because each symbol carries k bits,
    % -10log10(sps) because symbol energy is spread over sps samples.
    snrSample_dB = EbNo_dB + 10*log10(p.k) - 10*log10(p.sps);

    nErr  = zeros(1, array.M);
    nBits = zeros(1, array.M);

    for iFrame = 1:p.numFrames

        % --- TX
        dataBits = randi([0 1], p.payloadLengthBits, 1);
        [~, ~, txWave] = txChain(dataBits, preambleBits, p, txFilter);

        % --- Channel + noise
        chanOut = applyChannelSIMO(txWave, ch, array, env);
        rxWave  = addNoiseSIMO(chanOut, snrSample_dB);

        % --- RX front end: matched filter + integer timing per branch
        rx = rxFrontEndSIMO(rxWave, preambleBits, p, rxFilter);

        % --- Per-branch demodulation and BER
        for m = 1:array.M

            frameSym = rx.frameSymbolsCell{m};

            if length(frameSym) < numFrameSym
                % Frame not fully captured: skip this branch/frame.
                continue;
            end

            frameSym = frameSym(1:numFrameSym);
            frameSym = phaseGainNormalize(frameSym, preSym);

            payloadSym = frameSym(preLenSym+1:end);
            rxBits = qpskSymbolsToBits(payloadSym, p);

            L = min(length(dataBits), length(rxBits));
            nErr(m)  = nErr(m)  + sum(dataBits(1:L) ~= rxBits(1:L));
            nBits(m) = nBits(m) + L;
        end
    end

    berBranches(iSNR,:) = nErr ./ max(nBits, 1);

    fprintf('Eb/N0 = %2d dB | branch BER = %.3e %.3e %.3e %.3e | theory = %.3e\n', ...
        EbNo_dB, berBranches(iSNR,:), berTheorySISO(iSNR));
end

%% ------------------------------------------------------------------
%% Validation plot
%% ------------------------------------------------------------------

figure('Name','Etapa A - per-branch BER vs SISO theory');
semilogy(p.EbNo_dB, berTheorySISO, 'k-', 'LineWidth', 1.8); hold on;

markers = {'o--','s--','^--','d--'};
for m = 1:array.M
    semilogy(p.EbNo_dB, berBranches(:,m), markers{m}, 'LineWidth', 1.1);
end

grid on;
xlabel('E_b/N_0 (dB)');
ylabel('BER');
title(sprintf('Etapa A | direct channel | pulse = %s', p.pulseShape));
legend('Theory SISO', 'Branch 1', 'Branch 2', 'Branch 3', 'Branch 4', ...
    'Location', 'southwest');

save('results_etapaA.mat', 'env', 'array', 'p', 'ch', ...
    'berBranches', 'berTheorySISO');
