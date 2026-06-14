function rxWave = addNoiseSIMO(chanOut, snrSample_dB)
% addNoiseSIMO
% Adds independent complex AWGN to each hydrophone branch.
%
% This is the ONLY noise model of the rebuilt project.
% Physical justification: thermal/ambient noise at different sensors
% is statistically independent. Spatially correlated noise (the old
% "common noise" test) is NOT a valid nominal model and was removed.
%
% Inputs:
%   chanOut      : N x M matrix, channel output (one column per hydrophone)
%   snrSample_dB : SNR per sample in dB (Eb/N0 already converted by caller)
%
% Output:
%   rxWave       : N x M matrix, noisy received waveform

    [~, M] = size(chanOut);
    rxWave = zeros(size(chanOut));

    for m = 1:M
        % 'measured' makes awgn estimate the signal power of THIS branch,
        % so every branch gets the same SNR regardless of channel gain.
        rxWave(:,m) = awgn(chanOut(:,m), snrSample_dB, 'measured');
    end

end
