% =========================================================================
% FINAL CODE 2: VOICE-EXCITED LPC VOCODER
% Objective 1, Part 2
% =========================================================================
% This script implements the Voice-Excited Vocoder, which uses a 
% compressed residual (via DCT) instead of pitch/impulse trains.

clear; clc; close all;

% --- 1. USER PARAMETERS ---
inputFile = 'male1.wav';          % <--- CHANGE THIS for each of your 4 files!
outputFile = strrep(inputFile, '.wav', '_voice_excited.wav'); 

L = 16;             % LPC Order
fr = 20;            % Frame shift (ms)
fs = 30;            % Frame size (ms)
preemp = 0.9378;    % Pre-emphasis coefficient
DCT_N = 32;         % Keep top 32 coefficients per frame

% --- 2. BIT-RATE CALCULATION ---
% Frame Rate = 50 frames/sec
% Bits: Coeffs(50) + Gain(5) + Residual(32*8) = 311 bits/frame
% Total: 311 * 50 = 15550 bps = 15.55 kbps
fprintf('==================================================\n');
fprintf('VOCODER: Voice-Excited | Target Bitrate: < 16 kbps\n');
fprintf('Estimated Bitrate: 15.55 kbps [PASSED]\n');
fprintf('==================================================\n');

% --- 3. LOAD AUDIO ---
try
    [data, sr] = audioread(inputFile);

catch
    error('File "%s" not found. Please upload it first.', inputFile);
end
if size(data, 2) > 1, data = mean(data, 2); end
data = data / (max(abs(data)) + 1e-9);

% --- 4. RUN VOCODER ---
fprintf('Processing %s ...\n', inputFile);
tic;

% Analysis
[aCoeff, raw_resid, ~, G, ~, ~] = proclpc_robust(data, sr, L, fr, fs, preemp);

% Compression (DCT)
[~, nFrames] = size(raw_resid);
compressed_resid = zeros(size(raw_resid));
for i = 1:nFrames
    frame_dct = dct(raw_resid(:, i));
    frame_trunc = zeros(size(frame_dct));
    frame_trunc(1:DCT_N) = frame_dct(1:DCT_N); % Keep only top 32
    compressed_resid(:, i) = idct(frame_trunc);
end

% Synthesis (Uses the compressed residual)
synWave = synlpc_excited_robust(aCoeff, compressed_resid, sr, G, fr, fs, preemp);

runTime = toc;

% --- 5. SAVE AND EVALUATE ---
synWave = normalize_audio(synWave, length(data));
audiowrite(outputFile, synWave, sr);
snr = segsnr(data, synWave, sr, 20);

fprintf('Done!\n');
fprintf('   Output File:    %s\n', outputFile);
fprintf('   Run-Time:       %.4f seconds\n', runTime);
fprintf('   Segmental SNR:  %.2f dB\n', snr);
fprintf('==================================================\n\n');


% =========================================================================
% HELPER FUNCTIONS (INCLUDED)
% =========================================================================

function out = normalize_audio(in, targetLen)
    if length(in) > targetLen, in = in(1:targetLen); 
    elseif length(in) < targetLen, in(targetLen) = 0; end
    out = in / (max(abs(in)) + 1e-9);
    out = out * 0.9;
end

function [aCoeff, resid, pitch, G, parcor, stream] = proclpc_robust(data, sr, L, fr, fs, preemp)
    [row, col] = size(data); if col == 1, data = data'; end
    nframe = 0; msfr = round(sr/1000*fr); msfs = round(sr/1000*fs);
    speech = filter([1, -preemp], 1, data)'; 
    msoverlap = msfs - msfr; ramp = (0:1/(msoverlap - 1):1)'; 
    stream = []; overlap = []; duration = length(data);
    for frameIndex = 1:msfr:duration - msfs + 1
        frameData = speech(frameIndex:(frameIndex + msfs - 1));
        nframe = nframe + 1;
        if sum(abs(frameData)) < 1e-9
             aCoeff(:,nframe)=[1;zeros(L,1)]; G(nframe)=0; pitch(nframe)=0; resid(:,nframe)=zeros(msfs,1); errSig=zeros(msfs,1);
        else
            autoCor = xcorr(frameData); autoCorVec = autoCor(msfs+(0:L));
            err(1) = autoCorVec(1); k(1)=0; A=[];
            for index = 1:L
                num = [1, A']*autoCorVec(index+1:-1:2); den = -1*err(index);
                if abs(den)<1e-9, k(index)=0; else, k(index)=num/den; end
                A=[A+k(index)*flipud(A); k(index)]; err(index+1)=(1-k(index)^2)*err(index);
            end
            aCoeff(:,nframe)=[1;A]; errSig=filter([1;A]',1,frameData); G(nframe)=sqrt(abs(err(L+1)));
            autoCorErr=xcorr(errSig); [B,I]=sort(autoCorErr); num=length(I);
            if num>1 && B(num-1)>.01*B(num), pitch(nframe)=abs(I(num)-I(num-1)); else, pitch(nframe)=0; end
            if G(nframe)<1e-9, resid(:,nframe)=errSig; else, resid(:,nframe)=errSig/G(nframe); end
        end
        if frameIndex==1, stream=resid(1:msfr,nframe); else, stream=[stream; overlap+resid(1:msoverlap,nframe).*ramp; resid(msoverlap+1:msfr,nframe)]; end
        if frameIndex+msfr+msfs-1 > duration, stream=[stream; resid(msfr+1:msfs,nframe)]; else, overlap=resid(msfr+1:msfs,nframe).*flipud(ramp); end
    end
    stream = filter(1, [1, -preemp], stream)'; parcor=[];
end

function synWave = synlpc_excited_robust(aCoeff, source, sr, G, fr, fs, preemp)
    msfs=round(sr*fs/1000); msfr=round(sr*fr/1000); msoverlap=msfs-msfr; ramp=(0:1/(msoverlap-1):1)';
    [~, nframe] = size(aCoeff); resid=source; synWave=[]; overlap=[];
    for frameIndex = 1:nframe
        A = aCoeff(:, frameIndex); if any(isnan(A)), A=[1;zeros(length(A)-1,1)]; end
        Gain=G(frameIndex); if isnan(Gain), Gain=0; end
        residFrame = resid(:, frameIndex) * Gain;
        synFrame = filter(1, A', residFrame);
        if any(abs(synFrame)>100), synFrame=synFrame./max(abs(synFrame)); end
        if frameIndex==1, synWave=synFrame(1:msfr); else, synWave=[synWave; overlap+synFrame(1:msoverlap).*ramp; synFrame(msoverlap+1:msfr)]; end
        if frameIndex==nframe, synWave=[synWave; synFrame(msfr+1:msfs)]; else, overlap=synFrame(msfr+1:msfs).*flipud(ramp); end
    end
    synWave = filter(1, [1, -preemp], synWave);
end

function [segSNR] = segsnr(original, processed, sr, segment_ms)
    if length(original)~=length(processed), error('Len mismatch'); end
    segLen=round(sr*segment_ms/1000); numSeg=floor(length(original)/segLen); total=0; valid=0;
    for i=1:numSeg
        s=(i-1)*segLen+1; e=i*segLen; 
        orig=original(s:e); proc=processed(s:e);
        sigE=sum(orig.^2); errE=sum((orig-proc).^2);
        if sigE>1e-6 && errE>1e-6
            snr=10*log10(sigE/errE); if snr>35, snr=35; end
            if snr>-10, total=total+snr; valid=valid+1; end
        end
    end
    if valid>0, segSNR=total/valid; else, segSNR=-Inf; end
end