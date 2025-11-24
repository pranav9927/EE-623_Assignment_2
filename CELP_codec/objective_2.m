% ==============================================================================
% EE623 Assignment 2 - Objective 2: CELP Codec Implementation
% Target Bitrate: ~14 kbps
% Platform: MATLAB Online Compatible
% ==============================================================================
clc; clear; close all;

% --- CONFIGURATION ---
files = {'male1.wav', 'male2.wav', 'female1.wav', 'female2.wav'};
target_bitrate_kbps = 14; 

% Storage for average metrics
total_snr = 0;
total_lsd = 0;
total_runtime = 0;
valid_files_count = 0;

fprintf('==========================================================\n');
fprintf('EE623 CELP Codec Implementation\n');
fprintf('Target Bitrate: %d kbps\n', target_bitrate_kbps);
fprintf('Metric: Log-Spectral Distance (LSD) used in place of PESQ\n');
fprintf('==========================================================\n');

% --- MAIN PROCESSING LOOP ---
for i = 1:length(files)
    fname = files{i};
    
    % Check if file exists
    if ~isfile(fname)
        fprintf('\n[WARNING] File %s not found. Skipping.\n', fname);
        continue;
    end
    
    valid_files_count = valid_files_count + 1;
    fprintf('\nProcessing File %d/%d: %s ...\n', i, length(files), fname);
    
    % 1. READ AUDIO
    [s, fs] = audioread(fname);
    
    % 2. PRE-PROCESSING (Resample to 8kHz narrowband standard)
    if fs ~= 8000
        s = resample(s, 8000, fs);
        fs = 8000;
    end
    s = s(:); % Ensure column vector
    
    % 3. RUN CELP CODEC
    tic;
    [s_synth, stats] = run_celp_14kbps(s, fs); 
    run_time = toc;
    
    % 4. CALCULATE METRICS
    % A. Segmental SNR
    val_snr = calc_segsnr(s, s_synth, fs);
    
    % B. Log-Spectral Distance (LSD)
    val_lsd = calc_lsd(s, s_synth, fs);
    
    % C. Actual Bitrate Calculation
    duration_sec = length(s)/fs;
    actual_bps = stats.total_bits / duration_sec;
    actual_kbps = actual_bps / 1000;
    
    % 5. SAVE OUTPUT
    out_name = strrep(fname, '.wav', '_celp_14kbps.wav');
    audiowrite(out_name, s_synth, fs);
    
    % 6. DISPLAY RESULTS
    fprintf('  > Status:        Completed\n');
    fprintf('  > Runtime:       %.4f s\n', run_time);
    fprintf('  > Bitrate:       %.2f kbps (Target: 14.00)\n', actual_kbps);
    fprintf('  > SegSNR:        %.2f dB\n', val_snr);
    fprintf('  > LSD Score:     %.2f dB (Lower is better)\n', val_lsd);
    fprintf('  > Saved as:      %s\n', out_name);
    
    % Accumulate
    total_snr = total_snr + val_snr;
    total_lsd = total_lsd + val_lsd;
    total_runtime = total_runtime + run_time;
    
    % Optional: Plot Comparison for the first file only
    if i == 1 || i==2 || i==3 || i==4
        figure('Name', ['Comparison: ' fname]);
        subplot(2,1,1); plot(s); title(['Original: ' fname]); grid on; axis tight;
        subplot(2,1,2); plot(s_synth); title('Synthesized (CELP)'); grid on; axis tight;
        drawnow;
    end
end

% --- FINAL REPORT SUMMARY ---
if valid_files_count > 0
    fprintf('\n==========================================================\n');
    fprintf('FINAL AVERAGE RESULTS (Over %d files)\n', valid_files_count);
    fprintf('==========================================================\n');
    fprintf('Avg Runtime:       %.4f s\n', total_runtime / valid_files_count);
    fprintf('Avg Segmental SNR: %.2f dB\n', total_snr / valid_files_count);
    fprintf('Avg LSD Score:     %.2f dB\n', total_lsd / valid_files_count);
    fprintf('----------------------------------------------------------\n');
else
    fprintf('\n[ERROR] No valid audio files found. Please add male1.wav, etc.\n');
end


% ==============================================================================
% LOCAL FUNCTIONS
% ==============================================================================

function [s_synth, stats] = run_celp_14kbps(s, fs)
    % CELP Logic Tuned for High Bitrate (14 kbps)
    
    % --- PARAMETERS ---
    frame_ms = 10; 
    L_frame = round((frame_ms/1000) * fs); % 80 samples
    LPC_ORDER = 10;
    
    N = length(s);
    s_synth = zeros(N, 1);
    exc_mem = zeros(N, 1); % Adaptive codebook memory
    total_bits = 0;
    
    % Pre-emphasis
    preemp = 0.9378;
    s_pre = filter([1 -preemp], 1, s);
    
    num_frames = floor(N / L_frame);
    
    for n = 1:num_frames
        idx_start = (n-1)*L_frame + 1;
        idx_end = idx_start + L_frame - 1;
        
        frame_curr = s_pre(idx_start:idx_end);
        
        % 1. LPC Analysis [FIX ADDED: Check for silence to avoid NaNs]
        if sum(abs(frame_curr)) < 1e-9
            A = [1 zeros(1, LPC_ORDER)]; % Default filter for silence
        else
            [A, ~] = lpc(frame_curr, LPC_ORDER);
        end
        
        % Check for invalid A (NaNs) just in case
        if any(isnan(A))
            A = [1 zeros(1, LPC_ORDER)];
        end
        
        % 2. Calculate Residual
        residual = filter(A, 1, frame_curr);
        
        % 3. Adaptive Codebook (Pitch) Search
        start_search = max(1, idx_start - 145);
        best_lag = 20; 
        
        if idx_start > 20
            prev_exc = exc_mem(max(1, idx_start-150) : idx_start-1);
            if length(prev_exc) > L_frame && ~any(isnan(prev_exc))
                 [c, lags] = xcorr(residual, prev_exc);
                 valid_mask = (lags >= 20) & (lags <= 140);
                 if any(valid_mask)
                     [~, max_idx] = max(abs(c(valid_mask)));
                     valid_lags = lags(valid_mask);
                     best_lag = abs(valid_lags(max_idx));
                 end
            end
        end
        
        % Construct Adaptive Vector
        if idx_start > best_lag
            past_idx = idx_start - best_lag;
            adapt_vec = exc_mem(past_idx : past_idx + L_frame - 1);
        else
            adapt_vec = zeros(L_frame, 1);
        end
        
        % Calculate Pitch Gain
        num_g = sum(residual .* adapt_vec);
        den_g = sum(adapt_vec.^2) + 1e-6;
        g_pitch = max(-1.5, min(1.5, num_g / den_g));
        
        % 4. Fixed Codebook (Residual Quantization)
        target_fixed = residual - (g_pitch * adapt_vec);
        
        g_fixed = mean(abs(target_fixed)); 
        if isnan(g_fixed), g_fixed = 0; end % Safety
        
        fixed_signs = sign(target_fixed);
        fixed_signs(fixed_signs == 0) = 1;
        
        fixed_vec = fixed_signs * g_fixed;
        
        % 5. Update Excitation
        exc_curr = (g_pitch * adapt_vec) + fixed_vec;
        
        % [FIX ADDED] Sanitize excitation to prevent error propagation
        exc_curr(isnan(exc_curr)) = 0; 
        
        exc_mem(idx_start:idx_end) = exc_curr;
        
        % 6. Synthesis Filter
        s_synth_frame = filter(1, A, exc_curr);
        s_synth(idx_start:idx_end) = s_synth_frame;
        
        % 7. Bit Counting
        b_lpc = 30;     
        b_pitch = 7+5;  
        b_fixed_g = 5;  
        b_fixed_v = L_frame * 1; 
        total_bits = total_bits + (b_lpc + b_pitch + b_fixed_g + b_fixed_v);
    end
    
    % De-emphasis
    s_synth = filter(1, [1 -preemp], s_synth);
    
    % [FIX ADDED] Final safety check for NaNs in output
    s_synth(isnan(s_synth)) = 0;
    
    stats.total_bits = total_bits;
end

function snr_val = calc_segsnr(sig, ref, fs)
    % Calculates Segmental SNR
    len = min(length(sig), length(ref));
    sig = sig(1:len); ref = ref(1:len);
    
    % Handle NaNs in inputs
    sig(isnan(sig)) = 0;
    ref(isnan(ref)) = 0;
    
    N = floor(0.02 * fs); 
    n_frames = floor(len/N);
    acc_snr = 0;
    count = 0;
    
    for k = 1:n_frames
        idx = (k-1)*N + 1 : k*N;
        s_frame = sig(idx);
        noise = ref(idx) - s_frame;
        
        p_s = sum(s_frame.^2);
        p_n = sum(noise.^2);
        
        if p_n > 1e-9 && p_s > 1e-9
            val = 10*log10(p_s / p_n);
            val = max(-10, min(35, val)); 
            acc_snr = acc_snr + val;
            count = count + 1;
        end
    end
    if count > 0, snr_val = acc_snr/count; else, snr_val = 0; end
end

function lsd = calc_lsd(orig, synth, fs)
    % Calculates Log-Spectral Distance
    len = min(length(orig), length(synth));
    x = orig(1:len); y = synth(1:len);
    
    % [FIX ADDED] Fix for "Expected X to be non-NaN" error
    x(isnan(x)) = 0;
    y(isnan(y)) = 0;
    
    nfft = 256; window = hanning(nfft); overlap = floor(nfft/2);
    
    [Sx, ~] = spectrogram(x, window, overlap, nfft, fs);
    [Sy, ~] = spectrogram(y, window, overlap, nfft, fs);
    
    Px = abs(Sx).^2 + 1e-10;
    Py = abs(Sy).^2 + 1e-10;
    
    diff = 10*log10(Px) - 10*log10(Py);
    lsd = mean(sqrt(mean(diff.^2, 1)));
end