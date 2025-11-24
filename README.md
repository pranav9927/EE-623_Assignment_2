# EE-623 Assignment 2 – Speech Coding
**Sriram Pranav Gumpalli** | Roll No: 220102038

## Repository Structure
├── CELP_codec/               # CELP @ ~12.7 kbps
├── Plain_LPC_Vocoder/        # Plain LPC @ 3.1 kbps
├── Voice_Excited_Vocoder/    # Voice-Excited LPC @ 15.55 kbps
├── female1.wav, female2.wav
├── male1.wav, male2.wav
└── README.md

## Objective 1 – Wideband LPC Vocoders (16 kHz)

| Vocoder                | Bitrate   | Excitation                  | Avg SegSNR | Avg Time |
|-----------------------|-----------|-----------------------------|------------|----------|
| Plain LPC (LPC-10)    | 3.1 kbps  | Impulse train / Noise       | -2.21 dB   | 0.089 s  |
| Voice-Excited LPC     | 15.55 kbps| DCT-compressed residual (32 coeffs × 8 bit) | 7.12 dB | 0.155 s  |

- Frame: 30 ms, shift 20 ms, LPC order 16  
- Voice-Excited LPC gives ~10 dB better quality for 5× bitrate

## Objective 2 – CELP Codec (Narrowband 8 kHz)

- Target: ~14 kbps → achieved **12.7 kbps**  
- Frame: 10 ms (100 frames/s), 127 bits/frame  
- Avg SegSNR: **3.12 dB**  
- Avg LSD: **9.49 dB** (fair-to-good quality)  
- Avg encode/decode time: **0.07 s**
