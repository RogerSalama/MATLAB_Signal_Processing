# Signal Processing Project: Multi-Tone Analysis and Filtering

## Overview
This MATLAB project demonstrates signal generation, frequency analysis, and filtering techniques. It creates a multi-tone signal with frequencies at 500Hz, 1000Hz, 1500Hz, and 2000Hz, then processes it through low-pass and high-pass Butterworth filters. Key features include:
- Time-domain signal visualization
- Frequency spectrum analysis using FFT
- Energy calculation and Parseval's theorem verification
- Filter design and application
- Audio file generation

## Prerequisites
- MATLAB (R2019a or later recommended)
- Signal Processing Toolbox

## Files
1. `SignalsProject.m` - Main MATLAB script
2. `signal_1.wav` - Original multi-tone signal
3. `signal_2.wav` - Low-pass filtered output
4. `signal_3.wav` - High-pass filtered output
5. `README.md` - This documentation

## How to Run
1. Clone the repository:
2. Open MATLAB and navigate to the project directory
3. Run the main script:


## Code Structure
```matlab
%% Parameters
f1 = 500; f2 = 1000; f3 = 1500; f4 = 2000;  % Signal frequencies
fs = 10000;  % Sampling frequency
t = 0:1/fs:10;  % Time vector

%% Signal Generation
x = cos(2*pi*f1*t) + cos(2*pi*f2*t) + cos(2*pi*f3*t) + cos(2*pi*f4*t);

%% File Operations
audiowrite('signal_1.wav', x/4, fs);  % Save normalized signal

%% Visualization
plot(t(1:100), x(1:100));  % Plot first 10ms

%% Frequency Analysis
X_f = fft(x(1:round(1/gcd([f1,f2,f3,f4])*fs));  % FFT of one period
f = (-fs/2):(fs/length(X_f)):(fs/2 - fs/length(X_f));
plot(f, abs(fftshift(X_f)/length(X_f)));

%% Energy Verification (Parseval's Theorem)
energy_time = sum(x.^2)/length(x);
energy_freq = sum(abs(fftshift(X_f)/length(X_f)).^2);

%% Filter Design (Low-pass)
fc = 1250; order = 20;
[B, A] = butter(order, fc/(fs/2), 'low');
y = filter(B, A, x);
audiowrite('signal_2.wav', y/4, fs);

%% Filter Design (High-pass)
[B1, A1] = butter(order, fc/(fs/2), 'high');
y2 = filter(B1, A1, x);
audiowrite('signal_3.wav', y2/4, fs);
```

## Expected Outputs
1. **Time-domain Plots**:
   - Original multi-tone signal
   - Low-pass filtered signal (y1)
   - High-pass filtered signal (y2)

2. **Frequency-domain Plots**:
   - Original signal spectrum (4 peaks at 500, 1000, 1500, 2000 Hz)
   - y1 spectrum (only ≤1250 Hz components)
   - y2 spectrum (only >1250 Hz components)

3. **Energy Verification**:
```
X Energy in Time Domain: 2.00000
X Energy in Frequency Domain: 2.00000
Parseval's theorem is verified!

Y1 Energy in Time Domain: 1.00000
Y1 Energy in Frequency Domain: 1.00000
Parseval's theorem is verified!

Y2 Energy in Time Domain: 1.00000
Y2 Energy in Frequency Domain: 1.00000
Parseval's theorem is verified!
```

4. **Audio Files**:
   - `signal_1.wav`: Original composite signal
   - `signal_2.wav`: Low-pass filtered (bass-enhanced)
   - `signal_3.wav`: High-pass filtered (treble-enhanced)
  
## Plots

![Generated Signal](https://github.com/user-attachments/assets/34f49423-ea69-4514-9d25-965154d9357a)
![X(f)](https://github.com/user-attachments/assets/16e0af7b-f00c-49a4-9d7a-b7bba9e845ca)
![y1 Response](https://github.com/user-attachments/assets/fc08b871-2a10-4262-a97b-baf8f7bdade4)
![y1 Signal](https://github.com/user-attachments/assets/b310ec86-e2d3-41f0-a43a-c825f80954db)
![y1(f)](https://github.com/user-attachments/assets/eaf1d521-a8b4-4155-9784-65010dc5c6ac)
![y2 Response](https://github.com/user-attachments/assets/a3df8582-eec6-471c-b0cf-755bafefce94)
![y2 Signal](https://github.com/user-attachments/assets/ce06d53f-03b6-487e-a839-1517857ef06a)
![y2(f)](https://github.com/user-attachments/assets/d4f72dea-df93-4c7a-8e79-205d54cd009a)


## Key Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| `fs` | 10000 Hz | Sampling frequency |
| `fc` | 1250 Hz | Filter cutoff frequency |
| Order | 20 | Filter order |
| Duration | 10 sec | Signal length |

## Technical Notes
1. **Normalization**: Signals are divided by 4 before saving to prevent clipping
2. **Period Selection**: Fundamental period calculated using GCD of frequencies (500Hz period)
3. **FFT Normalization**: Spectra scaled by `1/N` for accurate magnitude representation
4. **Energy Calculation**: Time-domain energy normalized by signal length

