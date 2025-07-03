% Parameters
f1 = 500;
f2 = 1000;
f3 = 1500;
f4 = 2000;

fs = 10000;  % Sampling frequency
ts = 1/fs;
t = 0:ts:10;  % Time vector
t_plot = 0:ts:0.01;  % Short time vector for plotting
x = cos(2*pi*f1*t) + cos(2*pi*f2*t) + cos(2*pi*f3*t) + cos(2*pi*f4*t);
x_plot = cos(2*pi*f1*t_plot) + cos(2*pi*f2*t_plot) + cos(2*pi*f3*t_plot) + cos(2*pi*f4*t_plot);



% Step 2: Save the signal to a WAV file
filename = 'signal_1.wav';
x = x / 4;  % Normalize to prevent clipping
audiowrite(filename, x, fs);  % Save the normalized signal
%sound(x, fs);

x = x * 4;  % Restore original amplitude

% Step 3: plot the signal
figure;
plot(t_plot, x_plot);
title('Generated Signal');
xlabel('Time (s)');
ylabel('Amplitude');

% Step 4: Compute the energy of the signal in the time domain
energy_time = sum(x.^2) / length(x);

%Extract one fundamental period of the signal
gcd_freq = gcd(gcd(f1, f2), gcd(f3, f4));  % GCD of frequencies
T_fundamental = 1 / gcd_freq;  % Fundamental period
N0 = round(T_fundamental * fs);  % Number of samples in one period
x_periodic = x(1:N0);  % Extract one period

% Step 5: Compute the frequency spectrum using the FFT
X_f = fft(x_periodic);

% Frequency vector
f = (-fs/2):(fs/length(x_periodic)):(fs/2 - fs/length(x_periodic));

% Shift the spectrum to center the zero frequency component
X_f_shifted = fftshift(X_f) / length(x_periodic);

% Step 6: Plot the magnitude of the frequency spectrum
figure;
stem(f, abs(X_f_shifted), 'marker', 'none');
title('Magnitude of Frequency Spectrum X(f)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Step 7:
% Compute the energy in the frequency domain
energy_freq = sum(abs(X_f_shifted).^2);

% Verify Parseval's theorem
fprintf('X Energy in Time Domain: %.5f\n', energy_time);
fprintf('X Energy in Frequency Domain: %.5f\n', energy_freq);

if abs(energy_time - energy_freq) < 0.001
    disp('Parseval’s theorem is verified!');
else
    disp('Parseval’s theorem is NOT verified.');
end

%===================================================
% Step 8: Designing a butterworth lpf:
fc = 1250 %cut-off freq
Wn = fc / (fs / 2);
% normalize cut-off w.r.t to half samlpling freq
order = 20;
%here I design the filter
[B, A] = butter(order, Wn, 'low');


% Frequency response
[H, w] = freqz(B, A, 4000, fs); % Increase points for better resolution

% Step 9:
% Magnitude plot
figure;
subplot(2, 1, 1);
plot(w, abs(H));  % Plot the magnitude of H
title('y1 Magnitude Response of the Filter');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Phase plot without unwrapping to see the raw phase
subplot(2, 1, 2);
plot(w, unwrap(angle(H)));  % Plot the phase of H (in radians)
title('y1 Phase Response of the Filter');
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
grid on;

% Step 10:
y = filter(B, A, x);
y_plot = filter(B, A, x_plot);

% Step 11: Save the signal to a WAV file
filename = 'signal_2.wav';
y = y / 4;  % Normalize to prevent clipping
audiowrite(filename, y, fs);  % Save without normalization to preserve amplitude
%sound(y, fs);

y=y*4;

% Step 12: Plot the original and filtered signal
figure;
subplot(2, 1, 1);
plot(t_plot, x_plot);
title('y1 Original Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(t_plot, y_plot);
title('y1 Filtered Signal');
xlabel('Time (s)');
ylabel('Amplitude');



% Step 13: Energy of y1:
energy_y1_time = sum(y.^2) / length(y);

%Fundamental Period Calculation
gcd_freq = gcd(gcd(f1, f2), gcd(f3, f4));  % GCD of frequencies
T_fundamental = 1 / gcd_freq;  % Fundamental period
N0 = round(T_fundamental *fs);  % Number of samples in one period
y_periodic = y(10*N0:11*N0-1);  % Extract one period
%======================

% Step 14: Compute the frequency spectrum using the FFT
Y_f = fft(y_periodic);

% Frequency vector
f = (-fs/2):(fs/length(y_periodic)):(fs/2 - fs/length(y_periodic));

% Shift the spectrum to center the zero frequency component
Y_f_shifted = fftshift(Y_f)/length(y_periodic);

% Step 15: Plot the magnitude of the frequency spectrum
figure;
stem(f, abs(Y_f_shifted), 'marker', 'none');
title('Magnitude of Frequency Spectrum y1(f)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Step 16:
% Compute the energy in the frequency domain
energy_y1_freq = sum(abs(Y_f_shifted).^2);

% Verify Parseval's theorem
fprintf('Y1 Energy in Time Domain: %.5f\n', energy_y1_time);
fprintf('Y1 Energy in Frequency Domain: %.5f\n', energy_y1_freq);

if abs(energy_y1_time - energy_y1_freq) < 0.001
    disp('Parseval’s theorem is verified!');
else
    disp('Parseval’s theorem is NOT verified.');
end

% Step 17:Designing a butterworth :
fc = 1250 %cut-off freq
Wn = fc / (fs / 2);
% normalize cut-off w.r.t to half samlpling freq
order = 20;
%here I design the filter
[B1, A1] = butter(order, Wn, 'high');

% Frequency response
[H1, w1] = freqz(B1, A1, 4000, fs); % Increase points for better resolution

%Step 18:
% Magnitude plot
figure;
subplot(2, 1, 1);
plot(w1, abs(H1));  % Plot the magnitude of H
title('y2 Magnitude Response of the Filter');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Phase plot without unwrapping to see the raw phase
subplot(2, 1, 2);
plot(w1, unwrap(angle(H1)));  % Plot the phase of H (in radians)
title('y2 Phase Response of the Filter');
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
grid on;

% Step 19:
y2 = filter(B1, A1, x);
y2_plot = filter(B1, A1, x_plot);

% Step 20: Save the signal to a WAV file
filename = 'signal_3.wav';
y2 = y2 / 4;  % Normalize to prevent clipping
audiowrite(filename, y2, fs);  % Save without normalization to preserve amplitude
%sound(y2, fs);

y2=y2*4;

% Step 21: Plot the original and filtered signal
figure;
subplot(2, 1, 1);
plot(t_plot, x_plot);
title('y2 Original Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(t_plot, y2_plot);
title('y2 Filtered Signal');
xlabel('Time (s)');
ylabel('Amplitude');

% Step 22: Energy of y2:
energy_y2_time = sum(y2.^2) / length(y2);

%Fundamental Period Calculation
gcd_freq = gcd(gcd(f1, f2), gcd(f3, f4));  % GCD of frequencies
T_fundamental = 1 / gcd_freq;  % Fundamental period
N0 = round(T_fundamental *fs);  % Number of samples in one period
y2_periodic = y2(10*N0:11*N0-1);  % Extract one period

% Step 23: Compute the frequency spectrum using the FFT
Y2_f = fft(y2_periodic);

% Frequency vector
f = (-fs/2):(fs/length(y2_periodic)):(fs/2 - fs/length(y2_periodic));

% Shift the spectrum to center the zero frequency component
Y2_f_shifted = fftshift(Y2_f)/length(y2_periodic);

% Step 24: Plot the magnitude of the frequency spectrum
figure;
stem(f, abs(Y2_f_shifted), 'marker', 'none');
title('Magnitude of Frequency Spectrum y2(f)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Step 25:
% Compute the energy in the frequency domain
energy_y2_freq = sum(abs(Y2_f_shifted).^2);

% Verify Parseval's theorem
fprintf('Y2 Energy in Time Domain: %.5f\n', energy_y2_time);
fprintf('Y2 Energy in Frequency Domain: %.5f\n', energy_y2_freq);

if abs(energy_y2_time - energy_y2_freq) < 0.001
    disp('Parseval’s theorem is verified!');
else
    disp('Parseval’s theorem is NOT verified.');
end
