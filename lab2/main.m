lab2;
k = 4; %l subplotów
% 1
zad1_sekund = 1/0.0033
Fs = 360;
sig_len_hours = sum(signal)/60/60

signal = signal(1:1e2);
% 2)
samples = length(signal);
x = cumsum(signal); 
y = signal;
figure(12), subplot(k,1,1), plot(x,y), axis("tight"), title("Raw signal");

tn = linspace(0, max(x), samples);
T = max(x)/samples;
fs = 1/T;
t = [0:samples-1]/fs;
x1 = interp1(x,y,tn);
subplot(k,1,2), plot(t, x1), title("Interpolated");
x1(1:2) = [0.7 0.703];
mustBeFinite(x1)
subplot(k,1,3), periodogram(x1);

F = linspace(0,fs/2, length(signal));

L = length(F)
F1 = []
F1(L) = 0;
for f=1:L
   F1(f) = pls(f, x, signal);
end
subplot(k,1,4), plot(F1);

function output = pf(f,t, x)
    part1 = sum(x .* cos(2*pi*f*t)).^2;
    part2 = sum(x .* sin(2*pi*f*t)).^2;
    output = (part1 + part2)/length(x);
end

