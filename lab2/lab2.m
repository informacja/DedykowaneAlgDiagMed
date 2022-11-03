% interpolacja - funkcja interp1
clc; clear all;

signal = importdata("data/nsr001.dat"); % diff(R)
% signal = importdata("data/chf206.dat"); 
% figure(1), plot(signal); title("raw signal")

% 4)
samples = 128;
f0 = 5;
t_k = [0, cumsum(rand(1, samples-1))]; % jest random, czyli jest niejednorodne

x = sin(2*pi*f0*t_k);
tn = linspace(0, max(t_k), samples);
fs = max(t_k)/samples;

%plot(x,y);
% y = sin(2*pi*fs*t_k)

% zaimplementowano metodę Lomb-Scargle
x1 = interp1(t_k,x,tn);
% figure(2), plot(x1);
% 5)
figure(5), periodogram(x1);

% periodogram(x1,rectwin(length(x1)),length(x1),fs)
% plot(w,10*log10(pxx))

F = linspace(0,fs/2, 128);

L = length(F); F1 = []; F2 = []; F1(L) = 0; F2(L) = 0;
for f=1:L
   F1(f) = pls(f, t_k, x);
   F2(f) = pf(f, t_k, x);
end
figure(3), plot(F, F1);
figure(7), plot(F, F2);

function output = pf(f, t, x) %power of frequency

    part1 = 0;
    k = 1;
    for i = t
        s = x(k) * cos(2*pi*f*i)^2;
        part1 = s+part1;
        k = k+1;
    end
    k = 1;
    part2 = 0;
    i = 0;
    for i = t
        s = x(k) * sin(2*pi*f*i)^2;
        part2 = s+part2;
    end

    part1 = sum(x .* cos(2*pi*f*t)).^2;
    part2 = sum(x .* sin(2*pi*f*t)).^2;
    output = (part1 + part2)/length(t);
end

function output = pls(f, t_k, x) %Periodogram Lomb-Scargle
    tau = 1/(4*pi*f) * atan(sum(sin(4*pi*f*t_k))/sum(cos(4*pi*f*t_k)));

    part1 = sum(x .* cos(2*pi*f*(t_k - tau))).^2;
    part2 = sum(cos(2*pi*f*(t_k - tau)).^2);

    part3 = sum(x .* sin(2*pi*f*(t_k - tau))).^2;
    part4 = sum(sin(2*pi*f*(t_k - tau)).^2);

    output = 0.5 * ( part1/part2 + part3/part4 );

end
