% interpolacja - funkcja interp1
clc; clear all;

signal = importdata("data/nsr001.dat"); 


plot(signal);

samples = 128;
f0 = 0.1;
t_k = [0, cumsum(rand(1, samples-1))]; % jest random, czyli jest niejednorodne

x = sin(2*pi*f0*t_k);
tn = linspace(0, max(t_k),samples);
fs = max(t_k)/samples;

%plot(x,y);


%y = sin(2*PI*fs*t_k)

% zaimplementować metodę Lomb-Scargle
x1 = interp1(t_k,x,tn);
figure(2);
plot(x1);
figure(4);
periodogram(x1);

% interp1
% periodogram

F = linspace(0,fs/2, 128);



F1 = []
for f=1:length(F)
   F1(end +1) = pls(f, t_k, x);
end
figure(3);
plot(F1);


function output = pf(f,t, x)

    part1 = sum(x .* cos(2*pi*f*t)).^2;
    part2 = sum(x .* sin(2*pi*f*t)).^2;
    output = (part1 + part2)/length(x);
end


function output = pls(f,t_k, x)
    tau = 1/(4*pi*f) * atan(sum(sin(4*pi*f*t_k))/sum(cos(4*pi*f*t_k)));

    part1 = sum(x .* cos(2*pi*f*(t_k - tau))).^2;
    part2 = sum(cos(2*pi*f*(t_k - tau)).^2);

    part3 = sum(x .* sin(2*pi*f*(t_k - tau))).^2;
    part4 = sum(sin(2*pi*f*(t_k - tau)).^2);

    output = 0.5 * ( part1/part2 + part3/part4 );

end
