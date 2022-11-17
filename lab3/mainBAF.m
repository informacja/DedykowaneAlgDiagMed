clear all; clc; close all;
filename = "data/nsr001.dat";
% filename = "data/chf206.dat";
signal = importdata(filename);  %odstępy czasowe

figure(1), plot(signal); title("Sygnał surowy")

% całkowanie
y = cumsum(signal - mean(signal));
x = cumsum(signal);
figure(2), plot(x,y); title("Scałkowany sygnał")

vec = 4:64;
len_signal = length(y);
n = 4; F = [];
to_miss = rem(len_signal , n); % reszta z dzielenia 
for n = vec 
    s = 0; 
    for i = 1: n : len_signal - ((n+1) - to_miss)
        % zebrac 4 probki
        c = polyfit(x(i:i+n),y(i:i+n), 1);
        yn = polyval(c, x(i:i+n));
        y2 = yn - y(i:i+n);
        s = s + sum(y2.^2);
    end
    F = [ F sqrt(s/len_signal) ];
end

% F
figure(3), plot(log10(vec), log10(F)); hold on; ylabel("Dopasowanie prostych w skali log10")
vec = 4:16;
a1 = polyfit(log10(vec), log10(F(vec)), 1)
plot(log10(vec), polyval(a1,log10(vec)))
vec = 17:61;
a2 = polyfit(log10(vec), log10(F(vec)), 1)
plot(log10(vec), polyval(a2,log10(vec)))
xlabel(sprintf("a1 = %f, a2 = %f (współczynniki kierunkowe)", a1(1), a2(1)));
title(filename);
return
% fluktuacje to roznice pomiędzy próbkami
% przedział n 4
% polyval
% kwadrat odległosci 

% 1:16 polyfit, 17:64
% wsp kiern prostej dla zdrowego i chorego

% PSW save all figures