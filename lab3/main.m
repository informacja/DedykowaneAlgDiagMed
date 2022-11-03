clear all;
signal = importdata("data/chf206.dat");

figure(1), plot(signal);

% n - [4;64]
signal = signal - mean(signal);
figure(6), plot(abs(fft(signal)));

x = signal;
L = length(x);
L=2000;
y(L) = 0;
s = 0;
for k = 1:L



    
    n = floor(L/k);
    Nfloor = n;
    D = x(1:k);
    
    t = cumsum(D-mean(D));


%     tmp = 0;
%     s = s + x(k);
%     for i = 1:k
%         tmp = tmp + x(i);
% %         y(k) = s - 30;
%     end
%     y(k) = t ;
end

figure(2), plot(t)