% clear all;

d = importdata("data/100_MLII.dat");
x= d(1:5e3);% only for preview
p = 4;% subplots
M=50;
n=-M:1:M;
h_low=zeros(2*M+1,1);
h_high=zeros(2*M+1,1);

fs=360;
fc_lowpass=15;
fc_highpass=5;

fc_cut_lowpass=fc_lowpass/(fs/2);
fc_cut_highpass=fc_highpass/(fs/2);


for i= -M:1:M
   if i==0
       h_low(i+M+1)=2*fc_cut_lowpass;
   else 
       h_low(i+M+1)=sin(2*pi*fc_cut_lowpass*i)./(pi.*i);
   end
end

for i= -M:1:M
   if i==0
       h_high(i+M+1)=1-(2*fc_cut_highpass);
   else 
       h_high(i+M+1)=-(sin(2*pi*fc_cut_highpass*i)./(pi.*i));
   end
end

N=(M*2)+1;
w = hann(N);
% w=zeros(N,1);

% for i = 1:1:N
%     w(i)=(0.54-(0.46*cos((2*pi*i)/(N-1))));
% end
 
wynik_low=w.*h_low;
wynik_high=w.*h_high;

splot_low=conv(x,wynik_low);
splot_high=conv(wynik_high,x);

splot_final=conv(splot_low,wynik_high);
subplot(p,1,2),plot(splot_final)

subplot(p,1,1),plot(x)

x = splot_final;

% b) RÓŻNICZKOWANIE

for n = 3:length(x)-2
    y(n)=(1/8)*(-x(n-2)-(2*x(n-1))+(2*x(n+1))+x(n+2));
end
subplot(p,1,3), plot(y); title("Różniczkowanie")

% c) Potęgowanie

for n = 1:length(y)
    z(n)= y(n)^2;
end

% d) Całkowanie

C = 15;

for n = 2*C+1:length(x)-C
    tmp = 0;
    for c = 1:C
        tmp = tmp + z(n-C-c);
    end
    v(n) = tmp/C;
end

subplot(p,1,4), plot(v); title("Całkowanie")