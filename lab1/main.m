 clear all, close all
d = importdata("data/100_MLII.dat");
% d = importdata("data/228_MLII.dat");
x = d(720*55:720*56);% only for preview
x = d;

    M = 59;
    C = 24;

Fs = 360;
DADM_Lab1 % pierwszy podpunkt
% x = v;
% figure(1),plot(x)

T = 1/Fs;

L= length(x);
t = (0:(L-1))*T;
f = Fs*(0:(L/2))/L;

% win = hann(L);
for i = 1:L
    win(i)=(0.54-(0.46*cos((2*pi*i)/(N-1))));
end
xm = x - mean(x);
X = abs(fft(xm.*win'));
% [length(f),length(X)/2]

% X = fftshift(X);
figure, 
subplot(211); plot(f, X(1:L/2+1)/L); xlabel("Hz");
subplot(212);
% semilogx(f, 20*log10(X(1:L/2+1)./max(X(1:L/2+1))),'b'); title("semiLog"), xlabel("Hz");
semilogx(2:Fs/2, 20*log10(X(2:Fs/2)./max(X(2:Fs/2))),'b');
% figure(222), plot(X(L/2:end-1)); 

% powyższe przyda się na lab 2

% moving average to cure diff sharpnes
for i = 5:length(v)-5
    tmp = 0;
    for j = -4:5
        tmp = tmp + v(i+j);
    end
    g(i)= tmp;
end

%progowanie
xd = diff(g);
% figure, plot(xd); title("diff(g)");
xdd(length(xd)) = 0;

xd = xd+ max(xd)*0.2;
last = sign(xd(1));
for i = 2:length(xd)
    if(sign(xd(i)) ~= last)
        xdd(i) = 1;
    else
        xdd(i) = 0;
    end
    last = sign(i);
end
xddd = 0;
anomaly = [];
xddd(length(xddd)) = 0;
last = sign(xdd(1));
for i = 2:length(xdd)
    if(xdd(i)>last)
        if(x(i) < 0)
            xddd(i) = 1;
        else
            anomaly = [anomaly; i];
        end
    else
        xddd(i) = 0;
    end
    last = xdd(i);
end
% ps(i) = peaks_seconds;
t= [1:length(xdd)]/Fs;
figure, plot(t,xdd()); xlabel("sekundy")
   
    
% min(peaks_seconds)
sum(xddd)
QRSdetected = find(xddd);

t= [1:length(x)]/Fs;
 figure, plot(t,x); hold on;
 plot(t(QRSdetected),x(QRSdetected),'r*');
 xlabel("sekundy"); ylabel("mVolty"); title("EKG");
QRS
%  figure(3); plot(t,w,'r',t,yb,'k',t,y,'g'); axis('tight');legend('time','butter','MA');% odpowiedz skokowa   
% xlabel(sprintf('MSE(butter) = %g MSE(MA) = %g', MSEb, MSE));
FS=Fs;
order    =   5;
fcutlow  =   5;
fcuthigh =  25;
 
% onliczone wsp???‚czynnki filtra
[b1, a1] = butter(order,[fcutlow,fcuthigh]/(FS/2), 'bandpass');
%--------------------------------------------------------------------------
%Charakterystyka projektowanego filtra
[h,f]=freqz(b1,a1,2048,FS);
w = x;
N=length(w); t=[1:N];
%filtracja
yb=filter(b1,a1,w); 
M1=1; for(i=M1+1:N) y1(i)=mean(w(i-M1:i)); end,
M=32; for(i=M+1:N)  y(i) =mean(y1(i-M:i)); end;
y=y1-y;
figure
subplot(211); plot(yb); title("butter")
subplot(212); plot(y);  title("MA")
   
A=abs(h); fi=atan(imag(h)./real(h));
g1=[ones(1,M1)/M1 zeros(1,2048-M1)]; 
g=[ones(1,M)/M zeros(1,2*2048-M)]; 
h1=fft(g1); h1=h1(1:2048); A1=abs(h1); fi1=atan(imag(h1)./real(h1));
h2=fft(g);  h2=h2(1:2048); A2=abs(h2); fi2=atan(imag(h2)./real(h2));
fs=1./[1:2048]; As=A1-A2;
figure; subplot(2,1,1); plot(f, A,'m'); hold on; plot(fs,A1,'r',fs,A2,'k',f,As,'b'); hold off;
grid on; axis('tight'); legend("Butter",'',"Moving Average");
subplot(2,1,2); plot(2*f, fi/pi); grid on; axis('tight')


if (~isempty(anomaly))
    figure, plot(t, x, t(anomaly), x(anomaly), 'k*');
    if( length(anomaly) == 1)
        xlim([anomaly-360  anomaly+500]);
    end
    title("Anomalia");
end
% end
% figure,plot(indexy,peaks_seconds)
% figure, plot(xdd) 
return;
%Filtracja pasmowa
FS=2048;
order    =    5;
fcutlow  =   20;
fcuthigh =  400;
 
% onliczone wsp???‚czynnki filtra
[b1, a1] = butter(order,[fcutlow,fcuthigh]/(FS/2), 'bandpass');
%--------------------------------------------------------------------------
%Charakterystyka projektowanego filtra
[h,f]=freqz(b1,a1,2048,FS);
plot(f, 20*log10(abs(h)) );
set(gca, 'XScale', 'log')
ylim([-50 5])
xlim([20 40000])
grid on;
%--------------------------------------------------------------------------
X=zeros(FS, 1);
X(1)=1; %Delta kroneckera

% Y jest przefiltrowanym sygna?‚em X
Y = filter(b1, a1, X);    % Bandpass filter
%--------------------------------------------------------------------------
%Charakterystyka rzczywistego filtra na podstawie delty kronecker
figure(33)
Y=abs(fft(Y, FS));
semilogx(2:FS/2, 20*log10(Y(2:FS/2)./max(Y(2:FS/2))),'b');
xlim([0 30000])
ylim([-50 0])
grid on

A=abs(h); fi=atan(imag(h)./real(h));

 figure(2); subplot(2,1,1); plot(f, A); grid on; axis('tight'); 
 subplot(2,1,2); plot(f, fi/pi); grid on; axis('tight')
N=2100;t=[1:N];w=[zeros(1,100) randn(1,N-100)]; %ones(1,N-100)];

%  w=rawData(1,161:end); 
%  w = audioread('data\01108.wav'); w = w(:,1);
x = v;
 N=length(w); t=[1:N];
yb=filter(b1,a1,w); 
last = 1000;
% for j = 1:N/2
%   for k = 1:j/2  
    M1=1; for(i=M1+1:N) y1(i)=mean(w(i-M1:i)); end,
    M=32; for(i=M+1:N)  y(i)=mean(y1(i-M:i)); end; y=y1-y;
    MSEb = sum((w-yb).^2 );MSE= sum((w-y').^2);% wskaĹşnik bĹ‚Ä™u
%     figure(55), hold on; plot(j,MSEb,'b*',j,MSE,'g*'); hold off;
%     if(MSE < last) best = [k,j]; last = MSE; end    
%   end
% end
% best
figure(3); plot(t,w,'r',t,yb,'k',t,y,'g'); axis('tight');legend('time','butter','MA');% odpowiedz skokowa   
xlabel(sprintf('MSE(butter) = %g MSE(MA) = %g', MSEb, MSE));
g1=[ones(1,M1)/M1 zeros(1,2048-M1)]; g=[ones(1,M)/M zeros(1,2*2048-M)]; 
h1=fft(g1); h1=h1(1:2048); A1=abs(h1); fi1=atan(imag(h1)./real(h1)); h2=fft(g); h2=h2(1:2048); A2=abs(h2); fi2=atan(imag(h2)./real(h2));
fs=1./[1:2048]; As=A1-A2;
figure(2); subplot(2,1,1); plot(f, A,'m'); hold on; plot(fs,A1,'r',fs,A2,'k',f,As,'b'); hold off;
grid on; axis('tight'); subplot(2,1,2); plot(2*f, fi/pi); grid on; axis('tight')


