x = importdata("data/100_MLII.dat");

figure(1),plot(x)

Fs = 360;
T = 1/Fs;

L= 512;
t = (0:(L-1))*T;

win = hann(L);
xm = x - mean(x);
X = fft(xm(1:L).*,Fs);

X = fftshift(X/L);

% figure(2), plot(t,abs(X(L/2:end-1))); xlabel("Hz");

% powyższe przyda się na lab 2

DADM_Lab1 % pierwszy podpunkt

%progowanie
xd = diff(v);
figure, plot(xd); title("diff(v)");
xdd(length(xd)) = 0;

last = sign(xd(1));
for i = 2:length(xd)
    if(sign(xd(i)) ~= last)
        xdd(i) = 1;
    else
        xdd(i) = 0;
    end
    last = sign(i);
end

for i = 2:length(xdd)

end
sum(xdd)
figure, plot(t,xdd(1:512)), xlabel("sekundy")

figure, plot(xdd) 
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
figure
Y=abs(fft(Y, FS));
semilogx(2:FS/2, 20*log10(Y(2:FS/2)./max(Y(2:FS/2))),'b');
xlim([0 30000])
ylim([-50 0])
grid on

A=abs(h); fi=atan(imag(h)./real(h));

% figure(2); subplot(2,1,1); plot(f, A); grid on; axis('tight'); 
% subplot(2,1,2); plot(f, fi/pi); grid on; axis('tight')
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
    figure(55), hold on; plot(j,MSEb,'b*',j,MSE,'g*'); hold off;
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


