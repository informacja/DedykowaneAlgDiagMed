load abr_signal2.mat
data2 = [];
for (i = 1:length(abr_signal2)) 
%         data = [data audioread(strcat(folder, fnames(i)))];
        data2 = [data2 abr_signal2{i}.data+abr_signal2{i}.dB];
%         figpos = {'north','south','east','west','northeast','northwest','southeast','southwest','center','onscreen'};
%         movegui(cell2mat( figpos(mod(nfig+i,length(figpos)-1)+1) ));
end
% figure(1); plot(data2); title("Data 2");
Ewykl=2; 
Takcji=5000; T1sek=2048; dt=1/T1sek; TNyq=2/T1sek;
Lk=0; nfig0=2; txtau='\tau'; txDelta='\Delta'; 
[Ldsz Lszer] = size(data2); ldC = Ldsz;
nrwykr=1:Lszer; 
% ====== Parametry MTF ===============
% Tu okres skĹadowej cyklicznej - szerokoĹÄ pasma filtru MTF
% Takcji - liczba prĂłbek w jednej akcji 
LFg=2; LFd=2; fd=8; %fg=1/(Tud*dt); Tud=1/(fg*dt);
Tud=Takcji/(12); Tud=Takcji/(13); %Tud=Takcji/(10);
Tud1=Tud/1.2; %turedniania trendu
TuG=Tud*[1/2 1/3 1/10]; TuG=Tud*[1/10 1/12 1/10]; TuG=Tud*[1/5 1/6 1/10]; TuG=Tud*[1/4 1/5 1/10]; %13.5]; % LFg=2 bÄdÄ dwa filtry Tug
%Tud=1/(fd*dt);
TR=Tud;% Truchu=Takcji/12;
% ................ #Wybor typu filtru .......................
typMTF=[5 5 5 3]; %[3 3 5 3]; % typ filtru: 3 - Z1 - filtr koncowy zwykly; 5 - z3 - trend 3.go rzÄdu 
typFf=1; % - typ filtru konc.: 0 nowy, 1 zwykly

% =============== Parametry szergu i filtra ==============
LSzeregow=length(nrwykr); % liczba badanych szeregĂłw
%close all, 
Lk=0; Losob=Lszer/8; LKanalow=round(LSzeregow/Losob); 
nfig=nfig0; nFig=nfig0+1+LKanalow;
% ----------- Okrelenie zakresu lT analizy widm (tak aby dobrze spróbkować widmo)
lTx=Tud*200; lTx=Tud*500; %1000;
nTx=ceil(lTx/4096); lT=nTx*4096; %500; %350; %75; %round(lA/10);
Ewykl=2; 
ntypZ=typMTF(1); rzadB=5; 
tSynt=tic;
if(0) figF=1; % Liczenie z rysowaniem wszystkich charakterystyk
    %[bfBd, afBd, tauhPfd, MTFd, Fzwc, LpFc, Amp, nA01, WcB, Fzwd, Fzwd2, fi,  Ffc, Ff1, Ff2]=...
    [bfBd, afBd, tauhPfd, MTFd, Fzwc, LpFc, Amp, nA01, WcB]=...
        desMTFcButter(ntypZ,[Tud Tud1],rzadB,lT,figF,'Synteza filtrow MTF (Fzws_b Fzwd_r Fzw2_m Ffc_g Ff1_c Ff2_{c--}) i Butter_k','kbrmgc');
    Fzwd=[]; Fzwd2=[]; fi=[];  Ffc=[]; Ff1=[]; Ff2=[];
else % Tylko synteza filtrow Fzwc i Butter (szybka wersja)
    nA01=0;
    [bfBd, afBd, tauhPfd, MTFd, Fzwc, LpFc, Amp, nA01, WcB]=desMTFcButter(ntypZ,[Tud Tud1],rzadB,lT);
    %[bfB, afB, tauhPf, MTFd, Fzwc, LpFc, Amp, nA01, WcB, Fzwd, Fzwd2]=desMTFcButter(ntypZ,[Tud Tud1],rzadB,lT);
end
ihP(1)=WcB*lT/2; 
if(length(nA01)>1)  
    AmpB=zeros(4,length(Amp(1,:))); Ampd=zeros(4,length(Amp(1,:)));
    iA01(1)=nA01(2); AmpB(1,:)=Amp(1,:);  Ampd(1,:)=Amp(2,:); 
end; 

% -------------- Synteza filtrów dla energii ----------------
if(0) % Liczenie z rysowaniem wszystkich charakterystyk
    figF=200;
    [bfBE, afBE, tauhPfE, MTFE, FzwcE, LpFcE, Amp, nA01E, WcB, FzwdE, Fzwd2E, fi,  Ffc, Ff1, Ff2]=...
        desMTFcButter(ntypZ,[Tud Tud1]/Ewykl,rzadB,lT,figF,'Synteza filtrow MTF (Fzws_b Fzwd_r Fzw2_m Ffc_g Ff1_c Ff2_{c--}) i Butter_k','kbrmgc');
else % Tylko synteza filtrow Fzwc i Butter (szybka wersja)
    [bfBE, afBE, tauhPfE, MTFE, FzwcE, LpFcE, Amp, nA01E, WcB]=desMTFcButter(ntypZ,[Tud Tud1]/Ewykl,rzadB,lT);
end
% ............. Koniec syntezy filtrów dolnoprzepustowych MTFd i Butter5 ............
% ============= Synteza filtrów górnoprzepustowych % ========================   
%[bfx, afx, tauhPfEg, MTFEg, FzwcEg, FzwgE, Fzwg2E, LpFcEg, Amg, nA01g, WcBg]=desMTFcButter(ntypZ,[Tud Tud1]/Ewykl); %,0,lT);
[bfEg, afEg, tauhPfEg, MTFEg, FzwcEg, LpFcEg, Amp, nA01Eg,WcEg]=desMTFcButter(ntypZ,TuG(1:2)/Ewykl,rzadB,lT);
ihP(4)=WcEg*lT/2; TuGg(2)=TuG(2)/Ewykl;
if(length(nA01Eg)>1)  iA01(4)=nA01Eg(2); AmpB(4,:)=Amp(1,:); Ampd(4,:)=Amp(2,:); end; %wPol=2*ihP(2)/lT;
[bfg, afg, tauhPfg, MTFg, Fzwcg, LpFcg, Amp, nA01g, Wcg]=desMTFcButter(ntypZ,TuG(1:2),rzadB,lT);
ihP(3)=Wcg*lT/2; TuGg(1)=TuG(2);
if(length(nA01g)>1)  iA01(3)=nA01g(2); AmpB(3,:)=Amp(1,:);  Ampd(3,:)=Amp(2,:); end, 
MTud=MTFd(1).M; 
% ======= Koniec syntezy - podstawiamy wyniki ==================
nA01=max(iA01); lwAm=nA01; Amp=AmpB(1:4,1:lwAm); AmpB=Amp; Amp=[]; Amp=Ampd(1:4,1:lwAm); Ampd=Amp; Amp=[]; 
% ======================================================================
afB(1,:)=afBd;  bfB(1,:)=bfBd; 
Lzwc(1)=LpFc(2); Lzwd=LpFc(3); Lzw2=LpFc(4); % długości filtrów Fzwc, Fzwd i Fzw2
tauhP(1)=-tauhPfd(1); wPol=2*ihP(1)/lT;
% .................................................................
afB(2,:)=afBE; bfB(2,:)=bfBE;
Lzwc(2)=LpFcE(2); LzwdE=LpFcE(3); Lzw2E=LpFcE(4); 
tauhP(2)=-tauhPfE(1); %ihPB(2)=-tauhPfE(1);% ihP(2)=-tauhPfE(2);
% .................................................................
afB(3,:)=afg; bfB(3,:)=bfg;
Lzwc(3)=LpFcg(2); Lzwdg=LpFcg(3); Lzw2g=LpFcg(4); 
tauhP(3)=-tauhPfg(1); %ihPB(3)=-tauhPfg(1); ihP(3)=-tauhPfg(2);
% .................................................................
afB(4,:)=afEg; bfB(4,:)=bfEg;
Lzwc(4)=LpFcEg(2); LzwdEg=LpFcEg(3); Lzw2Eg=LpFcEg(4); 
tauhP(4)=-tauhPfEg(1); %ihPB(4)=-tauhPfEg(1); ihP(4)=-tauhPfEg(2);
if(length(nA01Eg)>1)  iA01(4)=nA01Eg(2); end, %AmpB(2,:)=Amp(1,:); end; 
% ================= Synteza filtru wygładazającego MTF ============
lfA=lT/MTud; Tu=Tud/5; lfA=lT/Tud;
Tu=Tud/MTud*lfA; % PrzedziaĹ‚ wygĹ‚adzania widma
for(nfw=1:2)
    [M, Fzw]=MTFdesign(ntypZ, Tu); %[M, Fzw, Ff, F0]=MTFdesign(ntypZ, Tu);
    MTF(nfw).Tu=Tu; MTF(nfw).M=M; MTF(nfw).Fzw=Fzw; MTF(nfw).F0=[]; MTF(nfw).Ff=[];
    Tu=6;
end
[M, Fzwg0]=MTFdesign(ntypZ, TuG(3)); 
fprintf(1,'\nCzas syntezy %.4f sek.',toc(tSynt));
% ========================================================
typS2=typS+2;
            nYTr1(typS)=(Lzwc(typS)-1)/2;        m=nYTr1(typS);
            nYTr1(typS2)=(Lzwc(typS2)-1)/2+m;   mg=nYTr1(typS2); %Lzwc(typS); 
            yTrc(typS2,Nf)=0; yTrc(typS,Nf)=0;
Lzw=Lzwc(typS); Lzw1=Lzw-1;
            Lzwg=Lzwc(typS2); Lzwg1=Lzwg-1; Lzwg=Lzwg+nYTr1(typS);
            nYTr1(typS)=nYTr1(typS)+1; nYTr1(typS2)=nYTr1(typS2)+1;

Nd=MTFE(1).M-1; Ng=MTFEg(1).M-1; 
for(n=Lzw:Nf) m=m+1; yTrc(typS,m)=Yoryg(n-Lzw1:n)'*FzwcE; end
Yreszt=Yoryg-yTrc(typS,:)';
for(n=Lzwg:m) mg=mg+1; yTrc(typS2,mg)=Yreszt(n-Lzwg1:n)'*FzwcEg; end

nYTrf(typS)=m; nYTrf(typS2)=mg;