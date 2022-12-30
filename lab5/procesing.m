
figure(10+fnr); hold on
subplot(411); plot(data(:,1)); title("raw");
subplot(412); plot(x); title("x (filtered)");
subplot(413); plot(v); title("v (wavelet)");
subplot(414);  hold on; plot(A); plot(przejsciaPrzezZero,A(przejsciaPrzezZero), 'g^'); title("A (przejsciaPrzezZero)");
hold off;
tabelka = ["fala V-ta","III-cia", "I-wsza"];
for i = 1:size(data,2)
    signal = data(:,i);
    t = 0:length(signal)/fs;
    % butter
%     x = bandpass(signal,[100 3000],fs);
    order = 4;
%     [b,a]=butter(order,[100, 3000]/(fs/2),'bandpass'); 
    b = fir1(48,[100, 3000]/(fs/2));
    a = ones(size(b));
    x = filtfilt(b,1, signal);
%     x = signal;
    [SWA SWD] = swt(x, 3, "rbio3.1");
    
    figure(4), %plot(signal); hold on;
    plot(SWD(3,:)'-i*30); hold on; title("SWD"); hold off;
    
    SWD_ = circshift(SWD(3,:),1);
    
    A = SWD(3,:).*SWD_;
    v = SWD(3,:);
    v = [ zeros(1,100) v(100:950)];
    A = [ zeros(1,100) A(100:950)];
    przejsciaPrzezZero = find(A<0);
%     plot(przejsciaPrzezZero,v(przejsciaPrzezZero),'r*')
    %% Detecja I III V
    
    for i = 1:length(przejsciaPrzezZero)
        rank = x(przejsciaPrzezZero);
    end
    nowePrzejsciaPrzezZero = przejsciaPrzezZero;

    indexyFal(3) = 0;
    for i = 1:3
        [M,I] = max(rank); % znajdz szczyt
        indexyFal(i) = nowePrzejsciaPrzezZero(I);
        % zdefiniuj nowy przedział
        wieksze = find(nowePrzejsciaPrzezZero> nowePrzejsciaPrzezZero(I)+MARGINES);% marginesy
        mniejsze = find(nowePrzejsciaPrzezZero< nowePrzejsciaPrzezZero(I)-MARGINES);% marginesy
        nowePrzejsciaPrzezZero = [nowePrzejsciaPrzezZero(mniejsze) nowePrzejsciaPrzezZero(wieksze)];
        rank = [];
        for i = 1:length(nowePrzejsciaPrzezZero) % ponownie policz ranking
            rank = x(nowePrzejsciaPrzezZero);
        end
    end
    % próg słyszalności
    figure(7+fnr), plot(x); hold on; xlabel("ms")
    p2p = peak2peak(x(indexyFal));
    if(p2p >= PROG_SLYSZALNOSCI)
        plot(indexyFal, x(indexyFal),'m*')
        tabelka = [ x(indexyFal)'; tabelka];
    end
%     
%     p2p = peak2peak(A);
%     mp = p2p/100;
% 
%     p2p = peak2peak(A);
%     mp = p2p/100;
%     figure(fnr), title("FindPeaks(x)"); hold on; findpeaks(x-i*30,'MinPeakProminence',mp,'MinPeakDistance',100,"Annotate","peaks")
%     [pks,locs,width,prominence] = findpeaks(x,'MinPeakProminence',mp,'MinPeakDistance',100,"Annotate","peaks");
%     pks;

%     indexyFal(3) = 0;
end
tabelka