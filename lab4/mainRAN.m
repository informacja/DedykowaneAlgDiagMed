clear all
% dane sytetyczne, outlayers

fns = ["data1.mat" "data2.mat" "data3.mat"];

for filename = fns
    d1 = importdata(filename);
    zadanyP = 0.5;
    figure, syntetic
end

% dane rzeczywiste 
data = load("eng_signals.mat");
data1 = data.eng_signal1;
data2 = data.eng_signal2;
data3 = data.eng_signal3;

y = data1; filename = "eng_signal1";
p2p = peak2peak(y);
mp = p2p/2;
figure, findpeaks(y,'MinPeakProminence',mp,'MinPeakDistance',100,"Annotate","peaks")
[pks,locs,width,prominence] = findpeaks(y,'MinPeakProminence',mp,'MinPeakDistance',100,"Annotate","peaks",'Annotate','extents');
realData
 ynn = polyval([suma 4000],[238 353]);
 plot([238 353], ynn,'c');

y = data2; filename = "eng_signal2";
p2p = peak2peak(y);
mp = p2p/3;
figure, findpeaks(y,'MinPeakProminence',mp,'MinPeakDistance',100,"Annotate","peaks")
[pks,locs,width,prominence] = findpeaks(y,'MinPeakProminence',mp,'MinPeakDistance',100,"Annotate","peaks",'Annotate','extents');
realData

y = data3; filename = "eng_signal3";
p2p = peak2peak(y);
mp = p2p/3;
figure, findpeaks(y,'MinPeakProminence',mp,'MinPeakDistance',100,"Annotate","peaks")
[pks,locs,width,prominence] = findpeaks(y,'MinPeakProminence',mp,'MinPeakDistance',100,"Annotate","peaks",'Annotate','extents');
realData
xlabel({"Średni wsp. kierunkowy = " num2str(-suma)})


h2 = findall(groot,'Type','figure');
h3 = findobj('Type','figure');

urlwrite ('https://raw.githubusercontent.com/informacja/MTF/main/figPW.m', 'figPW.m');

for i = 1:length(h2)
    figure(h2(i))
    figPW;
end





return
tabindexes = 545:661; 
data = data(tabindexes ,:);
datax = tabindexes;
datax = datax';

pos = randi(length(data));
n1 = data(pos,:);
n1 = [pos, n1]

pos2 = randi(length(data));
n2 = data(pos2,:);
n2 = [pos2, n2]

n = [n1; n2];

% plot(data ,'o');

fitLineFcn = polyfit(n(:,1),n(:,2),1); % fit function using polyfit

hold on 
plot(n1(:,1), n1(:,2), 'o', Color="red")

plot(n2(:,1), n2(:,2), 'o', Color="green")
hold on

p = polyval(fitLineFcn, datax)

% plot(datax,p)
% 
% hold on

wsp = fitLineFcn(:,2);
wsp_gorny = wsp+0.1*wsp;

wsp_gorny = [fitLineFcn(:,1); wsp_gorny];

p_gorny = polyval(wsp_gorny, datax)

% plot(datax,p_gorny)
% hold on


wsp = fitLineFcn(:,2);
wsp_dolny = wsp-0.1*wsp;


wsp_dolny = [fitLineFcn(:,1); wsp_dolny];

p_dolny = polyval(wsp_dolny, datax)

% plot(datax,p_dolny)
% hold on

%punkty mniejsze od gory i wieksze od dolnej

to = (data < p_gorny) & (data > p_dolny);

wynik = sum(to);

pky = [];
pkx = [];


for i = 1:length(data)
    if to(i) == 1
        pkt = data(i);
        pktx = datax(i);

        pky = [pky,pkt];
        pkx = [pkx, pktx];
    end
end

plot(pkx , pky, 'o');

hold on

wspolczynniki = polyfit(pkx,pky,1); 
p2 = polyval(wspolczynniki, pkx)
hold on
plot(pkx,p2)
