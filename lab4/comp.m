
clear all.close all
data = load("eng_signals.mat");
data = data.eng_signal1;
plot(data);
data = data(545:661,:);
datax = 545:661;
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
hold on
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

to = data < p_gorny & data > p_dolny;

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