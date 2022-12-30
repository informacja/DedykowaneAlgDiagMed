%{ Konsultacje
% amplitudy do tabelki
% sprawdzić czy niższe natężenia występują coraz później
%} 
%% Const
MARGINES = 95; % 1/2 przedziału [ms]
PROG_SLYSZALNOSCI = 5; % minimal peak2peak
fs = 1e5;
%% Load data
load abr_signal2.mat
data2 = []; data3 = []; data4 = [];
for (i = 1:length(abr_signal2))
    data2 = [data2 abr_signal2{i}.data+abr_signal2{i}.dB];
end
figure(1); plot(data2); title("Data 2");

load abr_signal3.mat
for (i = 1:length(abr_signal3)) 
    data3 = [data3 abr_signal3{i}.data+abr_signal3{i}.dB];
end
figure(2); plot(data3);  title("Data 3");        

load abr_signal4.mat
for (i = 1:length(abr_signal4)) 
    data4 = [data4 abr_signal4{i}.data+abr_signal4{i}.dB];
end
figure(3); plot(data4); title("Data 4");   

%% Data procesing
data = data2; fnr = 6;
procesing; title("abr_signal2");

data = data3; fnr = fnr + 1;
procesing; title("abr_signal3");

data = data4; fnr = fnr + 1; 
procesing; title("abr_signal4"); 

%% Figures
figpos = {'north','south','east','west','northeast','northwest','southeast','southwest','center','onscreen'};
h2 = findall(groot,'Type','figure');
for i = 1:length(h2)
    figure(h2(i))
    movegui(cell2mat( figpos(mod(i,length(figpos)-1)+1) ));
end

if ~isfile('figPSW.m')
    urlwrite ('https://raw.githubusercontent.com/informacja/MTF/main/figPSW.m', 'figPSW.m');
end
figPSW;