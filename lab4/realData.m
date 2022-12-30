d = diff(locs);
s = 0;
hold on;
for i = 1:length(d)
    x = [locs(i) locs(i+1)];
    y = [pks(i) pks(i+1)-prominence(i)];
    fitLineFcn = polyfit(x, y, 1);
    yn = polyval(fitLineFcn, locs);
    plot(locs, yn, 'r*');
    k = 'm';
    if (fitLineFcn(1) < -200)
        k = 'g';
        
    end
    plot(x, y, k);
    s = s+fitLineFcn(1);
end
suma = s/length(d)
xlabel({"Średni wsp. kierunkowy = " num2str(suma)})
title("Dane rzeczywiste", filename );
