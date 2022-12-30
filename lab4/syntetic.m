
x = d1(1,:);
y = d1(2,:);
plot(x, y,'o')
% zadanyP = .6;
param = 0.15;
procent = 0;% startowy
while( procent < zadanyP )
    [a] = [randi(length(x)) randi(length(x))];
    fitLineFcn = polyfit(x(a), y(a), 1); % fit function using polyfit
    yn = polyval(fitLineFcn, x);
    m = mean(y);
    sigma = param*m;
    
    w = y < yn+sigma & y > yn-sigma;
   
    procent = sum(w)/length(x)
end

    hold on 
    plot(x, yn+sigma, "r")
    plot(x, yn, "g")
    plot(x, yn-sigma, "b")
     
    plot(x(w), y(w), 'om')
    title("Dane syntetyczne ", filename)
    xlabel(sprintf("Próg zadany %g%%, obliczony = %g%% sigma = %g", zadanyP, procent, sigma ))
