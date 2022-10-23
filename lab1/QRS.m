q = []; r = []; s = [];
% q(length(QRSdetected)) = 0; r(length(QRSdetected)) = 0; s(length(QRSdetected)) = 0;
for i = QRSdetected
    Q = i;
    next = Q+1;
    while(x(Q) > x(next))
        Q = next;
        next = Q+1;
    end
    q = [q; Q];

    R = Q+1;
    next = R+1;
    while(x(R) <= x(next))
        R = next;
        next = R+1;
    end
    r = [r; R];
    
    S = R+1;
    next = S+1;
    while(x(S) > x(next))
        S = next;
        next = S+1;
    end
    s = [s; S];
    
end

figure, hold on;
plot(t, x);
plot(t(q),x(q),'r^');
plot(t(r),x(r),'g^');
plot(t(s),x(s),'b^');
hold off;