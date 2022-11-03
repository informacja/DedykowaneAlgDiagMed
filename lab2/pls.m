function output = pls(f, t_k, x)%Periodogram Lomb-Scargle
    tau = 1/(4*pi*f) * atan(sum(sin(4*pi*f*t_k))/sum(cos(4*pi*f*t_k)));

    part1 = sum(x .* cos(2*pi*f*(t_k - tau))).^2;
    part2 = sum(cos(2*pi*f*(t_k - tau)).^2);

    part3 = sum(x .* sin(2*pi*f*(t_k - tau))).^2;
    part4 = sum(sin(2*pi*f*(t_k - tau)).^2);

    output = 0.5 * ( part1/part2 + part3/part4 );

end

