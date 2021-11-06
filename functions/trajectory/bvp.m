function s1 = bvp(dt, runtime, dest, s0)
    p0 = s0(1:3);
    v0 = s0(4:6);
    a0 = s0(7:9);
    pf = dest(1:3);
    vf = dest(4:6);
    af = dest(7:9);
    timeint = linspace(0,runtime,runtime/dt);

    % s0 = [p0, v0, a0];
    % sf = [pf, vf, af];
    T = runtime;

    delta =  [(pf - p0 -v0 * T - 0.5 * a0 * T^2) ; ...
              (vf - v0 -a0 * T) ; ...
              (af - a0)];
    m = [720, -360*T, 60*T^2 ; -360*T, 168*T^2, -24*T^3 ; 60*T^2, -24*T^3, 3*T^4];
    M = zeros(length(m)*width(m),length(m)*width(m));
    for i = 1:length(m)*width(m)
        M1 = eye(3);
        l = mod(i,length(m)) + length(m)*(~mod(i,length(m)));
        h = ceil(i/3);
        M1(M1==1) = m(l,h);
        %  fprintf('(%d) %d, %d\n',i,1+(l-1)*length(m),1+(h-1)*length(m));
        M(1+(l-1)*length(m):1+(l-1)*length(m)+2,1+(h-1)*length(m):1+(h-1)*length(m)+2) = M1;
    end

    abg = 1/T^5 * M * delta;
    alpha = abg(1:3); beta = abg(4:6); gamma = abg(7:9);

    for ts = 1:width(timeint)
        t = timeint(ts);
        sOpt(:,ts) = [(alpha/120 * t^5 + beta/24 * t^4 + gamma/6 * t^3 + a0/2 * t^2 + v0 * t + p0); ...
                         (alpha/24 * t^4 + beta/6 * t^3 + gamma/2 * t^2 + a0 * t + v0); ...
                         (alpha/6 * t^3 + beta/2 * t^2 + gamma * t + a0)];
    end
    s1 = sOpt;
end

