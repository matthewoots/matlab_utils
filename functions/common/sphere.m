function [x, y, z] = sphere(r,c,n)
    theta = (-n:2:n)/n*pi;
    phi = (-n:2:n)'/n*pi/2;
    cosphi = cos(phi); cosphi(1) = 0; cosphi(n+1) = 0;
    sintheta = sin(theta); sintheta(1) = 0; sintheta(n+1) = 0;
    x = r(1)*cosphi*cos(theta) + c(1);
    y = r(2)*cosphi*sintheta + c(2);
    z = r(3)*sin(phi)*ones(1,n+1) + c(3);
end

