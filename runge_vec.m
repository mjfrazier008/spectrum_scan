function [vout, V] = runge_vec(vin, By, y, x)
n = size(x, 2);
[k, m] = size(vin);
vout = vin;
V = zeros(k, m, n);
for j = 1:n-1
    for k = 1:m
        V(:, k, j) = vout(:, k);
    end
    dx = x(j+1)-x(j);
    B0 = By(y(x(j)));
    B1 = By(y(x(j) + dx/2));
    B2 = By(y(x(j+1)));
    k1 = B0*vout;
    k2 = B1*(vout + (dx/2)*k1);
    k3 = B1*(vout + (dx/2)*k2);
    k4 = B2*(vout + dx*k3);
    vout = vout + (dx/6)*(k1 + 2*k2 + 2*k3 + k4);
end
for k = 1:m
        V(:, k, n) = vout(:, k);
end
%for j = 1:m
%    vout(:, j) = vout(:, j)/norm(vout(:, j));
%end