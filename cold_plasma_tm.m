function B = cold_plasma_tm(E, kx)
op = 2;
s1 = [[0, 1];[1, 0]];
B = @(Om) 1j*s1*((1/(E*(-op^2-Om^2+E^2)))*[[op^2*(E^2-op^2), 1j*op^2*Om*kx];...
    [-1j*op^2*Om*kx, kx^2*(E^2-Om^2)]] - E*eye(2));