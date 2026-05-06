function [modes, V, Eleft, Eright] = test_modes_vec(By, y, x, discont)
if size(y, 2) ~= size(discont, 2) + 1
    error('Number of functions should equal number of segments.')
end
try
    yend = y{end}(x(end));
    y1 = y{1}(x(1));
catch
    y1 = y(x(1));
    yend = y(x(end));
end
tol = 10^(-8);
[vin, ein] = eig(By(y1));
[vout, eout] = eig(By(yend));
ein = diag(ein);
eout = diag(eout);
d = size(vin, 2);
decayleft = zeros([1, d]);
waveleft = zeros([1, d]);
% waveleft_out = zeros([1, d]);
decayright = zeros([1, d]);
waveright = zeros([1, d]);
%waveright_out = zeros([1, d]);
f1 = 0;
f2 = 0;
f3 = 0;
f4 = 0;
for n = 1:d
    if real(ein(n)) > tol
        f1 = f1 + 1;
        decayleft(f1) = n;
    elseif abs(real(ein(n))) < tol && abs(imag(ein(n))) > tol
        f2 = f2 + 1;
        waveleft(f2) = n;
    end
    if real(eout(n)) < -tol
        f3 = f3 + 1;
        decayright(f3) = n;
    elseif abs(real(eout(n))) < tol && abs(imag(eout(n))) > tol
        f4 = f4 + 1;
        waveright(f4) = n;
    end
end

decayleft = decayleft(1:f1);
waveleft = waveleft(1:f2);
decayright = decayright(1:f3);
waveright = waveright(1:f4);
if size(discont, 2) == 0
    h = round(size(y, 2)/2);
    vin = vin(:, [decayleft, waveleft]);
    vout = vout(:, [decayright, waveright]);
    [vleft, Vleft] = runge_vec(vin, By, yleft, x(1:h));
    [vright, Vright] = runge_vec(vout, By, yright, flip(x(h:end)));
    % modes = 2*ones([1, 4]);
    [Qleft, Rleft] = qr(vleft, "econ");
    [Qright, Rright] = qr(vright, "econ");
    [p, e] = eigs((Qleft*Qleft')*(Qright*Qright'), 1, 1);
    modes = 1-e;
    vt_left = Rleft\(Qleft'*p);
    vt_right = Rright\(Qright'*p);
    vec_left = zeros(size(Vleft, 1), size(Vleft, 3));
    vec_right = zeros(size(Vright, 1), size(Vright, 3));
    for j = 1:size(Vleft, 3)
        vec_left(:, size(Vleft, 3)+1-j) = Vleft(:, :, size(Vleft, 3)+1-j)*vt_left;
        %vec_left(:, size(Vleft, 3)+1-j) = vec_left(:, size(Vleft, 3)+1-j)/norm(vec_left(:, size(Vleft, 3)+1-j));
        vec_right(:, j) = Vright(:, :, size(Vright, 3) + 1-j)*vt_right;
        %vec_right(:, j) = vec_right(:, j)/norm(vec_right(:, j));
    end
else
    vleft = vin(:, decayleft);
    vright = vout(:, decayright);
    if mod(size(discont, 2), 2) == 1
        discont = [1, discont, size(x, 2)];
    else
        vleft = runge(vleft, By, y{1}, x(1:discont(1)));
        discont = [discont, size(x, 2)];
    end


Eleft = ein([decayleft, waveleft]);
Eright = eout([decayright, waveright]);
V = [vec_left, p, vec_right];
