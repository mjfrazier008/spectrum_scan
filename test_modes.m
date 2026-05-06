function modes = test_modes(By, y, x, discont)
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
%waveleft = zeros([1, d]);
decayright = zeros([1, d]);
%waveright = zeros([1, d]);
f1 = 0;
%f2 = 0;
f3 = 0;
%f4 = 0;
for n = 1:d
    if real(ein(n)) > tol
        f1 = f1 + 1;
        decayleft(f1) = n;
    %elseif abs(real(ein(n))) < tol && abs(imag(ein(n))) > tol
        %f2 = f2 + 1;
        %waveleft(f2) = n;
    end
    if real(eout(n)) < -tol
        f3 = f3 + 1;
        decayright(f3) = n;
    %elseif abs(real(eout(n))) < tol && abs(imag(eout(n))) > tol
        %f4 = f4 + 1;
        %waveright(f4) = n;
    end
end
decayleft = decayleft(1:f1);
%waveleft = waveleft(1:f2);
decayright = decayright(1:f3);
%waveright = waveright(1:f4);
if size(discont, 2) == 0
    h = round(size(x, 2)/2);
    vleft = runge(vin(:, decayleft), By, y, x(1:h));
    vright = runge(vout(:, decayright), By, y, flip(x(h:end)));
    % vleft = runge(vin(:, decayleft), By, y, x);
    % vright = vout(:, decayright);
    %modes = 2*ones([1, 4]);
    try
        Qdecayr = orth(vright);
        %Qwaver = orth(vright(:, waveright));
        Qdecayl = orth(vleft);
        %Qwavel = orth(vleft(:, waveleft));
        modes = abs(eigs((Qdecayr*Qdecayr')*(Qdecayl*Qdecayl'), 1, 1) -1);
    catch
        modes = 1;
    end
else
    vleft = vin(:, decayleft);
    vright = vout(:, decayright);
    if mod(size(discont, 2), 2) == 1
        discont = [1, discont, size(x, 2)];
        for j = 1:(size(discont, 2)-1)/2
            vleft = runge(vleft, By, y{j}, x(discont(j):discont(j+1)));
            vright = runge(vright, By, y{end-j+1}, flip(x(discont(end-j):discont(end-j+1))));
        end
    else
        % vright = runge(vright, By, y{end}, flip(x(discont(end):end)));
        % discont = [1, discont];
        vleft = runge(vleft, By, y{1}, x(1:discont(1)));
        discont = [discont, size(x, 2)];
        for j = 1:(size(discont, 2)-1)/2
            vleft = runge(vleft, By, y{j+1}, x(discont(j):discont(j+1)));
            vright = runge(vright, By, y{end-j+1}, flip(x(discont(end-j):discont(end-j+1))));
        end
    end
    % discont = [1, discont, size(x, 2)];
    % for j = 1:size(discont, 2)-1
    %     vleft = runge(vleft, By, y{j}, x(discont(j):discont(j+1)));
    % end
    try 
        Qdecayr = orth(vright);
        Qdecayl = orth(vleft);
        modes = abs(eigs((Qdecayr*Qdecayr')*(Qdecayl*Qdecayl'), 1, 1) -1);
    catch
        modes = 1;
    end
end


%modes(2) = abs(eigs((Qwaver*Qwaver')*(Qdecayl*Qdecayl'), 1, 1) -1);
%modes(3) = abs(eigs((Qdecayr*Qdecayr')*(Qwavel*Qwavel'), 1, 1) -1);
%modes(4) = abs(eigs((Qwaver*Qwaver')*(Qwavel*Qwavel'), 1, 1) -1);
% v1 = Qdecay*Qdecay'*vleft;
% v2 = Qwave*Qwave'*vleft;
% if f1 > 0
%     for i = 1:f1
%         if f3 > 0
%             modes(1) = norm(vleft(:, decayleft(i)) - v1(:, decayleft(i)));
%         end
%         if f4 > 0
%             modes(2) = norm(vleft(:, decayleft(i)) - v2(:, decayleft(i)));
%         end
%     end
% end
% if f2 > 0
%     for i = 1:f2
%         if f3 > 0
%             modes(3) = norm(vleft(:, waveleft(i)) - v1(:, waveleft(i)));
%         end
%         if f4 > 0
%             modes(4) = norm(vleft(:, waveleft(i)) - v2(:, waveleft(i)));         
%         end
%     end
% end
