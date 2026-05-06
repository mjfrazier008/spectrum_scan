function [Aedge] = residue_map(kbounds, Ebounds,...
    Nk, Ne, B, m, x, discont)
E = linspace(Ebounds(1), Ebounds(2), Ne);
k = linspace(kbounds(1), kbounds(2), Nk);
Aedge = zeros(Ne, Nk);
% AL = zeros(size(Aedge));
% AR = zeros(size(Aedge));
% ALR = zeros(size(Aedge));

% CHANGE TO FOR LOOP IF NO PARALLEL CAPABILITY:
parfor i = 1:Nk
    ky = k(i);
    for j = 1:Ne
        En = E(j);
        H = B(En, ky);
        %modes = test_modes(H, m, x, discont);
        modes = test_modes_1(H, m, x, discont);
        Aedge(j, i) = modes;
        % AR(j, i) = modes(2);
        % AL(j, i) = modes(3);
        % ALR(j, i) = modes(4);
    end
end

