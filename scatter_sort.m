function [ksc, Esc] = scatter_sort(Aedge, k, E, tol)
% Inputs: 
% --Aedge; output of spectral calculation by grid_adapt.m or
% residue_map.m. 
% --k; vector, size of 2nd dimension of Aedge, values
% corresponding to the k values of each column of Aedge. 
% --E; vector, size of 1st dimension of Aedge, values corresponding to E
% values of each the rows of Aedge.
% Outputs: ksc, Esc, pairs of k, E values to plot as a scatter plot
% (scatter(k, E)).
Esc = [];
ksc = [];
[Ne, Nk] = size(Aedge);
for i = 1:Nk
    ky = k(i);
    for j = 1:Ne
        En = E(j);
        if Aedge(j, i) < tol
            ksc(end+1) = ky;
            Esc(end + 1) = En;
        end
    end
end