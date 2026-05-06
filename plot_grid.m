function [] = plot_grid(Aedge, AL, AR, ALR, k, E, tol)
edge = [];
bulkR = [];
bulkL = [];
bulkLR = [];
xedge = [];
xL = [];
xR = [];
xLR = [];
[Ne, Nk] = size(AL);
for i = 1:Nk
    ky = k(i);
    for j = 1:Ne
        En = E(j);
        if Aedge(j, i) < tol
            xedge(end+1) = ky;
            edge(end + 1) = En;
        end
        if AL(j, i) < tol
            xL(end+1) = ky;
            bulkL(end+1) = En;
        end
        if AR(j, i) < tol
            xR(end+1) = ky;
            bulkR(end+1) = En;
        end
        if ALR(j, i) < tol
            xLR(end+1) = ky;
            bulkLR(end+1) = En;
        end
    end
end
figure();
hold on
scatter(xedge, edge, 6, 'r', '.');
scatter(xL, bulkL, 6, 'b', '.');
scatter(xR, bulkR, 6, 'g', '.');
scatter(xLR, bulkLR, 6, 'k', '.');
legend('Edge Mode', 'Bulk Left', 'Bulk Right', 'Bulk Both');

