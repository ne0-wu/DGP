function fig = visualizeCurvature(V, T, curvature)

fig = figure;

trisurf(T, V(:, 3), V(:, 1), V(:, 2), curvature, 'EdgeColor', 'none')
axis equal
shading interp

colormap jet
colorbar