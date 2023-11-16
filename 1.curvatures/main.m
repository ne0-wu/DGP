clear
close all

[V, T] = readObj('../mesh_samples/cathead.obj');

visualizeCurvature(V, T, gaussianCurvature(V, T));
visualizeCurvature(V, T, meanCurvature(V, T));