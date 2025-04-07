clc;
clear all;

results = load("results.mat")

figure;
quiver(flipud(results.u),flipud(results.v),10)
title("Velocity Plot")
axis equal tight