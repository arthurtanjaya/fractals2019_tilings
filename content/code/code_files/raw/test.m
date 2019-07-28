% Script for testing
% Currently very rudimentary

close all
clear
clc

x = Vertex([0 1 2])
y = Vertex([1 0 1])
z = Vertex([1 0 0])

G = GraphApprox(3)

figure
hold on
%x.draw(1)
%y.draw(1)
%z.draw(1)
G.draw(1)
hold off

u = [1  0 0;
     2  3 0;
     4  5 6];

v = [7  0 0;
     8 -9 0;
     3  2 1];

GVertexIP(u, v)
