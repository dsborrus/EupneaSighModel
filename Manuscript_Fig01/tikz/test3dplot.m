clear; close all; clc;
f = @(x,y) 100*sin(x)+y.^3;
g = @(x,y) 100*sin(y)+x.^3;

x = linspace(-2*pi,2*pi,1e1);
y = linspace(-2*pi,2*pi,1e1);

[a,b] = meshgrid(x,y);
x2 = zeros(numel(a),1);
y2 = zeros(numel(a),1);

for i = 1:numel(a)
    x2(i) = a(i);
    y2(i) = b(i);
end

z = f(x2,y2);
w = g(x2,y2);

A = [x2 y2 z];
A2 = A';

B = [x2 y2 w];
B2 = B';

%save('testdata.dat','A','-ascii')

% crazy crazy formatting here. Tikz wants a line break to be clear. But if
% they all have the same line breaks, it doesn't recognize the line breaks
% So only one line gets weird line breaks. Idk it works.
fileID = fopen('testdata.dat','w');
formatSpec1 = '%f %f %f \n';
formatSpec2 = '%f %f %f \n\n';
fprintf(fileID,formatSpec1,A2(1:3));
fprintf(fileID,formatSpec2,A2(4:end));
fclose('all');

fileID = fopen('testdata2.dat','w');
formatSpec1 = '%f %f %f \n';
formatSpec2 = '%f %f %f \n\n';
fprintf(fileID,formatSpec1,B2(1:3));
fprintf(fileID,formatSpec2,B2(4:end));
fclose('all');