clear; close all; clc;
% testing surface plots in tikz
% this script generates a data set for tikz to turn into a 3d surface plot

% define the two surfaces, f and g
f = @(x,y) 100*sin(x)+y.^3;
g = @(x,y) 100*sin(y)+x.^3;

% defining the "mesh grid"
x = linspace(-2*pi,2*pi,1e1);
y = linspace(-2*pi,2*pi,1e1);

[a,b] = meshgrid(x,y);
x2 = zeros(numel(a),1);
y2 = zeros(numel(a),1);

for i = 1:numel(a)
    x2(i) = a(i);
    y2(i) = b(i);
end

% calculating the surfaces (z and w) from the grid
z = f(x2,y2);
w = g(x2,y2);

% Get the matrixes ready for export to tikz
A = [x2 y2 z];
A2 = A';

B = [x2 y2 w];
B2 = B';

% normally I would use save, and save the data to a .dat, like so...

%save('testdata.dat','A','-ascii')

% but, the formatting for the 3d plot for tikz requires specific formatting
% crazy crazy formatting here. Tikz wants a line break to be clear. But if
% they all have the same line breaks, it doesn't recognize the line breaks
% So only one line gets a line break. Idk why, but it works.

% to underrstand a little more about the "empty line" thing, see the tikz
% documentation http://pgfplots.sourceforge.net/pgfplots.pdf, page 124...

% otherwise, here is my workaround for exporting the data.
% exporting the first surface...
fileID = fopen('testdata.dat','w');
formatSpec1 = '%f %f %f \n';
formatSpec2 = '%f %f %f \n\n';
fprintf(fileID,formatSpec1,A2(1:3));
fprintf(fileID,formatSpec2,A2(4:end));
fclose('all');

% exporting the second surface...
fileID = fopen('testdata2.dat','w');
formatSpec1 = '%f %f %f \n';
formatSpec2 = '%f %f %f \n\n';
fprintf(fileID,formatSpec1,B2(1:3));
fprintf(fileID,formatSpec2,B2(4:end));
fclose('all');