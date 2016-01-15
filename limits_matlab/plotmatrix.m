function plotmatrix(A)

%%% Plot matrix A as image with color scale
%%% and write the values in each cell.
%%% Use a red-blue color scale (code from Adam Auton)
%%% Example: plotmatrix(2*rand(10,15)-1)
%%%
%%% Author: Didier Gonze
%%% Created: 30/7/2014
%%% Updated: 6/8/2015 - Sophie de Buyl 


N=length(A(:,1));

Amin=min(min(A));
Amax=max(max(A));
Alim=max(abs(Amin),abs(Amax));

clim=[-Alim Alim];   % scale according to the data
%clim=[-1 1];       % if normalized

imagesc(A',clim)

% colorbar      % uncomment to display color scale

colormap(redblue)

% set(gca,'YDir','normal')   % reverse Y axis


%%% Text

if N<20 % to remove explicit values of matrix enteries when the matrix is too big
    
ndec=2;   % number of decimals
fsize=10; % fontsize

for i=1:size(A,1)
    for j=1:size(A,2)
        a=round(A(i,j)*(10^ndec))/(10^ndec);
        text(i,j,sprintf('%g',a),'fontsize',fsize,'HorizontalAlignment','center')
    end
end

else
    dddd=2;

end

% ============================================================
% Color map: redblue
% ============================================================

function c = redblue(m)
%REDBLUE    Shades of red and blue color map
%   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright blue, range through shades of
%   blue to white, and then through shades of red to bright red.
%   REDBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redblue)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.

%   Adam Auton, 9th October 2009

if nargin < 1, m = size(get(gcf,'colormap'),1); end

if (mod(m,2) == 0)
    % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
    m1 = m*0.5;
    r = (0:m1-1)'/max(m1-1,1);
    g = r;
    r = [r; ones(m1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
else
    % From [0 0 1] to [1 1 1] to [1 0 0];
    m1 = floor(m*0.5);
    r = (0:m1-1)'/max(m1,1);
    g = r;
    r = [r; ones(m1+1,1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end

c = [r g b];

c=flipud(c);

