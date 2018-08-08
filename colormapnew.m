
% color map for sequential data (orange/brown) 

C =[102    37     6
    153    52     4
    204    76     2
    236   112    20
    251   154    41
    254   197    79
    254   227   145
    255   247   147
    255   255   229]/255;


x = (1:10:length(C)*10);
n = x(end);
% convert to HSV for interpolation
%C_HSV = rgb2hsv(C);
% interpolate hue value
for i = 1:3
C_interp(i,:) = interp1(x, C(:, i)', 1:n)';
end

C = C_interp';
% set colormap
colormap(C)
colorbar