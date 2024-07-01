%Code subject: Simple demo script for drawing a smoothed Bezier curve
%Programmer: Aaron Wetzler, aaronwetzler@gmail.com
%Date:12/12/2009

k=getPoints(5);
[cp1,cp2]=findControlPoints(k);
B = getBezier(k,cp1,cp2);
scatter(k(:,1),k(:,2)); hold on
plot(B(:,1),B(:,2));