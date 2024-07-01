%Code subject: Simple display of a set of 2d cubic Bezier curve segments
%Programmer: Aaron Wetzler, aaronwetzler@gmail.com
%Date:12/12/2009
%
%displayBezier(knots,cp1,cp2)
%knots - Set of points with each row holding a point . In this case they must be 2-d. The code can
%easily be changed to display 3d points as well.
%[cp1,cp2]- the set of control points used to determine the characteristics
%of the Bezier segments
%
%Each segment will be drawn with 11 subsegments i.e. straight lines

function [B] = getBezier(knots,cp1,cp2)

n=size(cp1,1);
dim=size(cp1,2);

t=[0:0.1:1]';
t=repmat(t,1,dim);
lt=size(t,1);

cpnts=cat(3, knots(1:end-1,:), cp1, cp2, knots(2:end,:));
hold on;
B = zeros(1,2);
for i=1:n
    B = [B; repmat(cpnts(i,:,1),lt,1).*((1-t).^3) + 3*repmat(cpnts(i,:,2),lt,1).*(t.*(1-t).^2) + 3*repmat(cpnts(i,:,3),lt,1).*((t.^2).*(1-t))+repmat(cpnts(i,:,4),lt,1).*(t.^3)];
    %plot(B(:,1),B(:,2),'r');
    %plot(knots(:,1),knots(:,2),'.b');
end
B = B(2:end,:);
