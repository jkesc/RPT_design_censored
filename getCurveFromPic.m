%This function is supposed to obtain reference length and origin from a
%picture, to create a set of points with similar ratios, and te "actual"
%sizes of the picture.
%Until now it doesn't really work, I think. but that may be the pictures
%fault?
function points = getCurveFromPic(figname)
pic=imread(figname);
figure(1)
hold off
imshow(pic);
hold on
input('Select two points with known length between them');
reference=ginput;
length=input('What was the length?');
rate=length./(sqrt(((reference(1,1)-reference(1,2)).^2+reference(2,1)-reference(2,2).^2)));
input('select an origin')
origin = ginput;
input('mark the points')
points=ginput;
plot(points(:,1),points(:,2));
points=points-origin;
points=points.*rate;
end