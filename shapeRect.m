t0 = (pi/24:pi/12:2*pi)';
% x0 = 4*cos(t0);
% y0 = 4*sin(t0);
x0 = [-3,3].*cos(t0);
x0 = x0(:);
y0 = [-3,3].*sin(t0);
y0 = y0(:);
circShp0 = alphaShape(x0,y0);
% figure
% plot(circShp0)

%[xg, yg] = meshgrid(-3:0.25:3);
[xg, yg] = meshgrid(circShp0.Points(:,1),circShp0.Points(:,2));
xg = xg(:);
yg = yg(:);
t = (pi/24:pi/24:2*pi)';
x = cos(t);
y = sin(t);
circShp = alphaShape(x,y,2);

in = inShape(circShp,xg,yg);
xg = [xg(~in); x];
yg = [yg(~in); y];
zg = ones(numel(xg),1);
xg = repmat(xg,5,1);
yg = repmat(yg,5,1);
zg = zg*(0:.25:1);
zg = zg(:);

%[xg1, yg1] = meshgrid(-3:0.25:3);
[xg1, yg1] = meshgrid(circShp0.Points(:,1),circShp0.Points(:,2));
xg1 = xg1(:);
yg1 = yg1(:);
zg1 = ones(numel(xg1),1);
xg1 = repmat(xg1,5,1);
yg1 = repmat(yg1,5,1);
zg1 = zg1*(-1.25:.25:-0.25);
zg1 = zg1(:);

%[xg2, yg2] = meshgrid(-3:0.25:3);
[xg2, yg2] = meshgrid(circShp0.Points(:,1),circShp0.Points(:,2));
xg2 = xg2(:);
yg2 = yg2(:);
zg2 = ones(numel(xg2),1);
xg2 = repmat(xg2,5,1);
yg2 = repmat(yg2,5,1);
zg2 = zg2*(1.25:0.25:2.25);
zg2 = zg2(:);

xg3 = [xg1; xg; xg2];
yg3 = [yg1; yg; yg2];
zg3 = [zg1; zg; zg2];

shp = alphaShape(xg3,yg3,zg3);
figure
plot(shp)

shp1 = alphaShape([xg1;xg],[yg1; yg],[zg1; zg]);
figure
plot(shp1)

%Obtain a surface mesh of the alphaShape object.
[elements,nodes] = boundaryFacets(shp);

% Put the data in the correct shape for geometryFromMesh.
nodes = nodes';
elements = elements';

%Create a PDE model and import the surface mesh.
model = createpde();
geometryFromMesh(model,nodes,elements);

%View the geometry and face numbers.
figure
pdegplot(model,'FaceLabels','on','FaceAlpha',0.5)

generateMesh(model);