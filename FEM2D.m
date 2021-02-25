clear;
close all;

numberOfPDE = 1;
model = createpde;
gm = geometryFromEdges(model,@circleg);
model.Geometry = gm;
pdegplot(model,'CellLabels','on','FaceAlpha',0.5);
title 'Geometry With Edge Labels Displayed';
xlabel x
ylabel y

% Create and view a finite element mesh for the problem.
mes = generateMesh(model,'GeometricOrder','linear');
figure
pdemesh(model);
xlabel x
ylabel y

connArr = mes.Elements;
nodeIndex = mes.Nodes;

epsilon = 8.85E-12;
mu = 4*pi*1E-7;
alpha = 1/mu;  %new1
%f = 2.5E9;
f = 50;
c = 3E8;
lambda = c/f;
ksq = (2*pi/lambda)^2;
theta = 60*pi/180;
beta = ksq*(epsilon - (1/mu)*((sin(theta))^2));
fe = 10;
Kfinal = zeros(length(nodeIndex));
Bfinal = zeros(1,length(nodeIndex));

nodeBS = [];
for i = 1:length(nodeIndex)
    if nodeIndex(1,i) < -0.8 || nodeIndex(1,i) > 0.8
        if nodeIndex(2,i) < -0.8 || nodeIndex(2,i) > 0.8
            nodeBS = [nodeBS i];
        end
    end
end

for e = 1:length(connArr)
    node1 = connArr(1,e);
    node2 = connArr(2,e);
    node3 = connArr(3,e);
    x1 = nodeIndex(1, node1);  y1 = nodeIndex(2, node1);
    x2 = nodeIndex(1, node2);  y2 = nodeIndex(2, node2);
    x3 = nodeIndex(1, node3);  y3 = nodeIndex(2, node3);
    a1 = x2*y3 - y2*x3;
    a2 = x3*y1 - y3*x1;
    a3 = x1*y2 - y1*x2;
    b1 = y2 - y3;
    b2 = y3 - y1;
    b3 = y1 - y2;
    b = [b1, b2, b3];
    c1 = x3 - x2;
    c2 = x1 - x3;
    c3 = x2 - x1;
    c = [c1, c2, c3];
    delta = 0.5*(b1*c2 - b2*c1);
%     N1 = @(x,y) (0.5/delta)*(a1 + b1*x + c1*y);
%     N2 = @(x,y) (0.5/delta)*(a2 + b2*x + c2*y);
%     N3 = @(x,y) (0.5/delta)*(a3 + b3*x + c3*y);
%     N(i,:) = [N1; N2; N3];
    for i = 1:3
        nie = connArr(i,e);
        for j = 1:3
            if i == j
                del = 1;
            else
                del = 0;
            end
            K(i,j) = (0.25*alpha/delta)*((b(i)*b(j)) + (c(i)*c(j))) + ...
                delta*beta*(1 + del)/12;
            nje = connArr(j,e);
            Kfinal(nie,nje) = Kfinal(nie,nje) + K(i,j);
        end
        b(i) = delta*fe/3;
        Bfinal(nie) = Bfinal(nie) + b(i);
    end
    eleK(e) = {K};
    eleB(e) = {b};
end    

p = 0;

for n = 1:length(nodeBS)
    Kfinal(nodeBS(n),nodeBS(n)) = 1;
    for i = 1:length(Kfinal(nodeBS(n),:))
        if i ~= nodeBS(n)
            Kfinal(nodeBS(n),i) = 0;
        end
    end
    Bfinal(nodeBS(n)) = p;
    for i = 1:length(Bfinal)
        if i ~= nodeBS(n)
            Bfinal(i) = Bfinal(i) - p*Kfinal(i,nodeBS(n));
            Kfinal(i,nodeBS(n)) = 0;
        end
    end
end
        
phi = linsolve(Kfinal,Bfinal');

figure
pdeplot(mes,'XYData',phi)
title('u(1)')




