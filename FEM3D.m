clear;
close all;

numberOfPDE = 1;
model = createpde;
dLen = 1;
%rad = 0.06;
rad = 0.15;
gm = multicylinder([rad rad+0.01],dLen,'Void',[true,false]);
model.Geometry = gm;
pdegplot(model,'CellLabels','on','FaceAlpha',0.5);
title 'Geometry With Edge Labels Displayed';
xlabel x
ylabel y

% Create and view a finite element mesh for the problem.
mes = generateMesh(model,'GeometricOrder','linear');
% mes = generateMesh(model);
figure
pdemesh(model);
xlabel x
ylabel y

connArr = mes.Elements;
nodeIndex = mes.Nodes;

epsilon = 8.85E-12;
mu = 4*pi*1E-7;
alphaE = 1/mu;  %new1
alphaM = 1/epsilon;  %new1
f = 2.5E9;
% f = 60E9;
c = 3E8;
lambda = c/f;
k = (2*pi/lambda)^2;
betaE = k*k*epsilon;
betaM = k*k*mu;
%fe = 10;
feE = -1i*k*50*1;
feM = -1i*k*50*1;
KfinalE = zeros(length(nodeIndex));
BfinalE = zeros(1,length(nodeIndex));
KfinalM = zeros(length(nodeIndex));
BfinalM = zeros(1,length(nodeIndex));

nodeBS = [];
for i = 1:length(nodeIndex)
    r = sqrt(((nodeIndex(1,i))^2) + ((nodeIndex(2,i))^2));
    if r < -rad || r > rad
        nodeBS = [nodeBS i];
    end
end

for e = 1:length(connArr)
    node1 = connArr(1,e);
    node2 = connArr(2,e);
    node3 = connArr(3,e);
    node4 = connArr(4,e);
    x1 = nodeIndex(1, node1);  y1 = nodeIndex(2, node1);  z1 = nodeIndex(3, node1);
    x2 = nodeIndex(1, node2);  y2 = nodeIndex(2, node2);  z2 = nodeIndex(3, node2);
    x3 = nodeIndex(1, node3);  y3 = nodeIndex(2, node3);  z3 = nodeIndex(3, node3);
    x4 = nodeIndex(1, node4);  y4 = nodeIndex(2, node4);  z4 = nodeIndex(3, node4);
    a1 = x2*(y3*z4 - y4*z3) + x3*(y4*z2 - y2*z4) + x4*(y2*z3 - y3*z2);
    a2 = x3*(y4*z1 - y1*z4) + x4*(y1*z3 - y3*z1) + x1*(y3*z4 - y4*z3);
    a3 = x4*(y1*z2 - y2*z1) + x1*(y2*z4 - y4*z2) + x2*(y4*z1 - y1*z4);
    a4 = x1*(y2*z3 - y3*z2) + x2*(y3*z1 - y1*z3) + x3*(y1*z2 - y2*z1);
    b1 = y2*z3 - z2*y3;  b2 = y3*z4 - z3*y4;  b3 = y4*z1 - z4*y1;
    b4 = y1*z2 - z1*y2;  c1 = z3 - z4;  c2 = z4 - z1;  c3 = z1 - z2;
    c4 = z2 - z3; d1 = x4 - x3;  d2 = x1 - x4;  d3 = x2 - x1;  d4 = x3 - x2;
    bE = [b1, b2, b3, b4];  c = [c1, c2, c3, c4];  d = [d1, d2, d3, d4];
    V = (1/6)*det([1 1 1 1; x1 x2 x3 x4; y1 y2 y3 y4; z1 z2 z3 z4]);
    for i = 1:4
        nie = connArr(i,e);
        for j = 1:4
            if i == j
                del = 1;
            else
                del = 0;
            end
            KE(i,j) = (alphaE/(36*V))*((bE(i)*bE(j)) + (c(i)*c(j)) + (d(i)*d(j))) + ...
                V*betaE*(1 + del)/20;
            KM(i,j) = (alphaM/(36*V))*((bE(i)*bE(j)) + (c(i)*c(j)) + (d(i)*d(j))) + ...
                V*betaM*(1 + del)/20;
            nje = connArr(j,e);
            KfinalE(nie,nje) = KfinalE(nie,nje) + KE(i,j);
            KfinalM(nie,nje) = KfinalM(nie,nje) + KM(i,j);
        end
        bE(i) = V*feE/4;
        BfinalE(nie) = BfinalE(nie) + bE(i);
        bM(i) = V*feM/4;
        BfinalM(nie) = BfinalM(nie) + bM(i);
    end
    eleKE(e) = {KE};
    eleBE(e) = {bE};
    eleKM(e) = {KM};
    eleBM(e) = {bM};
end    

pE = 1;
pM = 0;

for n = 1:length(nodeBS)
    KfinalE(nodeBS(n),nodeBS(n)) = 1;
    KfinalM(nodeBS(n),nodeBS(n)) = 1;
    for i = 1:length(KfinalE(nodeBS(n),:))
        if i ~= nodeBS(n)
            KfinalE(nodeBS(n),i) = 0;
        end
    end
    for i = 1:length(KfinalM(nodeBS(n),:))
        if i ~= nodeBS(n)
            KfinalM(nodeBS(n),i) = 0;
        end
    end
    BfinalE(nodeBS(n)) = pE;
    BfinalM(nodeBS(n)) = pM;
    for i = 1:length(BfinalE)
        if i ~= nodeBS(n)
            BfinalE(i) = BfinalE(i) - pE*KfinalE(i,nodeBS(n));
            KfinalE(i,nodeBS(n)) = 0;
        end
    end
    for i = 1:length(BfinalM)
        if i ~= nodeBS(n)
            BfinalM(i) = BfinalM(i) - pM*KfinalM(i,nodeBS(n));
            KfinalM(i,nodeBS(n)) = 0;
        end
    end
end
        
phiE = linsolve(KfinalE,BfinalE');

realPhiE = real(phiE);
imgPhiE = imag(phiE);

phiM = linsolve(KfinalM,BfinalM');

phiZFE = [];
phiZME = [];
phiZFM = [];
phiZMM = [];

for i = 1:length(nodeIndex)
    if nodeIndex(3,i) == dLen
        phiZFE = [phiZFE abs(phiE(i))];
    end
    if nodeIndex(3,i) > ((dLen/2) - 0.1) && nodeIndex(3,i) < ((dLen/2) + 0.1)
        phiZME = [phiZME abs(phiE(i))];
    end
    if nodeIndex(3,i) == dLen
        phiZFM = [phiZFM abs(phiM(i))];
    end
    if nodeIndex(3,i) > ((dLen/2) - 0.1) && nodeIndex(3,i) < ((dLen/2) + 0.1)
        phiZMM = [phiZMM abs(phiM(i))];
    end
end

figure
pdeplot3D(mes,'ColorMapData',real(phiE))
title('E')

figure
pdeplot3D(mes,'ColorMapData',real(phiM))
title('M')






