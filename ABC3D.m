lambda = 1;
rCyl = 0.1*lambda; % Radius of the circular PEC cylinder in wavelengths
rABC = 1.5*lambda; % Radius of the outer circular ABC boundary in wavelengths
h = 0.04*lambda; % Discretization size in wavelengths
len = 2*lambda; % Length of the waveguide

numberOfPDE = 1;
model = createpde;
% gm = multicylinder([rCyl rABC],len,'Void',[true,false]);
gm = multicylinder(rCyl,len);
model.Geometry = gm;
pdegplot(model,'CellLabels','on','FaceAlpha',0.5);
title 'Geometry With Edge Labels Displayed';
xlabel x
ylabel y
zlabel z

% Create and view a finite element mesh for the problem.
%mes = generateMesh(model,'GeometricOrder','linear','Hmax',0.3);
mes = generateMesh(model,'GeometricOrder','linear');
figure
pdemesh(model);
xlabel x
ylabel y

N = mes.Elements;
nElements = length(N);
nodeIndex = mes.Nodes;

% Lets find boundary edges
p = [nodeIndex(1,:) ; nodeIndex(2,:) ; nodeIndex(3,:)]; 
x = nodeIndex(1,:);
y = nodeIndex(2,:);
z = nodeIndex(3,:);
N = N';
boundedges = surftri(p',N);
% boundedges = boundedges';
hold on;
for i=1:length(boundedges)
    x1=x(boundedges(i,1)); x2=x(boundedges(i,2)); x3=x(boundedges(i,3));
    y1=y(boundedges(i,1)); y2=y(boundedges(i,2)); y3=y(boundedges(i,3));
    z1=z(boundedges(i,1)); z2=z(boundedges(i,2)); z3=z(boundedges(i,3));
    line([x1 x2 x3],[y1 y2 y3],[z1 z2 z3],'Color','r','LineWidth',2); % Plot the boundary edges
end
nDBC = 0; % no of Dirichlet Boundary Condition nodes
nABC = 0; % no of ABC nodes
nInt = 0; % no of internal nodes
nNodes=length(nodeIndex); % total no of nodes
DBC(1) = 0;
% Determine the ABC and DBC nodes
for jj=1:1:nNodes
    r0 = sqrt(x(jj)^2+y(jj)^2);
    if(abs(rABC - r0)<rABC)
%     if(abs(rABC - r0)<0.0001*rABC) % ABC will be imposed on a circle always
        nABC=nABC+1;
        ABC(nABC)=jj;
    else % either an internal node or a node on the object (scatterer)
        % if jj is an element of the be (bound edges) array
        % then it is a node on the object hence a DBC node
        % otherwise it is an internal node
        index=ArrSearch(boundedges,jj);
        if(index ~= -1) % DBC node
            nDBC=nDBC+1;
            DBC(nDBC)=jj;
%         else % internal node
%             nInt=nInt+1;
%             INT(nInt)=jj;
        end
    end
end
% flag each node
% flag = 1 for ABC node
% = 0 for internal node
% = ?-1 for DBC node
flag = zeros(nElements,4);
for e=1:nElements
    for i=1:4
        index=ArrSearch(ABC,N(e,i));
        if (index ~= -1)
            flag(e,i) = 1; % ABC node
        else
            index=ArrSearch(DBC,N(e,i));
            if(index ~= -1)
                flag(e,i) = -1; % DBC node
%             else
%                 flag(e,i) = 0; % internal node
            end
        end
    end
end
% hold on;
% for i=1:nDBC
%     nodeNum = DBC(i);
%     x1 = nodeIndex(1,nodeNum);
%     y1 = nodeIndex(2,nodeNum);
%     z1 = nodeIndex(3,nodeNum);
%     plot3(x1,y1,z1,'-o','Color','b','MarkerSize',10); % Plot the boundary edges
% end

% Create NS(s,i) array
nSegABC = 0;
for e=1:1:nElements
    if flag(e,1) == 1 && flag(e,2) == 1 && flag(e,3) == 1
        nSegABC = nSegABC + 1;
        NS(nSegABC, 1) = N(e,1);
        NS(nSegABC, 2) = N(e,2);
        NS(nSegABC, 3) = N(e,3);
    elseif flag(e,1) == 1 && flag(e,3) == 1 && flag(e,4) == 1
        nSegABC = nSegABC + 1;
        NS(nSegABC, 1) = N(e,1);
        NS(nSegABC, 2) = N(e,3);
        NS(nSegABC, 3) = N(e,4);
    elseif flag(e,1) == 1 && flag(e,4) == 4 && flag(e,2) == 1
        nSegABC = nSegABC + 1;
        NS(nSegABC, 1) = N(e,1);
        NS(nSegABC, 2) = N(e,4);
        NS(nSegABC, 3) = N(e,2);
    elseif flag(e,2) == 1 && flag(e,3) == 1 && flag(e,4) == 1
        nSegABC = nSegABC + 1;
        NS(nSegABC, 1) = N(e,2);
        NS(nSegABC, 2) = N(e,3);
        NS(nSegABC, 3) = N(e,4);
    end
end
fprintf(' RESULTS \n\n');
fprintf('The number of nodes = %6i\n', nNodes);
fprintf('The number of elements = %6i\n', nElements);
fprintf('The number of nodes on a Dirichlet boundary = %6i\n', nDBC);
fprintf('The number of nodes on an ABC boundary = %6i\n', nABC);

E_0 = 1; % Amplitude of incident electric field (V/m)
J = complex(0.,1.);
mur = 1; % Relative permeability (assuming vaccum)
epsr = 1; % Relative permittivity (asumming vacuum)
k_0 = 2*pi/lambda; % wavenumber
alpha = -1/mur;
kappa = 1/rABC;
gamma_1 = alpha*(J*k_0 + kappa/2 - J*kappa^2/(8*(J*kappa - k_0)));
mu_r = zeros(1,nElements); % relative permeability
eps_r = zeros(1,nElements); % relative permittivity
alpha_x = zeros(1,nElements); % coefficient of second order x term in PDE
alpha_y = zeros(1,nElements); % coefficient of second order y term in PDE
alpha_z = zeros(1,nElements);
beta = zeros(1,nElements); % coefficient of term linear in field in PDE
A = zeros(1,nElements); % area of element e
for e=1:1:nElements
    mu_r(e) = mur;
    eps_r(e) = epsr;
    alpha_x(e) = -1/mur;
    alpha_y(e) = -1/mur;
    alpha_z(e) = -1/mur;
    beta(e) = (k_0)*(k_0)*eps_r(e);
end
% Initialize global K matrix and the right hand side vector b
Kglobal = zeros(nNodes, nNodes);
btilde = zeros(nNodes, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FORMATION OF ELEMENT MATRICES AND ASSEMBLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coefficients of the linear interpolation functions
a = zeros(nElements,4);
b = zeros(nElements,4);
c = zeros(nElements,4);
be_tilde=zeros(1,nElements);
Klocal=zeros(nElements,4,4);
for e=1:1:nElements
    node1 = N(e,1);
    node2 = N(e,2);
    node3 = N(e,3);
    node4 = N(e,4);
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
    a = [a1, a2, a3, a4];
    b = [b1, b2, b3, b4];  c = [c1, c2, c3, c4];  d = [d1, d2, d3, d4];
    V(e) = (1/6)*det([1 1 1 1; x1 x2 x3 x4; y1 y2 y3 y4; z1 z2 z3 z4]);
    delta(e) = 0.5*(b1*c2 - b2*c1);
    % Evaluate K^e matrix (element K matrix)
    for i = 1:1:4
        for j = 1:1:4
            if (i == j)
                delij = 1;
            else
                delij = 0;
            end
            Klocal(e,i,j) = (alpha_x(e)*b(i)*b(j)+ ...
            alpha_y(e)*c(i)*c(j) + alpha_z(e)*d(i)*d(j))/((36*V(e)))+ ...
            V(e)*beta(e)*(1+delij)/20;
        end
    end
    % Evaluate b_tilde^e vector
    for i = 1:1:4
        be_tilde(e,i) = E_0*k_0*k_0*(1/mu_r(e) - eps_r(e))*V(e)/4;
    end
end
for e=1:nElements
    % Assemble K^e into Kglobal and b^e into b
    for i = 1:1:4
        for j = 1:1:4
            Kglobal(N(e,i),N(e,j)) = Kglobal(N(e,i),N(e,j)) + Klocal(e,i,j);
        end
        btilde(N(e,i)) = btilde(N(e,i)) + be_tilde(e,i);
    end
end
% Construct Ks
Ks=zeros(nSegABC,3,3);
bs_tilde=zeros(nSegABC,3);
for s = 1:1:nSegABC
    x1 = x(NS(s,1)); x2 = x(NS(s,2)); x3 = x(NS(s,3));
    y1 = y(NS(s,1)); y2 = y(NS(s,2)); y3 = y(NS(s,3));
    z1 = z(NS(s,1)); z2 = z(NS(s,2)); z3 = z(NS(s,3));
    
    del(s) = 0.5*sqrt((x2*y1 - x3*y1 -x1*y2 + x3*y2 + x1*y3 - x2*y3)^2 ...
        +(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)^2 ...
        +(y2*z1 - y3*z1 - y1*z2 + y3*z2 + y1*z3 - y2*z3)^2);
    
    for ii = 1:1:3
        fun2 = @(ksi,ksj)bstilde1(ksi,ksj,ii,alpha,kappa,k_0,x1,x2,x3,y1,y2,y3,del(s),E_0);
        bs_tilde(s,ii) = integral2(fun2,0,1,0,1);
        for jj = 1:1:3
            if (ii == jj)
                Ks(s,ii,jj) = gamma_1*del(s)/6;
            else
                Ks(s,ii,jj) = gamma_1*del(s)/12;
            end
        end
    end
end
% Assemble Ks into K
for s = 1:1:nSegABC
    for i = 1:1:3
        for j = 1:1:3
            Kglobal(NS(s,i),NS(s,j)) = Kglobal(NS(s,i),NS(s,j)) + Ks(s,i,j);
        end
    end
end
% Assemble bs into b
for s=1:1:nSegABC
    for i=1:3
        btilde(NS(s,i)) = btilde(NS(s,i)) + bs_tilde(s,i);
    end
end
% Impose Dirichlet Boundary COndition
for i=1:1:nDBC
    for j=1:1:nNodes
        if (j ~= ABC(i))
            btilde(j) = btilde(j) - Kglobal(j, DBC(i))*0;
        end
    end
    Kglobal(:,DBC(i)) = 0;
    Kglobal(DBC(i),:) = 0;
    Kglobal(DBC(i),DBC(i)) = 1;
    btilde(DBC(i)) = 0;
end
% Solution of global matrix system
Ez = Kglobal\btilde;
figure
pdeplot3D(model,'ColorMapData',abs(Ez))

% Generate solution over a grid and plot it
% [xgrid, ygrid, zgrid] = meshgrid(-rABC:0.01*(2*rABC):rABC, -rABC:0.01*(2*rABC):rABC, 0:0.1:len);
[xgrid, ygrid, zgrid] = meshgrid(-rABC:0.1*(2*rABC):rABC, -rABC:0.1*(2*rABC):rABC, 0:0.099*len:len);
% Ezgrid = zeros(101,101,101);
Ezgrid = zeros(11,11,11);
for i = 1:1:11
    for j = 1:1:11
        for k = 1:1:11
            for e = 1:1:nElements
                x1p = x(N(e,1)) - xgrid(i,j,k);
                x2p = x(N(e,2)) - xgrid(i,j,k);
                x3p = x(N(e,3)) - xgrid(i,j,k);
                x4p = x(N(e,4)) - xgrid(i,j,k);
                y1p = y(N(e,1)) - ygrid(i,j,k);
                y2p = y(N(e,2)) - ygrid(i,j,k);
                y3p = y(N(e,3)) - ygrid(i,j,k);
                y4p = y(N(e,4)) - ygrid(i,j,k);
                z1p = z(N(e,1)) - zgrid(i,j,k);
                z2p = z(N(e,2)) - zgrid(i,j,k);
                z3p = z(N(e,3)) - zgrid(i,j,k);
                z4p = z(N(e,4)) - zgrid(i,j,k);
                A1 = 0.5*abs(x2p*(y3p*z4p - y4p*z3p) + x3p*(y4p*z2p - y2p*z4p) +...
                    x4p*(y2p*z3p - y3p*z2p));
                x1p = x(N(e,1)) - xgrid(i,j,k);
                y1p = y(N(e,1)) - ygrid(i,j,k);
                A2 = 0.5*abs(x1p*(y2p*z3p - y3p*z2p) + x2p*(y3p*z1p - y1p*z3p) +...
                    x3p*(y1p*z2p - y2p*z1p));
                z1p = z(N(e,1)) - zgrid(i,j,k);
                A3 = 0.5*abs(x3p*(y4p*z1p - y1p*z4p) + x4p*(y1p*z3p - y3p*z1p) +...
                    x1p*(y3p*z4p - y4p*z3p));
                A4 = 0.5*abs(x4p*(y1p*z2p - y2p*z1p) + x1p*(y2p*z4p - y4p*z2p) +...
                    x2p*(y4p*z1p - y1p*z4p));
                x21 = x(N(e,2)) - x(N(e,1));
                x31 = x(N(e,3)) - x(N(e,1));
                x41 = x(N(e,4)) - x(N(e,1));
                y21 = y(N(e,2)) - y(N(e,1));
                y31 = y(N(e,3)) - y(N(e,1));
                y41 = y(N(e,4)) - y(N(e,1));
                z21 = z(N(e,2)) - z(N(e,1));
                z31 = z(N(e,3)) - z(N(e,1));
                z41 = z(N(e,4)) - z(N(e,1));
                node1 = N(e,1);
                node2 = N(e,2);
                node3 = N(e,3);
                node4 = N(e,4);
                x1 = nodeIndex(1, node1);  y1 = nodeIndex(2, node1);  z1 = nodeIndex(3, node1);
                x2 = nodeIndex(1, node2);  y2 = nodeIndex(2, node2);  z2 = nodeIndex(3, node2);
                x3 = nodeIndex(1, node3);  y3 = nodeIndex(2, node3);  z3 = nodeIndex(3, node3);
                x4 = nodeIndex(1, node4);  y4 = nodeIndex(2, node4);  z4 = nodeIndex(3, node4);
                V(e) = (1/6)*det([1 1 1 1; x1 x2 x3 x4; y1 y2 y3 y4; z1 z2 z3 z4]);
                 if abs(V(e) - (A1+A2+A3)) < 0.01*V(e)
                     ksi = (x31*(ygrid(i,j,k) - y(N(e,1)))*(zgrid(i,j,k) - z(N(e,4))) - ...
                         x31*(ygrid(i,j,k) - y(N(e,4)))*(zgrid(i,j,k) - z(N(e,1))) + ...
                         y31*(xgrid(i,j,k) - x(N(e,1)))*(zgrid(i,j,k) - z(N(e,4))) - ...
                         y31*(xgrid(i,j,k) - x(N(e,4)))*(zgrid(i,j,k) - z(N(e,1))) + ...
                         z31*(xgrid(i,j,k) - x(N(e,1)))*(ygrid(i,j,k) - y(N(e,4))) - ...
                         z31*(xgrid(i,j,k) - x(N(e,4)))*(ygrid(i,j,k) - y(N(e,1))))*V(e)/4;
                     ita = (x21*(ygrid(i,j,k) - y(N(e,1)))*(zgrid(i,j,k)-z(N(e,4))) - ...
                         y21*(xgrid(i,j,k) - x(N(e,1)))*(zgrid(i,j,k) - z(N(e,4))) - ...
                         y21*(xgrid(i,j,k) - x(N(e,4)))*(zgrid(i,j,k) - z(N(e,1))) + ...
                         z21*(xgrid(i,j,k) - x(N(e,1)))*(ygrid(i,j,k) - y(N(e,4))) - ...
                         z21*(xgrid(i,j,k) - x(N(e,4)))*(ygrid(i,j,k) - y(N(e,1))))*V(e)/4;
                     zia = (x41*(ygrid(i,j,k) - y(N(e,2)))*(zgrid(i,j,k) - z(N(e,3))) - ...
                         x41*(ygrid(i,j,k) - y(N(e,3)))*(zgrid(i,j,k) - z(N(e,2))) + ...
                         y41*(zgrid(i,j,k) - z(N(e,3)))*(xgrid(i,j,k) - x(N(e,2))) - ...
                         y41*(zgrid(i,j,k) - z(N(e,2)))*(xgrid(i,j,k) - x(N(e,3))) + ...
                         z41*(ygrid(i,j,k) - y(N(e,3)))*(xgrid(i,j,k) - x(N(e,2))) - ...
                         z41*(ygrid(i,j,k) - y(N(e,2)))*(xgrid(i,j,k) - x(N(e,3))))*V(e)/4;
                     N1 = 1 - ksi - ita - zia;
                     N2 = ksi;
                     N3 = ita;
                     N4 = zia;
                     Ezgrid(i,j,k) = N1*Ez(N(e,1)) + N2*Ez(N(e,2)) + N3*Ez(N(e,3)) + N4*Ez(N(e,4));
                end
            end
        end
    end
end
% Display contour plot of FEM solution
xgrid1 = xgrid(:,:,1);
ygrid1 = ygrid(:,:,1);
Ezgrid1 = Ezgrid(:,:,1);
figure;
contourf(xgrid1,ygrid1,abs(Ezgrid1));
xlabel('x (wavelengths)');
ylabel('y (wavelengths)');
title('Total Electric Field ? Contour Plot');
colorbar;
% Evaluate exact solution at rho = (rCyl + rABC)/2
d2p=pi/180;
dist=lambda; % where to evaluate the field
% Ez_eval=zeros(1,1441);
% phi=zeros(1,1441);
Ez_eval=zeros(1,11);
phi=zeros(1,11);
for I=1:11 %721
    phi(I)=(I-1)*0.25; %0.5;
    xeval=dist*cos(phi(I)*d2p);
    yeval=dist*sin(phi(I)*d2p);
    zeval=dist;
    for e=1:nElements
        x1p = x(N(e,1)) - xeval;
        x2p = x(N(e,2)) - xeval;
        x3p = x(N(e,3)) - xeval;
        x4p = x(N(e,4)) - xeval;
        y1p = y(N(e,1)) - yeval;
        y2p = y(N(e,2)) - yeval;
        y3p = y(N(e,3)) - yeval;
        y4p = y(N(e,4)) - yeval;
        z1p = z(N(e,1)) - zeval;
        z2p = z(N(e,2)) - zeval;
        z3p = z(N(e,3)) - zeval;
        z4p = z(N(e,4)) - zeval;
        A1 = 0.5*abs(x2p*(y3p*z4p - y4p*z3p) + x3p*(y4p*z2p - y2p*z4p) +...
            x4p*(y2p*z3p - y3p*z2p));
        A2 = 0.5*abs(x1p*(y2p*z3p - y3p*z2p) + x2p*(y3p*z1p - y1p*z3p) +...
            x3p*(y1p*z2p - y2p*z1p));
        A3 = 0.5*abs(x3p*(y4p*z1p - y1p*z4p) + x4p*(y1p*z3p - y3p*z1p) +...
            x1p*(y3p*z4p - y4p*z3p));
        A4 = 0.5*abs(x4p*(y1p*z2p - y2p*z1p) + x1p*(y2p*z4p - y4p*z2p) +...
            x2p*(y4p*z1p - y1p*z4p));
        x21 = x(N(e,2)) - x(N(e,1));
        x31 = x(N(e,3)) - x(N(e,1));
        x41 = x(N(e,4)) - x(N(e,1));
        y21 = y(N(e,2)) - y(N(e,1));
        y31 = y(N(e,3)) - y(N(e,1));
        y41 = y(N(e,4)) - y(N(e,1));
        z21 = z(N(e,2)) - z(N(e,1));
        z31 = z(N(e,3)) - z(N(e,1));
        z41 = z(N(e,4)) - z(N(e,1));
        node1 = N(e,1);
        node2 = N(e,2);
        node3 = N(e,3);
        node4 = N(e,4);
        x1 = nodeIndex(1, node1);  y1 = nodeIndex(2, node1);  z1 = nodeIndex(3, node1);
        x2 = nodeIndex(1, node2);  y2 = nodeIndex(2, node2);  z2 = nodeIndex(3, node2);
        x3 = nodeIndex(1, node3);  y3 = nodeIndex(2, node3);  z3 = nodeIndex(3, node3);
        x4 = nodeIndex(1, node4);  y4 = nodeIndex(2, node4);  z4 = nodeIndex(3, node4);
        V(e) = (1/6)*det([1 1 1 1; x1 x2 x3 x4; y1 y2 y3 y4; z1 z2 z3 z4]);
        if abs(V(e) - (A1+A2+A3)) < 0.01*V(e)
            ksi = (x31*(yeval - y(N(e,1)))*(zeval - z(N(e,4))) - ...
                 x31*(yeval - y(N(e,4)))*(zeval - z(N(e,1))) + ...
                 y31*(xeval - x(N(e,1)))*(zeval - z(N(e,4))) - ...
                 y31*(xeval - x(N(e,4)))*(zeval - z(N(e,1))) + ...
                 z31*(xeval - x(N(e,1)))*(yeval - y(N(e,4))) - ...
                 z31*(xeval - x(N(e,4)))*(yeval - y(N(e,1))))*V(e)/4;
            ita = (x21*(yeval - y(N(e,1)))*(zeval - z(N(e,4))) - ...
                 y21*(xeval - x(N(e,1)))*(zeval - z(N(e,4))) - ...
                 y21*(xeval - x(N(e,4)))*(zeval - z(N(e,1))) + ...
                 z21*(xeval - x(N(e,1)))*(yeval - y(N(e,4))) - ...
                 z21*(xeval - x(N(e,4)))*(yeval - y(N(e,1))))*V(e)/4;
             zia = (x41*(yeval - y(N(e,2)))*(zeval - z(N(e,3))) - ...
                 x41*(yeval - y(N(e,3)))*(zeval - z(N(e,2))) + ...
                 y41*(zeval - z(N(e,3)))*(xeval - x(N(e,2))) - ...
                 y41*(zeval - z(N(e,2)))*(xeval - x(N(e,3))) + ...
                 z41*(yeval - y(N(e,3)))*(xeval - x(N(e,2))) - ...
                 z41*(yeval - y(N(e,2)))*(xeval - x(N(e,3))))*V(e)/4;
            N1 = 1 - ksi - ita - zia;
            N2 = ksi;
            N3 = ita;
            N4 = zia;
            Ez_eval(I)=N1*Ez(N(e,1))+N2*Ez(N(e,2))+N3*Ez(N(e,3))+N4*Ez(Ne,4);
        end
    end
end
% Plot the analytical solution and the FEM solution at a distance 'dist'
figure;
plot(phi,abs(Ez_eval),'b--'),legend('FEM (1^{st} order ABC)');
xlabel('Angle (degrees)');
ylabel('Electric Field (V/m)');
% axis([0 360 0 2*E_0]);
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evaluate bsilde1(e,i) for the FIRST ORDER ABC (1 = FIRST ORDER ABC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f1 = bstilde1(ksi,ksj,m,alpha,kappa,k_0,x1,x2,x3,y1,y2,y3,L_s,E_0)
Ni = ksi;
Nj = ksj;
if (m == 1)
    Ni = 1 - ksi;
    Nj = 1 - ksj;
end
J = complex(0.,1.);
x = (1 - ksi).*x1 + ksi.*x2 + x3;
y = (1 - ksj).*y1 + ksj.*y2 + y3;
expfac = exp(-J*k_0.*x -J*k_0.*y);
f1 = alpha*(J*k_0 + 0.5*kappa - J*k_0*kappa.*x - J*k_0*kappa.*y).*...
    Ni.*Nj.*L_s*E_0.*expfac;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evaluate btilde(e,i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f2 = betilde(u,v,w,a,b,c,d,k_0,x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4)
J = complex(0.,1.);
x = (1 - u).*x1 + u.*((1 - v).*x2 + v.*((1 - w).*x3 + w.*x4));
y = (1 - u).*y1 + u.*((1 - v).*y2 + v.*((1 - w).*y3 + w.*y4));
z = (1 - u).*z1 + u.*((1 - v).*z2 + v.*((1 - w).*z3 + w.*z4));
jacob = (-x1+(1 - v).*x2+(1-w).*x3+w.*x4).*(-u.*y2+u.*y3) - ...
        (-y1+(1 - v).*y2+v.*y3).*(-u.*x2+u.*x3);
f2 = (a+b.*x+c.*y+d.*z).*exp(-J*k_0.*x).*jacob;
end
% simple linear array search
function index = ArrSearch(arr,val)
res = -1;
for ii=1:length(arr)
    if (arr(ii)==val)
        res=ii;
    break;
    end
end
index=res;
end