lambda = 1;
rCyl = 0.5*lambda; % Radius of the circular PEC cylinder in wavelengths
rABC = 1.5*lambda; % Radius of the outer circular ABC boundary in wavelengths
h = 0.04*lambda; % Discretization size in wavelengths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCT TRIANGULAR MESH
% Constructs a mesh of the circular solution domain with as many
% equilateral triangles as possible. The rest of the domain is
% filled with skewed triangles.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nRing = ceil((rABC - rCyl)/(sqrt(3)*h/2)); % No. of annular rings
dr = (rABC - rCyl)/nRing; % Radial step size
% Define nodes on the periphery of each annular ring
nNodes = 0; % Total number of nodes
for j = 1:1:(nRing+1)
    r = rCyl + (j - 1)*dr; % determine the radius of the j^th annular ring
    nPoints = ceil(2*pi*r/h); % no. of nodes on the bdary of the j^th ring
    dphi = 2*pi/nPoints; % incremental angle from one node to the next
    for jj = 1:1:nPoints
        nNodes = nNodes + 1; % increment total number of nodes
        phi = (jj - 1)*dphi; % absolute angular coordinate of current
        x(nNodes) = r*cos(phi); % x?coordinate of current node
        y(nNodes) = r*sin(phi); % y?coordinate of current node
    end
end

% Triangulate using Delaunay Triangulation
tri = delaunay(x, y);
% At this point, triplot(tri,x,y) shows that the region within the PEC
% has also been meshed, so we must eliminate the triangles in this region.
% Eliminate triangles within PEC
nElements = 0; % No of elements
% length(tri) = size(tri,1) = number of elements
% (before removal of elements interior to PEC)
for j = 1:1:length(tri)
    rem = 0;
    for jj = 1:1:3
        r0 = sqrt(x(tri(j,jj))^2 + y(tri(j,jj))^2); % determine radius of pt
        if (abs(rCyl - r0) < 0.0001*rCyl) % check for interioricity
            rem = rem + 1;
        end
    end
    if (rem <= 2)
        nElements = nElements + 1; % increment number of (valid) elements
        N(nElements,:) = tri(j,:); % element connectivity array n(e, i)
    end
end
% N(e, i) is the connectivity array
% The global node coordinates are stored in x and y
% To access the global node coordinates of node i of element e, use
% X = x(N(e,i)), Y = y(N(e,i))
% Delaunay Triangulation does not order nodes in a CCW orientation as
% required by the FEM algorithm. So, we have to reorient the local
% node numbering. The criterion is that if a triangle has a negative area
% computed using the expression for area when the nodes are CCW oriented,
% then the nodes must be re?numbered.
for e=1:1:nElements
    x21 = x(N(e,2)) - x(N(e,1));
    x31 = x(N(e,3)) - x(N(e,1));
    y21 = y(N(e,2)) - y(N(e,1));
    y31 = y(N(e,3)) - y(N(e,1));
    Ae = 0.5*(x21*y31 - x31*y21);
    if Ae < 0
        % simple swap of node numbers
        temp = N(e, 2);
        N(e, 2) = N(e, 3);
        N(e, 3) = temp;
    end
end
% Plot the mesh
figure;
triplot(N, x, y);
title('Mesh');
xlabel('x (wavelengths)');
ylabel('y (wavelengths)');
axis([-rABC rABC -rABC rABC]);
axis square;


% This program solves for the total electric field in a region bounded
% by a circular boundary, using the Nodal 2D Finite Element Method.
% The 1st order ABC is imposed on the boundary.
% It is assumed on entering this code that the mesh is well defined and
% p(:,1) = x?coordinates of vertices
% p(:,2) = y?coordinates of vertices
% N = node connectivity matrix (properly ordered)
% Define x and y coordinate vectors for use below
% x = p(:,1);
% y = p(:,2);

% Lets find boundary edges
p = [x ; y]; 
boundedges = boundedges(p',N);
% boundedges = boundedges';
hold on;
for i=1:length(boundedges)
    x1=x(boundedges(i,1)); x2=x(boundedges(i,2));
    y1=y(boundedges(i,1)); y2=y(boundedges(i,2));
    line([x1 x2],[y1 y2],'Color','b','LineWidth',1.5); % Plot the boundary edges
end
nDBC = 0; % no of Dirichlet Boundary Condition nodes
nABC = 0; % no of ABC nodes
nInt = 0; % no of internal nodes
nNodes=length(x); % total no of nodes
DBC(1) = 0;
% Determine the ABC and DBC nodes
for jj=1:1:nNodes
    r0 = sqrt(x(jj)^2+y(jj)^2);
    if(abs(rABC - r0)<0.0001*rABC) % ABC will be imposed on a circle always
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
        else % internal node
            nInt=nInt+1;
            INT(nInt)=jj;
        end
    end
end
% flag each node
% flag = 1 for ABC node
% = 0 for internal node
% = ?-1 for DBC node
flag = zeros(nElements,3);
for e=1:nElements
    for i=1:3
        index=ArrSearch(ABC,N(e,i));
        if (index ~= -1)
            flag(e,i) = 1; % ABC node
        else
            index=ArrSearch(DBC,N(e,i));
            if(index ~= -1)
                flag(e,i) = -1; % DBC node
            else
                flag(e,i) = 0; % internal node
            end
        end
    end
end
% Create NS(s,i) array
nSegABC = 0;
for e=1:1:nElements
    if flag(e,1) == 1 && flag(e,2) == 1
        nSegABC = nSegABC + 1;
        NS(nSegABC, 1) = N(e,1);
        NS(nSegABC, 2) = N(e,2);
    elseif flag(e,1) == 1 && flag(e,3) == 1
        nSegABC = nSegABC + 1;
        NS(nSegABC, 1) = N(e,1);
        NS(nSegABC, 2) = N(e,3);
    elseif flag(e,2) == 1 && flag(e,3) == 1
        nSegABC = nSegABC + 1;
        NS(nSegABC, 1) = N(e,2);
        NS(nSegABC, 2) = N(e,3);
    end
end
fprintf(' RESULTS \n\n');
fprintf('The number of nodes = %6i\n', nNodes);
fprintf('The number of elements = %6i\n', nElements);
fprintf('The number of nodes on a Dirichlet boundary = %6i\n', nDBC);
fprintf('The number of nodes on an ABC boundary = %6i\n', nABC);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = 1;
rABC = 1.5*lambda; % Radius of the outer circular ABC boundary in wavelengths
h = 0.04*lambda; % Discretization size in wavelengths
% E_0 = 1; % Amplitude of incident electric field (V/m)
E_0 = 0.8*exp(-1i*0.8); % Amplitude of incident electric field (V/m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE CONSTANTS FOR EACH ELEMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
beta = zeros(1,nElements); % coefficient of term linear in field in PDE
A = zeros(1,nElements); % area of element e
for e=1:1:nElements
    mu_r(e) = mur;
    eps_r(e) = epsr;
    alpha_x(e) = -1/mur;
    alpha_y(e) = -1/mur;
    beta(e) = (k_0)*(k_0)*eps_r(e);
end
% Initialize global K matrix and the right hand side vector b
Kglobal = zeros(nNodes, nNodes);
btilde = zeros(nNodes, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FORMATION OF ELEMENT MATRICES AND ASSEMBLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coefficients of the linear interpolation functions
a = zeros(nElements,3);
b = zeros(nElements,3);
c = zeros(nElements,3);
be_tilde=zeros(1,nElements);
Klocal=zeros(nElements,3,3);
for e=1:1:nElements
    x_21 = x(N(e,2)) - x(N(e,1));
    x_31 = x(N(e,3)) - x(N(e,1));
    x_32 = x(N(e,3)) - x(N(e,2));
    x_13 = -x_31;
    y_12 = y(N(e,1)) - y(N(e,2));
    y_21 = -y_12;
    y_31 = y(N(e,3)) - y(N(e,1));
    y_23 = y(N(e,2)) - y(N(e,3));
    A(e) = 0.5*(x_21*y_31 - x_31*y_21);
    x1 = x(N(e,1)); y1 = y(N(e,1));
    x2 = x(N(e,2)); y2 = y(N(e,2));
    x3 = x(N(e,3)); y3 = y(N(e,3));
    a(e,1) = x2*y3 - x3*y2;
    b(e,1) = y2 - y3;
    c(e,1) = x3 - x2;
    a(e,2) = x3*y1 - x1*y3;
    b(e,2) = y3 - y1;
    c(e,2) = x1 - x3;
    a(e,3) = x1*y2 - x2*y1;
    b(e,3) = y1 - y2;
    c(e,3) = x2 - x1;
    % Evaluate K^e matrix (element K matrix)
    for i = 1:1:3
        for j = 1:1:3
            if (i == j)
                delij = 1;
            else
                delij = 0;
            end
            Klocal(e,i,j) = (alpha_x(e)*b(e,i)*b(e,j)+ ...
            alpha_y(e)*c(e,i)*c(e,j))/(4*A(e))+ ...
            A(e)*beta(e)*(1+delij)/12;
        end
    end
    % Evaluate b_tilde^e vector
    for i = 1:1:3
        fun = @(u,v)betilde(u,v,a(e,i),b(e,i),c(e,i),k_0,x1,x2,x3,y1,y2,y3);
        % use for scattered field formulation (change eps and
        % mu for dielectrics, etc.)
        be_tilde(e,i) = (E_0*k_0*k_0*(1/mu_r(e) - eps_r(e))/(2*A(e)))*...
        integral2(fun,0,1,0,1);
    end
end
for e=1:nElements
    % Assemble K^e into Kglobal and b^e into b
    for i = 1:1:3
        for j = 1:1:3
            Kglobal(N(e,i),N(e,j)) = Kglobal(N(e,i),N(e,j)) + Klocal(e,i,j);
        end
        btilde(N(e,i)) = btilde(N(e,i)) + be_tilde(e,i);
    end
end
% Construct Ks
Ks=zeros(nSegABC,2,2);
bs_tilde=zeros(nSegABC,2);
for s = 1:1:nSegABC
    xdiff = x(NS(s,1)) - x(NS(s,2));
    ydiff = y(NS(s,1)) - y(NS(s,2));
    L_s = sqrt(xdiff^2 + ydiff^2);
    x1 = x(NS(s,1)); x2 = x(NS(s,2));
    y1 = y(NS(s,1)); y2 = y(NS(s,2));
    for ii = 1:1:2
        fun2 = @(ksi)bstilde1(ksi,ii,alpha,kappa,k_0,x1,x2,L_s,E_0);
        bs_tilde(s,ii) = integral(fun2,0,1);
        for jj = 1:1:2
            if (ii == jj)
                Ks(s,ii,jj) = gamma_1*L_s/3;
            else
                Ks(s,ii,jj) = gamma_1*L_s/6;
            end
        end
    end
end
% Assemble Ks into K
for s = 1:1:nSegABC
    for i = 1:1:2
        for j = 1:1:2
            Kglobal(NS(s,i),NS(s,j)) = Kglobal(NS(s,i),NS(s,j)) + Ks(s,i,j);
        end
    end
end
% Assemble bs into b
for s=1:1:nSegABC
    for i=1:2
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
% Generate solution over a grid and plot it
[xgrid, ygrid] = meshgrid(-rABC:0.01*(2*rABC):rABC, -rABC:0.01*(2*rABC):rABC);
Ezgrid = zeros(101,101);

for i = 1:1:101
    for j = 1:1:101
        for e = 1:1:nElements
            x2p = x(N(e,2)) - xgrid(i,j);
            x3p = x(N(e,3)) - xgrid(i,j);
            y2p = y(N(e,2)) - ygrid(i,j);
            y3p = y(N(e,3)) - ygrid(i,j);
            A1 = 0.5*abs(x2p*y3p - x3p*y2p);
            x2p = x(N(e,2)) - xgrid(i,j);
            x1p = x(N(e,1)) - xgrid(i,j);
            y2p = y(N(e,2)) - ygrid(i,j);
            y1p = y(N(e,1)) - ygrid(i,j);
            A2 = 0.5*abs(x2p*y1p - x1p*y2p);
            x1p = x(N(e,1)) - xgrid(i,j);
            x3p = x(N(e,3)) - xgrid(i,j);
            y1p = y(N(e,1)) - ygrid(i,j);
            y3p = y(N(e,3)) - ygrid(i,j);
            A3 = 0.5*abs(x1p*y3p - x3p*y1p);
            x21 = x(N(e,2)) - x(N(e,1));
            x31 = x(N(e,3)) - x(N(e,1));
            y21 = y(N(e,2)) - y(N(e,1));
            y31 = y(N(e,3)) - y(N(e,1));
            Ae = 0.5*abs(x21*y31 - x31*y21);
            if abs(Ae - (A1+A2+A3)) < 0.00001*Ae
                ksi = (y31*(xgrid(i,j) - x(N(e,1))) - x31*...
                    (ygrid(i,j) - y(N(e,1))))/(2*Ae);
                ita = (-y21*(xgrid(i,j) - x(N(e,1)))+...
                    x21*(ygrid(i,j) - y(N(e,1))))/(2*Ae);
                N1 = 1 - ksi - ita;
                N2 = ksi;
                N3 = ita;
                Ezgrid(i,j) = N1*Ez(N(e,1)) + N2*Ez(N(e,2)) + N3*Ez(N(e,3));
            end
        end
    end
end
% Display contour plot of FEM solution
figure;
contourf(xgrid,ygrid,abs(Ezgrid));
xlabel('x (wavelengths)');
ylabel('y (wavelengths)');
axis([-rABC rABC -rABC rABC]);
title('Total Electric Field ? Contour Plot');
axis square;
colorbar;
% Evaluate exact solution at rho = (rCyl + rABC)/2
d2p=pi/180;
dist=lambda; % where to evaluate the field
Ez_eval=zeros(1,1441);
phi=zeros(1,1441);
for I=1:1441 %721
    phi(I)=(I-1)*0.25; %0.5;
    xeval=dist*cos(phi(I)*d2p);
    yeval=dist*sin(phi(I)*d2p);
    for e=1:nElements
        x2p=x(N(e,2)) - xeval;
        x3p=x(N(e,3)) - xeval;
        y2p=y(N(e,2)) - yeval;
        y3p=y(N(e,3)) - yeval;
        A1=0.5*abs(x2p*y3p - x3p*y2p);
        x2p=x(N(e,2)) - xeval;
        x1p=x(N(e,1)) - xeval;
        y2p=y(N(e,2)) - yeval;
        y1p=y(N(e,1)) - yeval;
        A2=0.5*abs(x2p*y1p - x1p*y2p);
        x1p=x(N(e,1)) - xeval;
        x3p=x(N(e,3)) - xeval;
        y1p=y(N(e,1)) - yeval;
        y3p=y(N(e,3)) - yeval;
        A3=0.5*abs(x1p*y3p - x3p*y1p);
        x21=x(N(e,2)) - x(N(e,1));
        x31=x(N(e,3)) - x(N(e,1));
        y21=y(N(e,2)) - y(N(e,1));
        y31=y(N(e,3)) - y(N(e,1));
        Ae=0.5*(x21*y31 - x31*y21);
        if abs(Ae - (A1+A2+A3)) < 0.00001*Ae
            ksi=(y31*(xeval - x(N(e,1))) - x31*(yeval - y(N(e,1))))/(2*Ae);
            ita=(-y21*(xeval - x(N(e,1)))+x21*(yeval - y(N(e,1))))/(2*Ae);
            N1=1 - ksi - ita;
            N2=ksi;
            N3=ita;
            Ez_eval(I)=N1*Ez(N(e,1))+N2*Ez(N(e,2))+N3*Ez(N(e,3));
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

figure;
plot(phi,(((abs(Ez_eval)).^2)/377)*10^3,'b--'),legend('FEM (1^{st} order ABC)');
xlabel('Angle (degrees)');
ylabel('Power Density (mW/m2)');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evaluate bsilde1(e,i) for the FIRST ORDER ABC (1 = FIRST ORDER ABC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f1 = bstilde1(ksi,m,alpha,kappa,k_0,x1,x2,L_s,E_0)
N = ksi;
if (m == 1)
    N = 1 - ksi;
end
J = complex(0.,1.);
x = (1 - ksi).*x1 + ksi.*x2;
expfac = exp(-J*k_0.*x);
f1 = alpha*(J*k_0 + 0.5*kappa - J*k_0*kappa.*x).*N.*L_s*E_0.*expfac;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evaluate btilde(e,i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f2 = betilde(u,v,a,b,c,k_0,x1,x2,x3,y1,y2,y3)
J = complex(0.,1.);
x = (1 - u).*x1 + u.*((1 - v).*x2 + v.*x3);
y = (1 - u).*y1 + u.*((1 - v).*y2 + v.*y3);
jacob = (-x1+(1 - v).*x2+v.*x3).*(-u.*y2+u.*y3) - ...
        (-y1+(1 - v).*y2+v.*y3).*(-u.*x2+u.*x3);
f2 = (a+b.*x+c.*y).*exp(-J*k_0.*x).*jacob;
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