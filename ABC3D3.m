lambda = .5/6.35; %1;
rCyl = 1*lambda; % Radius of the circular PEC cylinder in wavelengths
rABC = 2; %1.5*lambda; % Radius of the outer circular ABC boundary in wavelengths
len = 3*lambda; % Length of the waveguide
%radius =1 %radius =6.35
%height =1
%lambda =.5/6.35

% numberOfPDE = 1;
% [xg, yg] = meshgrid(-2:0.5:2);
% xg = xg(:);
% yg = yg(:);
% t = (pi/24:pi/12:2*pi)';
% x = cos(t);
% y = sin(t);
% circShp = alphaShape(x,y,2);
% in = inShape(circShp,xg,yg);
% xg = [xg(~in); cos(t)];
% yg = [yg(~in); sin(t)];
% zg = ones(numel(xg),1);
% xg = repmat(xg,5,1);
% yg = repmat(yg,5,1);
% zg = zg*(0:.25:1);
% zg = zg(:);
% 
% [xg1, yg1] = meshgrid(-2:0.5:2);
% xg1 = xg1(:);
% yg1 = yg1(:);
% zg1 = ones(numel(xg1),1);
% xg1 = repmat(xg1,5,1);
% yg1 = repmat(yg1,5,1);
% zg1 = zg1*(-1.25:.25:-0.25);
% zg1 = zg1(:);
% [xg2, yg2] = meshgrid(-2:0.5:2);
% xg2 = xg2(:);
% yg2 = yg2(:);
% zg2 = ones(numel(xg2),1);
% xg2 = repmat(xg2,5,1);
% yg2 = repmat(yg2,5,1);
% zg2 = zg2*(1.25:0.25:2.25);
% zg2 = zg2(:);
% xg3 = [xg1; xg; xg2];
% yg3 = [yg1; yg; yg2];
% zg3 = [zg1; zg; zg2];
% 
% shp = alphaShape(xg3,yg3,zg3);
% shp1 = alphaShape(xg,yg,zg);
% figure
% plot(shp)
% 
% [elements,nodes] = boundaryFacets(shp);
% nodes = nodes';
% elements = elements';
% model = createpde();
% geometryFromMesh(model,nodes,elements);
% 
% %View the geometry and face numbers.
% figure
% pdegplot(model,'FaceLabels','on','FaceAlpha',0.5)
% mes = generateMesh(model,'GeometricOrder','linear');
% 
% N = mes.Elements;
% nElements = length(N);
% nodeIndex = mes.Nodes;
% 
% % Lets find boundary edges
% p = [nodeIndex(1,:) ; nodeIndex(2,:) ; nodeIndex(3,:)]; 
% x = nodeIndex(1,:);
% y = nodeIndex(2,:);
% z = nodeIndex(3,:);
% N = N';
% boundedges = surftri(p',N);
% % boundedges = boundedges';
% hold on;
% for i=1:length(boundedges)
%     x1=x(boundedges(i,1)); x2=x(boundedges(i,2)); x3=x(boundedges(i,3));
%     y1=y(boundedges(i,1)); y2=y(boundedges(i,2)); y3=y(boundedges(i,3));
%     z1=z(boundedges(i,1)); z2=z(boundedges(i,2)); z3=z(boundedges(i,3));
%     plot3(x1,y1,z1,'-o','Color','b','MarkerSize',5); % Plot the boundary edges
%     %line([x1 x2 x3],[y1 y2 y3],[z1 z2 z3],'Color','r','LineWidth',2); % Plot the boundary edges
% end
% nDBC = 0; % no of Dirichlet Boundary Condition nodes
% nABC = 0; % no of ABC nodes
% nInt = 0; % no of internal nodes
% nNodes=length(nodeIndex); % total no of nodes
% DBC(1) = 0;
% % Determine the ABC and DBC nodes
% for jj=1:1:nNodes
%     if z(jj) < 0 || z(jj) > 1.25
%         nABC=nABC+1;
%         ABC(nABC)=jj;
%     else
%         r0 = sqrt(x(jj)^2+y(jj)^2);
%         if(abs(r0) > 1) %% && (z(jj) >= 0.5)
%             nABC=nABC+1;
%             ABC(nABC)=jj;
%         else % either an internal node or a node on the object (scatterer)
%         % if jj is an element of the be (bound edges) array
%         % then it is a node on the object hence a DBC node
%         % otherwise it is an internal node
%             index=ArrSearch(boundedges,jj);
%             if(index ~= -1) % DBC node
%                 nDBC=nDBC+1;
%                 DBC(nDBC)=jj;
%             else % internal node
%                 nInt=nInt+1;
%                 INT(nInt)=jj;
%             end
%         end
%     end
% end
% % flag each node
% % flag = 1 for ABC node
% % = 0 for internal node
% % = ?-1 for DBC node
% flag = zeros(nElements,4);
% for e=1:nElements
%     for i=1:4
%         index=ArrSearch(ABC,N(e,i));
%         if (index ~= -1)
%             flag(e,i) = 1; % ABC node
%         else
%             index=ArrSearch(DBC,N(e,i));
%             if(index ~= -1)
%                 flag(e,i) = -1; % DBC node
%              else
%                 flag(e,i) = 0; % internal node
%             end
%         end
%     end
% end
% % hold on;
% % for i=1:nDBC
% %     nodeNum = DBC(i);
% %     x1 = nodeIndex(1,nodeNum);
% %     y1 = nodeIndex(2,nodeNum);
% %     z1 = nodeIndex(3,nodeNum);
% %     plot3(x1,y1,z1,'-o','Color','b','MarkerSize',10); % Plot the boundary edges
% % end
% 
% % Create NS(s,i) array
% nSegABC = 0;
% for e=1:1:nElements
%     if flag(e,1) == 1 && flag(e,2) == 1 && flag(e,3) == 1
%         nSegABC = nSegABC + 1;
%         NS(nSegABC, 1) = N(e,1);
%         NS(nSegABC, 2) = N(e,2);
%         NS(nSegABC, 3) = N(e,3);
%     elseif flag(e,1) == 1 && flag(e,3) == 1 && flag(e,4) == 1
%         nSegABC = nSegABC + 1;
%         NS(nSegABC, 1) = N(e,1);
%         NS(nSegABC, 2) = N(e,3);
%         NS(nSegABC, 3) = N(e,4);
%     elseif flag(e,1) == 1 && flag(e,4) == 4 && flag(e,2) == 1
%         nSegABC = nSegABC + 1;
%         NS(nSegABC, 1) = N(e,1);
%         NS(nSegABC, 2) = N(e,4);
%         NS(nSegABC, 3) = N(e,2);
%     elseif flag(e,2) == 1 && flag(e,3) == 1 && flag(e,4) == 1
%         nSegABC = nSegABC + 1;
%         NS(nSegABC, 1) = N(e,2);
%         NS(nSegABC, 2) = N(e,3);
%         NS(nSegABC, 3) = N(e,4);
%     end
% end
% fprintf(' RESULTS \n\n');
% fprintf('The number of nodes = %6i\n', nNodes);
% fprintf('The number of elements = %6i\n', nElements);
% fprintf('The number of nodes on a Dirichlet boundary = %6i\n', nDBC);
% fprintf('The number of nodes on an ABC boundary = %6i\n', nABC);
% 
% E_0 = 10; % Amplitude of incident electric field (V/m)
% J = complex(0.,1.);
% mur = 1; % Relative permeability (assuming vaccum)
% epsr = 1; % Relative permittivity (asumming vacuum)
% k_0 = 2*pi/lambda; % wavenumber
% alpha = -1/mur;
% kappa = 1/rABC;
% gamma_1 = alpha*(J*k_0 + kappa/2 - J*kappa^2/(8*(J*kappa - k_0)));
% mu_r = zeros(1,nElements); % relative permeability
% eps_r = zeros(1,nElements); % relative permittivity
% alpha_x = zeros(1,nElements); % coefficient of second order x term in PDE
% alpha_y = zeros(1,nElements); % coefficient of second order y term in PDE
% alpha_z = zeros(1,nElements);
% beta = zeros(1,nElements); % coefficient of term linear in field in PDE
% A = zeros(1,nElements); % area of element e
% for e=1:1:nElements
%     mu_r(e) = mur;
%     eps_r(e) = epsr;
%     alpha_x(e) = -1/mur;
%     alpha_y(e) = -1/mur;
%     alpha_z(e) = -1/mur;
%     beta(e) = (k_0)*(k_0)*eps_r(e);
% end
% 
% Kglobal = zeros(nNodes, nNodes);
% btilde = zeros(nNodes, 1);
% a = zeros(nElements,4);
% b = zeros(nElements,4);
% c = zeros(nElements,4);
% be_tilde=zeros(1,nElements);
% Klocal=zeros(nElements,4,4);
% for e=1:1:nElements
%     node1 = N(e,1);
%     node2 = N(e,2);
%     node3 = N(e,3);
%     node4 = N(e,4);
%     x1 = nodeIndex(1, node1);  y1 = nodeIndex(2, node1);  z1 = nodeIndex(3, node1);
%     x2 = nodeIndex(1, node2);  y2 = nodeIndex(2, node2);  z2 = nodeIndex(3, node2);
%     x3 = nodeIndex(1, node3);  y3 = nodeIndex(2, node3);  z3 = nodeIndex(3, node3);
%     x4 = nodeIndex(1, node4);  y4 = nodeIndex(2, node4);  z4 = nodeIndex(3, node4);
%     a1 = x2*(y3*z4 - y4*z3) + x3*(y4*z2 - y2*z4) + x4*(y2*z3 - y3*z2);
%     a2 = x3*(y4*z1 - y1*z4) + x4*(y1*z3 - y3*z1) + x1*(y3*z4 - y4*z3);
%     a3 = x4*(y1*z2 - y2*z1) + x1*(y2*z4 - y4*z2) + x2*(y4*z1 - y1*z4);
%     a4 = x1*(y2*z3 - y3*z2) + x2*(y3*z1 - y1*z3) + x3*(y1*z2 - y2*z1);
%     b1 = (y3*z4 - z3*y4) + (y2*z4 - y4*z2) + (y2*z3 - z2*y3);  
%     b2 = (y3*z4 - z3*y4) + (y1*z4 - z1*y4) + (y1*z3 - z1*y3);  
%     b3 = (y2*z4 - z2*y4) + (y1*z4 - z1*y4) + (y1*z2 - z1*y2);
%     b4 = (y2*z3 - z2*y3) + (y1*z3 - z1*y3) + (y1*z2 - z1*y2);  
%     c1 = (x3*z4 - z3*x4) + (x2*z4 - x4*z2) + (x2*z3 - z2*x3);  
%     c2 = (x3*z4 - z3*x4) + (x1*z4 - z1*x4) + (x1*z3 - z1*x3);  
%     c3 = (x2*z4 - z2*x4) + (x1*z4 - z1*x4) + (x1*z2 - z1*x2);
%     c4 = (x2*z3 - z2*x3) + (x1*z3 - z1*x3) + (x1*z2 - z1*x2);   
%     d1 = (x3*y4 - y3*x4) + (x2*y4 - x4*y2) + (x2*y3 - y2*x3);  
%     d2 = (x3*y4 - y3*x4) + (x1*y4 - y1*x4) + (x1*y3 - y1*x3);  
%     d3 = (x2*y4 - y2*x4) + (x1*y4 - y1*x4) + (x1*y2 - y1*x2);
%     d4 = (x2*y3 - y2*x3) + (x1*y3 - y1*x3) + (x1*y2 - y1*x2);
%     a = [a1, a2, a3, a4]; b = [b1, b2, b3, b4];  
%     c = [c1, c2, c3, c4];  d = [d1, d2, d3, d4];
%     V(e) = (1/6)*det([1 1 1 1; x1 x2 x3 x4; y1 y2 y3 y4; z1 z2 z3 z4]);
%     % Evaluate K^e matrix (element K matrix)
%     for i = 1:1:4
%         for j = 1:1:4
%             if (i == j)
%                 delij = 1;
%             else
%                 delij = 0;
%             end
%             Klocal(e,i,j) = (alpha_x(e)*b(i)*b(j)+ ...
%             alpha_y(e)*c(i)*c(j) + alpha_z(e)*d(i)*d(j))/((36*V(e)))+ ...
%             V(e)*beta(e)*(1+delij)/20;
%         end
%     end
%     % Evaluate b_tilde^e vector
%     for i = 1:1:4
%         be_tilde(e,i) = E_0*k_0*k_0*(1/mu_r(e) - eps_r(e))*V(e)/4;
%     end
% end
% for e=1:nElements
%     % Assemble K^e into Kglobal and b^e into b
%     for i = 1:1:4
%         for j = 1:1:4
%             Kglobal(N(e,i),N(e,j)) = Kglobal(N(e,i),N(e,j)) + Klocal(e,i,j);
%         end
%         btilde(N(e,i)) = btilde(N(e,i)) + be_tilde(e,i);
%     end
% end
% % Construct Ks
% Ks=zeros(nSegABC,3,3);
% bs_tilde=zeros(nSegABC,3);
% for s = 1:1:nSegABC
%     x1 = x(NS(s,1)); x2 = x(NS(s,2)); x3 = x(NS(s,3));
%     y1 = y(NS(s,1)); y2 = y(NS(s,2)); y3 = y(NS(s,3));
%     z1 = z(NS(s,1)); z2 = z(NS(s,2)); z3 = z(NS(s,3));
%     
%     del(s) = (1/2)*det([x1 x2 x3; y1 y2 y3; z1 z2 z3]);
%    
%     for ii = 1:1:3
%         fun2 = @(ksi,ksj)bstilde1(ksi,ksj,ii,alpha,kappa,k_0,x1,x2,x3,y1,y2,y3,del(s),E_0);
%         bs_tilde(s,ii) = integral2(fun2,0,1,0,1);
%         bs_tilde(s,ii) = del(s)/3;
%         for jj = 1:1:3
%             if (ii == jj)
%                 Ks(s,ii,jj) = gamma_1*del(s)/6;
%             else
%                 Ks(s,ii,jj) = gamma_1*del(s)/12;
%             end
%         end
%     end
% end
% % Assemble Ks into K
% for s = 1:1:nSegABC
%     for i = 1:1:3
%         for j = 1:1:3
%             Kglobal(NS(s,i),NS(s,j)) = Kglobal(NS(s,i),NS(s,j)) + Ks(s,i,j);
%         end
%     end
% end
% % Assemble bs into b
% for s=1:1:nSegABC
%     for i=1:3
%         btilde(NS(s,i)) = btilde(NS(s,i)) + bs_tilde(s,i);
%     end
% end
% % Impose Dirichlet Boundary COndition
% p0 = 0.5;
% for i=1:1:nDBC
%     for j=1:1:nNodes
%         if (j ~= ABC(i))
%             btilde(j) = btilde(j) - Kglobal(j, DBC(i))*p0;
%         end
%     end
%     Kglobal(:,DBC(i)) = p0;
%     Kglobal(DBC(i),:) = p0;
%     Kglobal(DBC(i),DBC(i)) = 1;
%     btilde(DBC(i)) = p0;
% end
% % Solution of global matrix system
% Ez = Kglobal\btilde;
% 
% figure
% pdeplot3D(model,'ColorMapData',abs(Ez),'FaceAlpha',.3);

[xgrid, ygrid] = meshgrid(-1:0.01:1, -1:0.01:1);
Ezgrid = zeros(201,201);

zgrid=1;
for i = 1:1:201
    for j = 1:1:201
        for e = 1:nElements
            x1p = x(N(e,1)) - xgrid(i,j);
            x2p = x(N(e,2)) - xgrid(i,j);
            x3p = x(N(e,3)) - xgrid(i,j);
            x4p = x(N(e,4)) - xgrid(i,j);
            y1p = y(N(e,1)) - ygrid(i,j);
            y2p = y(N(e,2)) - ygrid(i,j);
            y3p = y(N(e,3)) - ygrid(i,j);
            y4p = y(N(e,4)) - ygrid(i,j);
            z1p = z(N(e,1)) - zgrid;
            z2p = z(N(e,2)) - zgrid;
            z3p = z(N(e,3)) - zgrid;
            z4p = z(N(e,4)) - zgrid;
            a1p = x2p*(y3p*z4p - y4p*z3p) + x3p*(y4p*z2p - y2p*z4p) + x4p*(y2p*z3p - y3p*z2p);
            a2p = x3p*(y4p*z1p - y1p*z4p) + x4p*(y1p*z3p - y3p*z1p) + x1p*(y3p*z4p - y4p*z3p);
            a3p = x4p*(y1p*z2p - y2p*z1p) + x1p*(y2p*z4p - y4p*z2p) + x2p*(y4p*z1p - y1p*z4p);
            a4p = x1p*(y2p*z3p - y3p*z2p) + x2p*(y3p*z1p - y1p*z3p) + x3p*(y1p*z2p - y2p*z1p);

            x1 = nodeIndex(1, N(e,1));  y1 = nodeIndex(2, N(e,1));  z1 = nodeIndex(3, N(e,1));
            x2 = nodeIndex(1, N(e,2));  y2 = nodeIndex(2, N(e,2));  z2 = nodeIndex(3, N(e,2));
            x3 = nodeIndex(1, N(e,3));  y3 = nodeIndex(2, N(e,3));  z3 = nodeIndex(3, N(e,3));
            x4 = nodeIndex(1, N(e,4));  y4 = nodeIndex(2, N(e,4));  z4 = nodeIndex(3, N(e,4));
            V(e) = (1/6)*det([1 1 1 1; x1 x2 x3 x4; y1 y2 y3 y4; z1 z2 z3 z4]);
            if abs(V(e) - (a1p+a2p+a3p+a4p)) < 0.01*V(e)
                a1 = x2*(y3*z4 - y4*z3) + x3*(y4*z2 - y2*z4) + x4*(y2*z3 - y3*z2);
                b1 = (y3*z4 - z3*y4) + (y2*z4 - y4*z2) + (y2*z3 - z2*y3); 
                c1 = (x3*z4 - z3*x4) + (x2*z4 - x4*z2) + (x2*z3 - z2*x3);
                d1 = (x3*y4 - y3*x4) + (x2*y4 - x4*y2) + (x2*y3 - y2*x3); 
                a2 = x3*(y4*z1 - y1*z4) + x4*(y1*z3 - y3*z1) + x1*(y3*z4 - y4*z3);
                b2 = (y3*z4 - z3*y4) + (y1*z4 - z1*y4) + (y1*z3 - z1*y3); 
                c2 = (x3*z4 - z3*x4) + (x1*z4 - z1*x4) + (x1*z3 - z1*x3); 
                d2 = (x3*y4 - y3*x4) + (x1*y4 - y1*x4) + (x1*y3 - y1*x3); 
                a3 = x4*(y1*z2 - y2*z1) + x1*(y2*z4 - y4*z2) + x2*(y4*z1 - y1*z4);
                b3 = (y2*z4 - z2*y4) + (y1*z4 - z1*y4) + (y1*z2 - z1*y2);
                c3 = (x2*z4 - z2*x4) + (x1*z4 - z1*x4) + (x1*z2 - z1*x2);
                d3 = (x2*y4 - y2*x4) + (x1*y4 - y1*x4) + (x1*y2 - y1*x2);
                a4 = x1*(y2*z3 - y3*z2) + x2*(y3*z1 - y1*z3) + x3*(y1*z2 - y2*z1);
                b4 = (y2*z3 - z2*y3) + (y1*z3 - z1*y3) + (y1*z2 - z1*y2);  
                c4 = (x2*z3 - z2*x3) + (x1*z3 - z1*x3) + (x1*z2 - z1*x2); 
                d4 = (x2*y3 - y2*x3) + (x1*y3 - y1*x3) + (x1*y2 - y1*x2);            
                N2 = (b2*(xgrid(i,j) - x1) + c2*(ygrid(i,j) - y1) + d2*(zgrid - z1))/(6*V(e));
                N3 = (b3*(xgrid(i,j) - x1) + c3*(ygrid(i,j) - y1) + d3*(zgrid - z1))/(6*V(e));
                N4 = (b4*(xgrid(i,j) - x1) + c4*(ygrid(i,j) - y1) + d4*(zgrid - z1))/(6*V(e));
                N1 = 1 - N2 - N3 - N4;
                Ezgrid(i,j) = N1*Ez(N(e,1))+N2*Ez(N(e,2))+N3*Ez(N(e,3))+N4*Ez(N(e,4));
            end
        end
    end
end
% Display contour plot of FEM solution
figure;
contourf(xgrid,ygrid,abs(Ezgrid));
xlabel('x (wavelengths)');
ylabel('y (wavelengths)');
title('Total Electric Field Contour Plot');
colorbar;

figure;
contourf(xgrid(10:20,10:20),ygrid(10:20,10:20),abs(Ezgrid(10:20,10:20)));
xlabel('x (wavelengths)');
ylabel('y (wavelengths)');
title('Total Electric Field Contour Plot');
colorbar;

zgridt=[0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25];

for i = 1:length(zgridt)
    EzAtaDistance(Ez(10:20,1),x,y,z,N,nElements,nodeIndex,zgridt(i));
end

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