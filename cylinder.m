clear;
close all;


%geometryFromEdges(model,@circleg);

numberOfPDE = 1;
model = createpde;
gm = geometryFromEdges(model,@circleg);
%gm = multicylinder([0.5 0.6],20,'Void',[true,false]);
model.Geometry = gm;
pdegplot(model,'CellLabels','on','FaceAlpha',0.5);
title 'Geometry With Edge Labels Displayed';
xlabel x
ylabel y

% Create and view a finite element mesh for the problem.
mes = generateMesh(model);
figure
pdemesh(model);
xlabel x
ylabel y

% Specify PDE Coefficients

epsilon = 8.85E-12;
mu = 4*pi*1E-7;
c = 1/mu;  %new1
%initial c = 1;
%initial a = 0;
%f = 2.5E9;
f = 50;
omega = 2*pi*f;
%sigma = 1E4;
sigma = 57E6;
a = (omega^2)*epsilon - 1i*omega*sigma;   %new1
f = 10; %new1
%initial f = 0;
%initial m = 1;
%initial d = 0;
m = 0;  %new1
d = 0;  %new1
specifyCoefficients(model,'m',m,'d',d,'c',c,'a',a,'f',f);


% Set zero Dirichlet boundary conditions on the left (edge 4) and right (edge
% 2) and zero Neumann boundary conditions on the top (edge 1) and bottom
% (edge 3).
applyBoundaryCondition(model,'dirichlet','Face',1,'u',0);
%applyBoundaryCondition(model,'neumann','Edge',([1 3]),'g',0);

u0 = 10;
% u0 = @(location) atan(cos(pi/2*location.x));
%ut0 = @(location) 3*sin(pi*location.x).*exp(sin(pi/2*location.y));
 setInitialConditions(model,u0);

% Find the solution at 31 equally-spaced points in time from 0 to 5.
% n = 31;
%n=2;
% tlist = linspace(0,5,n);
%tlist = linspace(0,1,n);
model.SolverOptions.ReportStatistics ='on';
%result = solvepde(model,tlist);
result = solvepde(model);
u = result.NodalSolution;

figure
pdeplot3D(model,'ColorMapData',real(u(:,1)))
title('u(1)')

nodZ = result.Mesh.Nodes(3,:);
nodX = result.Mesh.Nodes(1,:);
nodY = result.Mesh.Nodes(2,:);
nodUnZ = unique(result.Mesh.Nodes(3,:));

for i = 1:length(nodZ)
    for j = 1:length(nodUnZ)
        if nodZ(i) == nodUnZ(j)
            tab(i,:) = [nodX(i); nodY(i); nodZ(i); u(i)];
        end
    end
end

%zlen = [1 4 9 14 19 24 29 34 39 44 49];
zlenL = [1 9 19];
[tab1, count1] = functab(tab, nodX, nodY, nodZ, u, 1);
[tab2, count2] = functab(tab, nodX, nodY, nodZ, u, 9);
[tab3, count3] = functab(tab, nodX, nodY, nodZ, u, 19);

function [tab1, count1]= functab(tab, nodX, nodY, nodZ, u, zlenL)
zlenU = zlenL + 1.1;
for i = 1:length(tab(:,3))
    if tab(i,3) > zlenL && tab(i,3) < zlenU
        tab1(i,:) = [nodX(i); nodY(i); nodZ(i); u(i)];
    end
end
figure
plot(tab1(:,1),tab1(:,4))
title('nodal solution at Z=1 versus X');
count1 = 0;
for i=1:length(tab1(:,3))
    if tab1(i,4) > 0.7
        count1 = count1+1;
    end
end
end


% for i=1:length(tab(:,3))
%     if tab(i,3) > 1 && tab(i,3) < 2
%         tab1(i,:) = [nodX(i); nodY(i); nodZ(i); u(i)];
%     end
%     if tab(i,3) > 9 && tab(i,3) < 10.1
%         tab2(i,:) = [nodX(i); nodY(i); nodZ(i); u(i)];
%     end
%     if tab(i,3) > 19 && tab(i,3) < 20.1
%         tab3(i,:) = [nodX(i); nodY(i); nodZ(i); u(i)];
%     end
%     if tab(i,3) > 29 && tab(i,3) < 30.1
%         tab4(i,:) = [nodX(i); nodY(i); nodZ(i); u(i)];
%     end
%     if tab(i,3) > 39 && tab(i,3) < 40.1
%         tab5(i,:) = [nodX(i); nodY(i); nodZ(i); u(i)];
%     end
%     if tab(i,3) > 49 && tab(i,3) < 50.1
%         tab6(i,:) = [nodX(i); nodY(i); nodZ(i); u(i)];
%     end
%     if tab(i,3) > 59 && tab(i,3) < 60.1
%         tab7(i,:) = [nodX(i); nodY(i); nodZ(i); u(i)];
%     end
%     if tab(i,3) > 69 && tab(i,3) < 70.1
%         tab8(i,:) = [nodX(i); nodY(i); nodZ(i); u(i)];
%     end
%     if tab(i,3) > 79 && tab(i,3) < 80.1
%         tab9(i,:) = [nodX(i); nodY(i); nodZ(i); u(i)];
%     end
% end
% 
% figure
% plot(tab1(:,1),tab1(:,4))
% title('nodal solution at Z=1 versus X');
% 
% figure
% plot(tab1(:,2),tab1(:,4))
% title('nodal solution at Z=1 versus Y');
% 
% figure
% plot(tab2(:,1),tab2(:,4))
% title('nodal solution at Z=10 versus X');
% 
% figure
% plot(tab2(:,2),tab2(:,4))
% title('nodal solution at Z=10 versus Y');
% 
% figure
% plot(tab3(:,1),tab3(:,4))
% title('nodal solution at Z=20 versus X');
% 
% figure
% plot(tab3(:,2),tab3(:,4))
% title('nodal solution at Z=20 versus Y');
% 
% figure
% plot(tab4(:,1),tab4(:,4))
% title('nodal solution at Z=30 versus X');
% 
% figure
% plot(tab5(:,1),tab5(:,4))
% title('nodal solution at Z=40 versus X');
% 
% figure
% plot(tab6(:,1),tab6(:,4))
% title('nodal solution at Z=50 versus X');
% 
% figure
% plot(tab7(:,2),tab7(:,4))
% title('nodal solution at Z=60 versus X');
% 
% figure
% plot(tab8(:,2),tab8(:,4))
% title('nodal solution at Z=70 versus X');
% 
% figure
% plot(tab9(:,2),tab9(:,4))
% title('nodal solution at Z=80 versus X');
% 
% count1 =0;  count2 =0;  count3 =0;  count4 =0; count5 = 0; count6 = 0;
% count7 =0;  count8 =0;  count9 =0;  count10 =0; count11 = 0;
% 
% for i=1:length(tab1(:,3))
%     if tab1(i,4) > 0.7
%         count1 = count1+1;
%        tab11(i,:) = [nodX(i); nodY(i); nodZ(i); u(i)];
%     end
% end
% for i=1:length(tab2(:,3))
%     if tab2(i,4) > 0.7
%         count2 = count2+1;
%        tab21(i,:) = [nodX(i); nodY(i); nodZ(i); u(i)];
%     end
% end
% for i=1:length(tab3(:,3))
%     if tab3(i,4) > 0.7
%         count3 = count3+1;
%        tab31(i,:) = [nodX(i); nodY(i); nodZ(i); u(i)];
%     end
% end
% for i=1:length(tab4(:,3))
%     if tab4(i,4) > 0.7
%         count4 = count4+1;
%        tab41(i,:) = [nodX(i); nodY(i); nodZ(i); u(i)];
%     end
% end
% for i=1:length(tab5(:,3))
%     if tab5(i,4) > 0.7
%         count5 = count5+1;
%        tab51(i,:) = [nodX(i); nodY(i); nodZ(i); u(i)];
%     end
% end
% for i=1:length(tab6(:,3))
%     if tab6(i,4) > 0.7
%         count6 = count6+1;
%        tab61(i,:) = [nodX(i); nodY(i); nodZ(i); u(i)];
%     end
% end
% for i=1:length(tab7(:,3))
%     if tab7(i,4) > 0.7
%         count7 = count7+1;
%        tab71(i,:) = [nodX(i); nodY(i); nodZ(i); u(i)];
%     end
% end
% for i=1:length(tab8(:,3))
%     if tab8(i,4) > 0.7
%         count8 = count8+1;
%        tab81(i,:) = [nodX(i); nodY(i); nodZ(i); u(i)];
%     end
% end
% for i=1:length(tab9(:,3))
%     if tab9(i,4) > 0.7
%         count9 = count9+1;
%        tab91(i,:) = [nodX(i); nodY(i); nodZ(i); u(i)];
%     end
% end

