function Ezgrid = EzGridFunc(Ez,x,y,z,N,nElements,nodeIndex,xgrid,ygrid,zgrid)
for i = 1:1:21
    for j = 1:1:21
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
            if abs(V(e) - (a1p+a2p+a3p+a4p)) < 0.00001*V(e)
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
%axis([-1 1 -1 1]);
title('Total Electric Field Contour Plot');
%axis square;
colorbar;