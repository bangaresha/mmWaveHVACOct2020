function Ez_eval = EzAtaDistance(Ez,x,y,z,N,nElements,nodeIndex,zeval)
d2p=pi/180;
dista=1; % where to evaluate the field
Ez_eval=zeros(1,1441);
phi=zeros(1,1441);
for I=1:1441 %721
    phi(I)=(I-1)*0.25; %0.5;
    xeval=dista*cos(phi(I)*d2p);
    yeval=dista*sin(phi(I)*d2p);
   % zeval=dist;
    for e=11:20 %1:nElements
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
            N2 = (b2*(xeval - x1) + c2*(yeval - y1) + d2*(zeval - z1))/(6*V(e));
            N3 = (b3*(xeval - x1) + c3*(yeval - y1) + d3*(zeval - z1))/(6*V(e));
            N4 = (b4*(xeval - x1) + c4*(yeval - y1) + d4*(zeval - z1))/(6*V(e));
            N1 = 1 - N2 - N3 - N4;
            Ez_eval = N1*Ez(N(e,1))+N2*Ez(N(e,2))+N3*Ez(N(e,3))+N4*Ez(N(e,4));
        end
    end
end
% Plot the analytical solution and the FEM solution at a distance 'dist'
figure;
plot(phi,abs(Ez_eval),'b--'),legend('FEM (1^{st} order ABC)');
xlabel('Angle (degrees)');
ylabel('Electric Field (V/m)');
grid on;
end