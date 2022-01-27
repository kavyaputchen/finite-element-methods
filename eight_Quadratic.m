% Constructing the Global Stiffness Matrix

a_e = 1;
c_e = -1;
h_e = 1/8;
k = a_e/(3*h_e)*[7,-8,1;-8,16,-8;1,-8,7]+(c_e*h_e)/30*[4,2,-1;2,16,2;-1,2,4];

K = zeros(17,17);
for i=1:3
    for j=1:3
        K(i,j) = k(i,j);
    end
end
for i=3:5
    for j=3:5
        K(i,j) =K(i,j)+k(i-2,j-2);
    end
end
for i=5:7
    for j=5:7
        K(i,j) =K(i,j)+k(i-4,j-4);
    end
end
for i=7:9
    for j=7:9
        K(i,j) =K(i,j)+k(i-6,j-6);
    end
end
for i=9:11
    for j=9:11
        K(i,j) =K(i,j)+k(i-8,j-8);
    end
end
for i=11:13
    for j=11:13
        K(i,j) =K(i,j)+k(i-10,j-10);
    end
end
for i=13:15
    for j=13:15
        K(i,j) =K(i,j)+k(i-12,j-12);
    end
end
for i=15:17
    for j=15:17
        K(i,j) =K(i,j)+k(i-14,j-14);
    end
end

% COnstructing the F Matrix

h = 1/8;
x = zeros(8,1);
x(1) = 0;
x(2) = h;
for i=3:8
    x(i) = x(i-1)+h;
end
A = zeros(8,1);
B = zeros(8,1);
C = zeros(8,1);
for i=1:8
A(i) = -(h/60)*(-h^2+10*x(i)^2);
B(i) = -(h/15)*(3*h^2+10*x(i)^2+10*x(i)^2*h);
C(i) = -(h/60)*(9*h^2+20*x(i)^2+20*x(i)*h);
end
F = zeros(17,1);
F(1) = A(1);
F(17) = C(8);
F(2) = B(1);
F(16) = B(8);
F(3) = C(1);
F(4) = A(2);
F(5) = B(2);
F(6) = C(2)+A(7);
F(7) = A(3)+C(6);
F(8) = B(3)+B(6);
F(9) = C(3)+A(6);
F(10) = A(4)+C(5);
F(11) = B(4)+B(5);
F(12) = C(4)+A(5);
F(13) = B(7);
F(14) = C(7);
F(15) = A(8);
K_new = zeros(15,15);
for i=2:16
    for j=2:16
        K_new(i-1,j-1)=K(i,j);
    end
end
F_new = zeros(15,1);
for i=2:16
    F_new(i-1)=F(i);
end
 U = (K_new\F_new);
 
 U_new = ones(17,1);
 U_new(1)=0;
 U_new(17)=0;
 for i=2:16
    U_new(i)=U(i-1);
 end
 
 % Exact Solution
 
syms u(x)
Du = diff(u);

ode = diff(u,x,2) == x^2-u;
cond1 = u(0) == 0;
cond2 = u(1) == 0;
conds = [cond1 cond2];
ySol(x) = dsolve(ode,conds);
ySol = simplify(ySol);
x=0:(1/16):1;
Exact_sol= 2.*cos(x)+x.^2-((sin(x)).*(2.*cos(1)-1)/(sin(1)))-2;

% y=[0:(1/16):1];plot(y,U_new);hold on;plot(y,Exact_sol)
 