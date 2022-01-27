% Constructing the Global Stiffness Matrix
a_e = 1;
c_e = -1;
h_e = 1/8;
k = a_e/h_e*[1,-1;-1,1]+(c_e*h_e/6)*[2,1;1,2];
K = zeros(9,9);

for i=1:2
    for j=1:2
    K(i,j)=k(i,j);
    end
end

for i=2:3
    for j=2:3
    K(i,j)=K(i,j)+k(i-1,j-1);
    end
end
for i=3:4
    for j=3:4
    K(i,j)=K(i,j)+k(i-2,j-2);
    end
end
for i=4:5
    for j=4:5
    K(i,j)=K(i,j)+k(i-3,j-3);
    end
end
for i=5:6
    for j=5:6
    K(i,j)=K(i,j)+k(i-4,j-4);
    end
end
for i=6:7
    for j=6:7
    K(i,j)=K(i,j)+k(i-5,j-5);
    end
end
for i=7:8
    for j=7:8
    K(i,j)=K(i,j)+k(i-6,j-6);
    end
end
for i=8:9
    for j=8:9
    K(i,j)=K(i,j)+k(i-7,j-7);
    end
end

% Constructing the F Matrix
h = 1/8;
x = zeros(9,1);
x(1)=0;
x(2)=h;
for i=3:9
   x(i) = x(i-1)+h;
end
A = zeros(8,1);
B = zeros(8,1);
for i=1:8
A(i,1) = vpa(-1/h*((x(i+1)/3)*(x(i+1)^3-x(i)^3)-0.25*(x(i+1)^4-x(i)^4)));
B(i,1) = vpa(-1/h*(-(x(i)/3)*(x(i+1)^3-x(i)^3)+0.25*(x(i+1)^4-x(i)^4)));
end
F = zeros(9,1);
F(1) = A(1);
F(9) = B(8);
for i=2:8
    F(i) = B(i-1)+A(i);
end

K_new = zeros(7,7);
for i=2:8
    for j=2:8
        K_new(i-1,j-1)=K(i,j);
    end
end
F_new = zeros(7,1);
for i=2:8
    F_new(i-1)=F(i);
end
 U = (K_new\F_new);
 
 U_new = ones(9,1);
 U_new(1)=0;
 U_new(9)=0;
 for i=2:8
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
x=0:(1/8):1;
Exact_sol= 2.*cos(x)+x.^2-((sin(x)).*(2.*cos(1)-1)/(sin(1)))-2;

% y=[0:(1/8):1];plot(y,U_new);hold on;plot(y,Exact_sol)
