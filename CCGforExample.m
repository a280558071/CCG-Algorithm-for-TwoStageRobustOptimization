%% CCG algorithm for Two-stage Robust Optimization
% Original Problem: min_{y} c'*y + max_{u} min_{x∈F(y,u)} b'*x s.t. A*y>=d, y∈S_y
% where F(y,u)={x∈S_x: G*x >= h - E*y - M*u}
% SP2: max_{u in U} min_{x} {b'*x: G*x>=h-E*s_y-Mu, G'*pi<=b, u∈U, pi>=0} is
% solved by its MIP recouse derived by KKT conditions.
% [1] Zeng, Bo, and Long Zhao. "Solving two-stage robust optimization problems using a column-and-constraint generation method." Operations Research Letters 41, no. 5 (2013): 457-461.

clear all
%% Constant and Parameter matrix setting
f = [400, 414, 326]';
a = [18, 25, 20]';
C = [22, 33, 24;
    33, 23, 30;
    20, 25, 27];
b=[22, 33, 24, 33, 23, 30, 20, 25, 27]';
dl = [206, 274, 220]';
du = [40, 40, 40]';
G=zeros(6,9);
G(1,1:3)=-ones(1,3);
G(2,4:6)=-ones(1,3);
G(3,7:9)=-ones(1,3);
G(4:6,1:3)=diag(ones(3,1));
G(4:6,4:6)=diag(ones(3,1));
G(4:6,7:9)=diag(ones(3,1));
E=zeros(6,6);
E(1:3,4:6)=diag(ones(3,1));
M=zeros(6,3);
M(4:6,:)=-diag(du.*ones(3,1));
h=zeros(6,1);
h(4:6)=dl;
BigM=1e5;

%% Variable statement
MaxIter=100; % Max iteration
% for i=1:MaxIter
%    x{i}=sdpvar(9,MaxIter,'full');
% end
x=sdpvar(9,MaxIter,'full');
y=binvar(3,1);
z=sdpvar(3,1);
pi=sdpvar(6,1);
d=sdpvar(3,1);
g=sdpvar(3,1);
u=binvar(9,1); % ancillary binary variables for BigM on (b-G'*pi).*x==0, u(j)==0 means x(j)==0, (b-G'*pi)(j) free
v=binvar(6,1); % ancillary binary variables for BigM on (G*x-(h-E*y-M*g)).*pi==0, v(i)==0 means pi(i)==0, (G*x-(h-E*y-M*g))(i) free

%% CCG algorithm
% Original Problem: min_{y} c'*y + max_{u} min_{x∈F(y,u)} b'*x s.t. A*y>=d,
% y∈S_y
% where F(y,u)={x∈S_x: G*x >= h - E*y - M*u}
% SP2: max_{u,pi} {(h-E*y-Mu)'*pi:G'*pi<=b,u∈U,pi>=0}
LB=-inf;
UB=inf;
k=1;
eta=sdpvar(1);
Obj_MP2=[f;a]'*[y;z]+eta;
Cons_MP2=[z<=800*y, x(:,1)>=0, z>=0, eta>=0, sum(z)>=772];
% sum(z)>=772 is added because we need to ensure Subproblem is feasible
% when s_y and s_z is produced in first-round of MP2.
% Cons_MP2=[Cons_MP2, eta>=b'*x(:,1)];
ops=sdpsettings('solver','gurobi','verbose',0);

Epsilon=0.01;
while UB-LB>=Epsilon
    % Solve MP2
    sol_MP2=optimize(Cons_MP2,Obj_MP2,ops);
    s_y=value(y);
    s_z=value(z);
    s_eta=value(eta);
    LB=value(Obj_MP2);
    
    % Solve SP2
    Obj_SP2 = -b'*x(:,k);
    Cons_SP2 = [pi>=0, x(:,k)>=0, G'*pi<=b, 1>=g>=0, sum(g)<=1.8, g(1)+g(2)<=1.2, d==dl+du.*g];
    Cons_SP2 = [Cons_SP2, G*x(:,k) >= h-E*[s_y;s_z]-M*g];
    Cons_SP2 = [Cons_SP2, G*x(:,k)-(h-E*[s_y;s_z]-M*g) <= BigM*(1-v), pi<=BigM*v];
    Cons_SP2 = [Cons_SP2, b-G'*pi <= BigM*(1-u), x(:,k)<=BigM*u];
    sol_SP2 = optimize(Cons_SP2,Obj_SP2,ops);
    s_g = value(g);
    
    % Add constraints and variables "x(:,k+1)" in MP2
    if sol_SP2.problem==0 % SP2 is solved
        UB=min(UB,[f;a]'*[s_y;s_z]+value(-Obj_SP2));
        display(['Iter ',num2str(k),' g = ',num2str(s_g')]);
        Cons_MP2 = [Cons_MP2, eta>=b'*x(:,k+1), E*[y;z]+G*x(:,k+1)>=h-M*s_g, x(:,k+1)>=0];
        % !!! "x(:,k+1)>=0" is important!!! Dont't forget it !!!
    else % SP2 is unbounded, not completed yet. Because still don't know how to identify scenario for which Q(y*)=inf., i.e. s_g could not be produced.
        Cons_MP2 = [Cons_MP2, eta>=b'*x(:,k+1), E*[y;z]+G*x(:,k+1)>=h-M*s_g, x(:,k+1)>=0];
    end
    
    k=k+1;
    display(['UB: ',num2str(UB),' LB: ',num2str(LB)]);
    if k>=100
        break
    end
end