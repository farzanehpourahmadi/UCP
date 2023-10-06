%load('data_118bus.mat')
tic
n_w=2; %number of wind farms
n_g=54; %number of generators
n_d=91; %number of loads
n_line=186; %number of transmission lines
n_t=24;   %number of operating hours                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  00;
n_b=118;  %number of buses
n_n=5000; %number of samples
mpc=case118; % upload case118 from Matpower
[Ybus, Yf, Yt] = makeYbus(mpc);
PTDF=makePTDF(mpc);
Hg=PTDF*busgen; % sensitivity of each line to each generation unit
Hw=PTDF*buswind; % sensitivity of each line to each wind farm
Hl=PTDF*busload;  % sensitivity of each line to each load
W=[300,300]; % Capacity of two wind farms
%% Solving for each sample
for n=1:n_n
    wp(1,1:n_t)=Windfarm1(1:n_t,1)'.*xx_new(n,1)*W(1);
    wp(2,1:n_t)=Windfarm2(1:n_t,1)'.*xx_new(n,2)*W(2);
    for i=1:n_d
    pd(i,1:n_t)=Demand1(1:n_t,1)'.*xx_new(n,i+2)*Ld(i);
    pq(i,1:n_t)=Demand1(1:n_t,1)'.*xx_new(n,i+2)*Lq(i);
    end  

% definig variables
p=sdpvar(n_g,n_t,'full');
q=sdpvar(n_g,n_t,'full');
c=sdpvar(n_b,n_b,n_t,'full');
s=sdpvar(n_b,n_b,n_t,'full');
ws=sdpvar(n_b,n_t,'full');
ls=sdpvar(n_b,n_t,'full');
wsq=sdpvar(n_b,n_t,'full');
lsq=sdpvar(n_b,n_t,'full');
x=binvar(n_g,n_t,'full');
u=binvar(n_g,n_t,'full');

Constraints=[];
Constraints=[Constraints,ls>=zeros(n_b,n_t)];
Constraints=[Constraints,ws>=zeros(n_b,n_t)];
Constraints=[Constraints,lsq>=zeros(n_b,n_t)];
Constraints=[Constraints,wsq>=zeros(n_b,n_t)];



for t=1:n_t
for i=1:n_b
    Constraints=[Constraints,c(i,i,t)<=Vmax(i)^2];
    Constraints=[Constraints,c(i,i,t)>=Vmin(i)^2];
end
end
for t=1:n_t
for i=1:n_b
    for j=i+1:n_b
        Constraints=[Constraints,c(i,j,t)==c(j,i,t)];
        Constraints=[Constraints,s(i,j,t)==-s(j,i,t)];
% %         Constraints=[Constraints,c(i,j,t)^2+s(i,j,t)^2<=c(i,i,t)*c(j,j,t)];
        Constraints=[Constraints,norm([2*c(i,j,t);2*s(i,j,t);c(i,i,t)-c(j,j,t)],2)<=c(i,i,t)+c(j,j,t)];

    end
end
end
% 
% 
for t=1:n_t
for i=1:n_b
    delta=[line_bus(find(line_bus(:,2)==i),1);line_bus(find(line_bus(:,1)==i),2)];
    if length(delta)>0
    Constraints=[Constraints,busgen(i,:)*p(:,t)+buswind(i,:)*wp(:,t)-busload(i,:)*pd(:,t)-ws(i,t)+ls(i,t)==(-real(Ybus(i,delta))*ones(1,length(delta))'*c(i,i,t)+real(Ybus(i,delta))*c(i,delta,t)'-imag(Ybus(i,delta))*s(i,delta,t)')*100];
    Constraints=[Constraints,busgen(i,:)*q(:,t)-busload(i,:)*qd(:,t)-wsq(i,t)+lsq(i,t)==(imag(Ybus(i,delta))*ones(1,length(delta))'*c(i,i,t)-imag(Ybus(i,delta))*c(i,delta,t)'-real(Ybus(i,delta))*s(i,delta,t)')*100];
    else
    Constraints=[Constraints,busgen(i,:)*p(:,t)+buswind(i,:)*wp(:,t)-busload(i,:)*pd(:,t)-ws(i,t)+ls(i,t)==0];
    Constraints=[Constraints,busgen(i,:)*q(:,t)-busload(i,:)*qd(:,t)-wsq(i,t)+lsq(i,t)==0];   
    end
end
end


% generating units
for t=1:n_t
    for j=1:n_g
        Constraints=[Constraints,p(j,t)<=pgmax(j)*x(j,t)];
        Constraints=[Constraints,p(j,t)>=pgmin(j)*x(j,t)];
        Constraints=[Constraints,q(j,t)<=qgmax(j)*x(j,t)];
        Constraints=[Constraints,q(j,t)>=qgmin(j)*x(j,t)];
    end
end
 for t=2:n_t
     for j=1:n_g
       Constraints=[Constraints,p(j,t)-p(j,t-1)<=ru(j)];
       Constraints=[Constraints,p(j,t-1)-p(j,t)<=rd(j)];
     end
 end

for t=1+1:n_t
    for j=1:n_g
        Constraints=[Constraints,u(j,t)>=x(j,t)-x(j,t-1)];
    end
end
% 
% % min. up/dn time
for j=1:n_g
    for t=lu(j)+1:n_t
        Constraints=[Constraints,u(j,t-lu(j)+1:t)*ones(length(t-lu(j)+1:t),1)<=x(j,t)];
    end
    for t=ld(j)+1:n_t
        Constraints=[Constraints,u(j,t-ld(j)+1:t)*ones(length(t-ld(j)+1:t),1)<=1-x(j,t-ld(j))];
    end
end

% %transmission constraints
for t=1:n_t
for i=1:n_line
    b_f=line_bus(i,1);
    b_t=line_bus(i,2);
    Constraints=[Constraints,(-real(Ybus(b_f,b_t))*c(b_f,b_f,t)+real(Ybus(b_f,b_t))*c(b_f,b_t,t)-imag(Ybus(b_f,b_t))*s(b_f,b_t,t))^2+((imag(Ybus(b_f,b_t)))*c(b_f,b_f,t)-imag(Ybus(b_f,b_t))*c(b_f,b_t,t)-real(Ybus(b_f,b_t))*s(b_f,b_t,t))^2<=fmax(i)^2];   
    Constraints=[Constraints,(-real(Ybus(b_t,b_f))*c(b_t,b_t,t)+real(Ybus(b_t,b_f))*c(b_t,b_f,t)-imag(Ybus(b_t,b_f))*s(b_t,b_f,t))^2+((imag(Ybus(b_t,b_f)))*c(b_t,b_t,t)-imag(Ybus(b_t,b_f))*c(b_t,b_f,t)-real(Ybus(b_t,b_f))*s(b_t,b_f,t))^2<=fmax(i)^2];   

%     Constraints=[Constraints,(-real(Ybus(b_f,b_t))*c(b_f,b_f,t)+real(Ybus(b_f,b_t))*c(b_f,b_t,t)-imag(Ybus(b_f,b_t))*s(b_f,b_t,t))^2+((imag(Ybus(b_f,b_t))-mpc.branch(i,5))*c(b_f,b_f,t)-imag(Ybus(b_f,b_t))*c(b_f,b_t,t)-real(Ybus(b_f,b_t))*s(b_f,b_t,t))^2<=fmax(i)^2];   
%     Constraints=[Constraints,(-real(Ybus(b_t,b_f))*c(b_t,b_t,t)+real(Ybus(b_t,b_f))*c(b_t,b_f,t)-imag(Ybus(b_t,b_f))*s(b_t,b_f,t))^2+((imag(Ybus(b_t,b_f))-mpc.branch(i,5))*c(b_t,b_t,t)-imag(Ybus(b_t,b_f))*c(b_t,b_f,t)-real(Ybus(b_t,b_f))*s(b_t,b_f,t))^2<=fmax(i)^2];   
end
end
% % objective function

Objective=0;

for t=1:n_t
    Objective=Objective+(cost_op(1:n_g)'*p(1:n_g,t)+(cost_st(1:n_g))'*u(1:n_g,t))+10000*ones(1,n_b)*(ls(1:n_b,t)+lsq(1:n_b,t))+10000*ones(1,n_b)*(ws(1:n_b,t)+wsq(1:n_b,t));
end
% Warm-start
assign(x,Active1{n}.x);
assign(u,Active1{n}.u);
assign(x,x1);

ops=sdpsettings('solver','gurobi','usex0',1,'savesolveroutput',1);

sol=optimize(Constraints,Objective,ops);
time_real3(n)=sol.solvertime;

% 
p=value(p);
x=value(x);
u=value(u);

result_ct=Objective;

   n
end

