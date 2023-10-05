   
%load('parameters_118bus.mat')
tic
n_w=2; %number of m plants
n_g=54; %number of generators
n_d=91; %number of loads
n_line=186; %number of transmission lines
n_t=24;   %number of operating hours                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  00;
n_b=118;  %number of buses
d_n=5000; %number of samples
mpc=case118;
[Ybus, Yf, Yt] = makeYbus(mpc);
PTDF=makePTDF(mpc);
Hg=PTDF*busgen; % sensitivity of each line to each generation unit
Hw=PTDF*buswind; % sensitivity of each line to each wind farm
Hl=PTDF*busload;  % sensitivity of each line to each load
n_n=5000;
W=[300,300];
%% Solving for each sample
for n=1:n_d
    wp(1,1:n_t)=Windfarm1(1:n_t,n)'*xx_new(n,1)*W(1);
    wp(2,1:n_t)=Windfarm2(1:n_t,n)'*xx_new(n,2)*W(2);
    for i=1:n_d
    pd(i,1:n_t)=Demand1(1:n_t,1)'*xx_new(n,i+2)*Ld(i);
    pq(i,1:n_t)=Demand1(1:n_t,1)'*xx_new(n,i+2)*Lq(i);
    end


% definig variables
p=sdpvar(n_g,n_t,'full');
ws=sdpvar(n_t,1,'full');
ls=sdpvar(n_t,1,'full');
x=binvar(n_g,n_t,'full');
u=binvar(n_g,n_t,'full');
tic

Constraints=[];
Constraints=[Constraints,ls>=zeros(n_t,1)];
Constraints=[Constraints,ws>=zeros(n_t,1)];
 for t=1:n_t
    Constraints=[Constraints,ones(1,n_g)*p(:,t)+ones(1,n_w)*(wp(:,t))-ws(t)+ls(t)==ones(1,n_d)*(pd(:,t))];
 end
%limitations of generating units
for t=1:n_t
    for j=1:n_g
        Constraints=[Constraints,p(j,t)<=pgmax(j)*x(j,t)];
        Constraints=[Constraints,p(j,t)>=pgmin(j)*x(j,t)];         
    end
end

 for t=2:n_t
     for j=1:n_g
     Constraints=[Constraints,p(j,t)-p(j,t-1)<=x(j,t-1)*ru(j)+(1-x(j,t-1))*QS(j)];
     Constraints=[Constraints,p(j,t-1)-p(j,t)<=x(j,t)*rd(j)+(1-x(j,t))*QS(j)];

     end
 end

for t=1+1:n_t
    for j=1:n_g
        Constraints=[Constraints,u(j,t)>=x(j,t)-x(j,t-1)];
    end
end

% min. up/dn time
for j=1:n_g
    for t=lu(j)+1:n_t
        Constraints=[Constraints,u(j,t-lu(j)+1:t)*ones(length(t-lu(j)+1:t),1)<=x(j,t)];
    end
    for t=ld(j)+1:n_t
        Constraints=[Constraints,u(j,t-ld(j)+1:t)*ones(length(t-ld(j)+1:t),1)<=1-x(j,t-ld(j))];
    end
end

%transmission constraints
for t=1:n_t
    for j=1:n_line
        Constraints=[Constraints,Hg(j,1:n_g)*p(:,t)+Hw(j,:)*(wp(:,t))-Hl(j,:)*(pd(:,t))<=fmax(j)];
        Constraints=[Constraints,Hg(j,1:n_g)*p(:,t)+Hw(j,:)*(wp(:,t))-Hl(j,:)*(pd(:,t))>=-fmax(j)];
    end
end


Objective=0;

for t=1:n_t
    Objective=Objective+(cost_op(1:n_g)'*p(1:n_g,t)+(cost_st(1:n_g))'*u(1:n_g,t))+1000*ls(t)+1000*ws(t);
end
% assign(x,Active1_Test{ii}{n}.x);
% % assign(u,Active1_Test{ii}{n}.u);
% assign(p,Active1_Test{ii}{n}.p);

% ops = sdpsettings('solver','gurobi','verbose',0,'savesolveroutput',1);
ops=sdpsettings('solver','gurobi','usex0',1,'savesolveroutput',1);


sol=optimize(Constraints,Objective,ops);
time=toc
time_real3(n)=sol.solvertime;
% 
p=value(p);
x=value(x);
u=value(u);

result_ct=Objective;

%% Specifying active constraints

Active{n}.x=x;
Active{n}.u=u;

end

