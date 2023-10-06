n_f=93; % number of features
n_n=5000; % number of training samples
d_n=1000; % number of validation samples

for idx = find(y==0)
   y(idx) = -1; % The unit is off when yy2 is -1
end

%% Cross Validation (no. of folds=4)
for ii=1:4
    xx_v(:,1:d_n)=xx_new(:,(ii-1)*d_n+1:ii*d_n);
    yy_v(:,1:d_n)=y(:,(ii-1)*d_n+1:ii*d_n);
    if ii>1
        xx_tr(:,1:(ii-1)*d_n)=xx_new(:,1:(ii-1)*d_n);
        yy_tr(:,1:(ii-1)*d_n)=y(:,1:(ii-1)*d_n);
    end
    xx_tr(:,(ii-1)*d_n+1:n_n-d_n)=xx_new(:,ii*d_n+1:n_n);
    yy_tr(:,(ii-1)*d_n+1:n_n-d_n)=y(:,ii*d_n+1:n_n);

%     Calculating Kernel Function
    KK=zeros(n_n-d_n,n_n-d_n);
    for j=1:n_n-d_n
        for n=1:n_n-d_n
            KK(n,j)=exp(-(xx_tr(:,n)-xx_tr(:,j))'*0.1*eye(length(xx_tr(:,n)))*(xx_tr(:,n)-xx_tr(:,j)));
        end
    end
KK=KK(1:n_n-d_n,1:n_n-d_n);
% Exploring different values for the regularization parameter
radius=[0.000001 0.000005 0.00001 0.00005 0.0001 0.0005 0.001 0.005 0.01 0.05];


for i=1:length(yy_tr(:,1))
    if sum(yy_tr(i,:))==length(yy_tr(1,:))
        alpha=ones((n_n-d_n),1);
    elseif sum(yy_tr(i,:))==-length(yy_tr(1,:))
        alpha=-ones((n_n-d_n),1);
    elseif sum(yy_tr(i,:))==0
        alpha=-ones((n_n-d_n),1);
    else
  % Kenrenalized suppor vector machine with regularization
  for m=1:length(radius)
        alpha_w=sdpvar((n_n-d_n),1,'full');
        vio=sdpvar((n_n-d_n),1,'full');
        lamda=sdpvar(1,1,'full');

        options = sdpsettings('solver','gurobi','gurobi.qcpdual',1);
%         options = sdpsettings('solver','Mosek');
        Constraints=[];  
        Constraints=[Constraints,1-diag(yy_tr(i,:))*KK*alpha_w<=vio];
        Constraints=[Constraints,norm(chol(KK)*alpha_w,2)<=lamda];
        Constraints=[Constraints,vio>=0];
        objective=lamda*radius(m)+(1/(n_n-d_n))*sum(vio);
        sol=optimize(Constraints,objective,options);
        alpha(:,m)=value(alpha_w);
        Objective_SVM{i}=value(objective);
 end
    i
    end
    Alpha{ii,i}=alpha;
end

%      Calculating Kernel Function_v
    KK_v=zeros(n_n-d_n,d_n);
    for j=1:d_n
        for n=1:n_n-d_n
            KK_v(n,j)=exp(-(xx_tr(:,n)-xx_v(:,j))'*0.1*eye(length(xx_tr(:,n)))*(xx_tr(:,n)-xx_v(:,j)));
        end
    end

    for i=1:length(yy_tr(:,1))
          for m=1:length(radius)
          for j=1:d_n
            if sum(yy_tr(i,:))~=length(yy_tr(1,:)) && sum(yy_tr(i,:))~=-length(yy_tr(1,:)) && sum(yy_tr(i,:))~=0

                c1(i,j)=0;
                for n=1:n_n-d_n
                    c1(i,j)=Alpha{ii,i}(n,m)*KK_v(n,j)+c1(i,j);
                end
                yy_t(i,j)=sign(c1(i,j));
            else
                yy_t(i,j)=sign(sum(Alpha{ii,i}(:,1)));
            end
          end
            k(m)=0;
            for j=1:d_n
                if sign(yy_t(i,j))-sign(yy_v(i,j))==0
                    k(m)=k(m)+1;
                end
            end
        radius_best_alpha{ii,i}=radius(find(k==max(k)));
        k_best_alpha{ii,i}=k;
          end
    end

    ii
end




