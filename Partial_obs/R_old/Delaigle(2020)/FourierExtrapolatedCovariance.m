function[Cov,k,lambda1,p]=FourierExtrapolatedCovariance(Y,t,xm,yn)
% PURPOSE: estimate covariance function from fragmented functional data
%--------------------------------------------------------------------------
% USAGE: [Cov,k,lambda1,p]=FourierExtrapolatedCovariance(Y,t,xm,yn)
% where: Y = (linearly interpolated) fragmented functional data 
%               (n x J, missing values are NaN)
%        t = equispaced time grid points corresponding to Y (J x 1)
%                where n = sample size, J = #equispaced time grid points
%       xm = time grids on x-axis for predicting the covariance (n1 x 1)
%       yn = time grids on y-axis for predicting the covairance (n2 x 1)
%--------------------------------------------------------------------------
% RETURNS: Cov = estimated covariance on entire domain (J x J)
%          k = #spline basis for smoothing empirical covariance on S
%          lambda1 = smoothing parameter for spline
%          p = #Fourier basis for extrapolating covariance function
%--------------------------------------------------------------------------
% Note: One can change the candidate values for choosing spline
% paramters in step 0.
%       Subroutines are all included in this file.
% -------------------------------------------------------------------------
% Reference: 
%   Aurore Delaigle, Peter Hall, Wei Huang, and Alois Kneip (2020). 
%       Estimating the covariance of fragmented and other related types of
%       functional data. Journal of American Statistical Association.
% -------------------------------------------------------------------------
% Written by:
%    Wei Huang
%    Lecturer
%    School of Mathematics and Statistics, University of Melbourne
%--------------------------------------------------------------------------
% Last updated:
%    November 19, 2020.
% -------------------------------------------------------------------------
%% Step 0: Set candidate values for choosing smoothing parameters for spline
kp = 1:6;                        %candidate values for #spline basis
lambdap=exp(linspace(-10,0,3));  %candidate values for spline smoothing parameter
%% Step 1: Calculate the empirical covariance on observed domain S
[n,J] = size(Y);

%Find starting and ending index of observed subintervals
idx1 = zeros(n,1);
idx2 = zeros(n,1);
for i =1:n
idx1(i) = find(~isnan(Y(i,:)),1);
idx2(i) = find(~isnan(Y(i,:)),1,'last');
end

%Calculate the #of observation on each grid of the observed domain
count = zeros(J);
  
  for i = 1:n
    count(idx1(i):idx2(i),idx1(i):idx2(i))=count(idx1(i):idx2(i),idx1(i):idx2(i))+1;
  end
  
%Raw empirical covariance (equation (2.2) of the paper)
 Z_P = zeros(J,J);
 for i = 1:J
     for j = i:J
         ind_temp= idx1 <= min(i,j) & idx2 >= max(i,j);
         Y_temp = Y(ind_temp,:);
         Ybar = nanmean(Y_temp,1);
         Y_tmp = Y_temp-Ybar;
         Y_tmp_i = Y_tmp(:,i);
         Y_tmp_j = Y_tmp(:,j);
         Z_P(i,j)= Y_tmp_i'*Y_tmp_j;
         Z_P(j,i) = Z_P(i,j);
     end
 end
 Z_P = Z_P./count;
 Z_P(isnan(Z_P))=0;
 
%Normalise time grid to interval [0,1]
grid = (t-t(1))/(t(end)-t(1));

%remove outliers (determine domian S)
fragment_length=grid(idx2)-grid(idx1);
delta = mean(fragment_length);
 q =  min(10, ceil(max(count(:))/20));
  for i = ceil(J*delta*0.25)+1:J
      for j = 1:ceil(i-J*delta*0.25)
          if count(i,j)<=q
             count(i,j)=nan;
             count(j,i)=nan;
          end
      end
  end
%Calculate empirical convariance on S
Z=Z_P.*count./count;

%% Step 2: Smooth empirical covariance on S using spline
%Vectorise the covariance and delete missing values for the ease of smoothing
int =1/length(t);
[sv,tv] = ndgrid(grid,grid);
sv = sv(:);
tv = tv(:);
z=reshape(Z,size(Z,1)*size(Z,2),1);
ind = isnan(z);
xs = tv(~ind);
ys = sv(~ind);
zs = z(~ind);
%set weight matrix (m(s,t) in equation (B.1))
nw = reshape(count,size(count,1)*size(count,2),1);
nw = nw(~ind);
n=length(xs);
  
%set spline order
order=4;
%set spline exterior knots
extKxl = zeros(1,order);
extKxu = ones(1,order);
extKyl = zeros(1,order);
extKyu = ones(1,order);

warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:SingularMatrix');

%set candidate matrix for choosing #spline basis k and smoothing parameter
[kp, lambdap] = ndgrid(kp, lambdap);
kk = kp(:);
lambdapp = lambdap(:);
pp = [kk'; lambdapp'];

%Choose the smoothing parameters in spline estimator
Cal = arrayfun(@(x)F(pp(1,x),pp(2,x),xs,ys,extKxl,extKxu,extKyl,extKyu,n,order,count,Y,ind,idx1,idx2,grid),1:size(pp,2));
paridx = find(Cal==min(Cal));
k = kk(paridx);
lambda1 = lambdapp(paridx);

%Calculate the spline estimator of covariance K on observed domain S
[N,allknots1,allknots2]= Nform(k,xs,ys,extKxl,extKxu,extKyl,extKyu,n,order);
D=Dform(allknots1,allknots2,n,xs,ys,k,order);
[~,K_hat] = betahat(lambda1,N,D,nw,zs,k,order);


%% Step 3: Extrapolate covariance to S_0 using Fourier basis augmented by x
%set weight for diagonal penalty
nw_diag = ones(J,J);
nw_diag(1:J+1:end) = 1/int;
nw_diag = reshape(nw_diag,J*J,1);
nw_diag = nw_diag(~ind);
diag_ind = tv==sv;
diag_ind = diag_ind(~ind);
  
%Refine indices of time grids for vectorised covariance prediction
 n1 = length(xm);
 n2 = length(yn);
 xms = (xm-xm(1))/(xm(length(xm))-xm(1));
 yns=(yn-yn(1))/(yn(length(yn))-yn(1));
 [xx,yy]=ndgrid(xms,yns);
 xx = xx(:);
 yy = yy(:);
 
%Choosing #Fourier basis p
fsize = floor((n^(1/2)-2)/2); %The largest K that we can choose, where $2*K+2 = p$.
KK= fsize-1:fsize; %Candidate grids for choosing K (in most of time, simply take the largest value)
[C,fi] = arrayfun(@(x)LS(KK(x),xs,ys,K_hat,int,diag_ind,nw_diag),1:length(KK),'UniformOutput',false);
fi = cell2mat(fi);
ind = find(fi==min(fi));
K_slect = KK(ind);
p = K_slect*2+2;
C = C{ind};
Cv = reshape(C,p^2,1);%Estimated coefficients for the Fourier estimator

%Plug in the estimated coeffients to calculate the final estiamtor.
if length(xx)>10000
    num = ceil(length(xx)/10000);
    fhat_vec = zeros(1,length(xx));
    for i = 1:num-1
        fhat_vec((i-1)*10000+1:i*10000) = f_hat(xx((i-1)*10000+1:i*10000),yy((i-1)*10000+1:i*10000),K_slect,Cv);
    end
    fhat_vec((num-1)*10000+1:length(xx))=f_hat(xx((num-1)*10000+1:length(xx)),yy((num-1)*10000+1:length(xx)),K_slect,Cv);
else
    fhat_vec= f_hat(xx,yy,K_slect,Cv);
end

fhat_vec = reshape(fhat_vec,1,length(xx));
f_f=reshape(fhat_vec,n1,n2);
Cov = (f_f+f_f')/2;     %Final symmetric positive semi-definite estimator      
          
end  

%% Subroutines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section 1: Spline smoothing
%Set Spline basis matrix
function [N,allknots1,allknots2]=Nform(K,xs,ys,extKxl,extKxu,extKyl,extKyu,n,order)
  p = linspace(0,1, K+2);
  intKx =  quantile(xs,p(2:K+1));
  intKy =  quantile(ys,p(2:K+1));
  
  allknots1=horzcat(extKxl,intKx,extKxu);
  allknots2=horzcat(extKyl,intKy,extKyu);
  
  k=K+order;
  Bdd1 = fnval(spmak(allknots1,eye(k)),xs)';
  Bdd2 = fnval(spmak(allknots2,eye(k)),ys)';
  
  N=zeros(n,k*k);
  
  for i = 1:n
  
    N_temp=kron(Bdd2(i,:),Bdd1(i,:));
    N(i,:)=N_temp;
  end
end

%Set differentiated spline basis matrix
function D=Dform(allknots1,allknots2,n,xs,ys,K,order)
  k=K+order;
  Bdd1t = fnval(spmak(allknots1,eye(k)),xs)';
  Bdd2t = fnval(spmak(allknots2,eye(k)),ys)';
  Bdd1q = fnval(fnder(spmak(allknots1,eye(k)),2),xs)';
  Bdd2q = fnval(fnder(spmak(allknots2,eye(k)),2),ys)';
  
  D= zeros(n,size(Bdd1q,2)*size(Bdd2t,2));
  for i = 1:n
    D_temp = kron(Bdd2t(i,:),Bdd1q(i,:))+kron(Bdd2q(i,:),Bdd1t(i,:));
    D(i,:) = D_temp;
  end
end

%coefficients estimation for the tensor product penalised spline
function [coef,estimator] = betahat(lambda1,N,D,nw,zs,K,order)
  lc = size(N,2);
  coef = zeros(1,lc);
  
  col.match=arrayfun(@(x)isequal(N(:,x),zeros(1,size(N,1))),1:lc); %delete zero basis
  N=N(:,~col.match);
  D=D(:,~col.match);
  
  WB = bsxfun(@times, N, nw);
  BTB = N'*WB;
  Dl = D'*D;
  Inverse = pinv(BTB+lambda1*Dl);

  for i=1:K+order
    for j=i:K+order
     Inverse((j-1)*(K+order)+i,:)=Inverse((i-1)*(K+order)+j,:);
    end
  end
  
  %Hat=WB*Inverse*WB';
  Cv=Inverse*(WB')*zs;
  estimator=N*Cv;
  coef(~col.match)=Cv;
end

%set leave one curve out cv for spline smoothing
function cv=F(K,lambda1,xs,ys,extKxl,extKxu,extKyl,extKyu,n,order,count,Y,ind,idx1,idx2,t)

  [N,allknots1,allknots2]=Nform(K,xs,ys,extKxl,extKxu,extKyl,extKyu,n,order);
  D=Dform(allknots1,allknots2,n,xs,ys,K,order);
  
  [n,J]=size(Y);
  cv=0;
  for i=1:n
    count_tmp = count;
    count_tmp(idx1(i):idx2(i),idx1(i):idx2(i))=count(idx1(i):idx2(i),idx1(i):idx2(i))-1;
    
    Ymi = Y;
    Ymi(i,:)=[];
    Yi = Y(i,:);
    
    %Calculate leave-one-curve-out empirical convariance and ith raw covariance
    Z_tmp = zeros(J,J);
    Z_emp = zeros(J,J);
    for k = 1:J
     for j = k:J
         ind_temp= idx1 <= min(k,j) & idx2 >= max(k,j);
         ind_temp(i)=[];
         
         Y_temp = Ymi(ind_temp,:);
         Ybar = nanmean(Y_temp,1);
         Y_tmp = Y_temp-Ybar;
         Y_tmp_k = Y_tmp(:,k);
         Y_tmp_j = Y_tmp(:,j);
         Z_tmp(k,j)= Y_tmp_k'*Y_tmp_j;
         Z_tmp(j,k) = Z_tmp(k,j);
         
         Yi_tmp = Yi - Ybar;
         Z_emp(k,j) = Yi_tmp(:,k)*Yi_tmp(:,j);
         Z_emp(j,k) = Z_emp(k,j);
     end
    end
    %q = min(10,ceil(max(count_tmp(:)/20)));
    %count_tmp(count_tmp<=q)=0; %exclude outliers if needed.
    Z_tmp=Z_tmp./count_tmp;  
    Z_emp = Z_emp(idx1(i):idx2(i),idx1(i):idx2(i));
    
    ztmp=reshape(Z_tmp,size(Z_tmp,1)*size(Z_tmp,2),1);
    indt = isnan(ztmp(~ind));
    ztmp=ztmp(~ind);
    ztmp=ztmp(~indt);
    
    Ntmp = N(~indt,:);
    Dtmp = D(~indt,:);
    
    nwtmp = reshape(count_tmp,size(count_tmp,1)*size(count_tmp,2),1);
    nwtmp = nwtmp(~ind);
    nwtmp = nwtmp(~indt);
    
    [coef,~]=betahat(lambda1,Ntmp,Dtmp,nwtmp,ztmp,K,order);
    C = reshape(coef,(K+order),(K+order));
    B  = fnval(spmak(allknots1,eye(K+order)),t);
    est_oco = B'*C*B;
 
    %cross-validation
    er = (Z_emp-est_oco(idx1(i):idx2(i),idx1(i):idx2(i))).^2;
    error = nansum(er(:));
    cv=cv+error;
  end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Section 2: Fourier extrapolation
%Set up Fourier basis (augmented by x) matrix
function N=Bform(K,xs,ys)

  n=length(xs);
  Bdd1 = Gram_Schmidt(xs,@(x)fourierx(x,K,1),0,1);
  Bdd2 = Gram_Schmidt(ys,@(x)fourierx(x,K,1),0,1);
    
  N=zeros(n,size(Bdd1,2)*size(Bdd2,2));
  for i = 1:n
    N_temp=kron(Bdd2(i,:),Bdd1(i,:));
    N(i,:) = N_temp;
  end
end


%Compute unconstraint coefficients as initial values
function coef=fit(WB,BTB,z,p,int)

  Inverse = inv(BTB);
  for i=1:p
    for j=i:p
     Inverse((j-1)*p+i,:)=Inverse((i-1)*p+j,:);
    end
  end

  coef=Inverse*(WB'*z*int^2);
end

%set Sp function for choosing #Fourier basis p
function [C,fval]=LS(K,xs,ys,K_hat,int,diag_ind,nw_diag)  
  p=2*K+2;
  
  N=Bform(K,xs,ys);
  WB = bsxfun(@times, N, nw_diag);
  BTB = N'*WB*int^2;
  inv(BTB);
  [~, msgid] = lastwarn;
  
  if strcmp(msgid,'MATLAB:singularMatrix')
      C=zeros(p,p);
      fval = Inf;
  elseif strcmp(msgid,'MATLAB:nearlySingularMatrix')
      lastwarn('');
      alpha = eps;
      BTB=BTB+alpha*eye(size(BTB,1));
      inv(BTB);
      [~, msgid] = lastwarn;
      if strcmp(msgid,'MATLAB:nearlySingularMatrix')
         while strcmp(msgid,'MATLAB:nearlySingularMatrix')
         lastwarn('');
         alpha = alpha*10;
         BTB=BTB+alpha*eye(size(BTB,1));
         inv(BTB);
         [~, msgid] = lastwarn;
         end
      end
  
  Cv=fit(WB,BTB,K_hat,p,int);
  
  C0 = reshape(Cv,p,p);
 
  fun = @(x)Obj(x,p,N,K_hat,int,diag_ind);

  [Ei_vec,Ei_val] = eig(C0);
  Ei_val(Ei_val<0)=min(Ei_val(Ei_val>0));
  par0 = sqrt(Ei_val)*(Ei_vec)';
  par0 = reshape(par0,p^2,1);
  options= optimoptions(@fminunc,'Algorithm','quasi-Newton','MaxFunctionEvaluations',1000000,'MaxIterations',10000,'SpecifyObjectiveGradient',true,'Display','off');
 
  [par,fval] = fminunc(fun,par0,options);
  V = reshape(par,p,p);
  C = V'*V;

  else
      Cv=fit(WB,BTB,K_hat,p,int);

      C0 = reshape(Cv,p,p);
  fun = @(x)Obj(x,p,N,K_hat,int,diag_ind);
 
 
  [Ei_vec,Ei_val] = eig(C0);
  Ei_val(Ei_val<0)=min(Ei_val(Ei_val>0));
  par0 = sqrt(Ei_val)*(Ei_vec)';
  par0 = reshape(par0,p^2,1);
  options= optimoptions(@fminunc,'Algorithm','quasi-Newton','MaxFunctionEvaluations',1000000,'MaxIterations',10000,'SpecifyObjectiveGradient',true,'Display','off');
 
           [par,fval] = fminunc(fun,par0,options);
            V = reshape(par,p,p); 
            C=V'*V;
  end
end

%Objective function for numerical minimisation of equation (2.5)
function [f,g]=Obj(par,p,N,K_hat,int,diag_ind)
  
    V = reshape(par,p,p);
    C = V'*V;
    coef = reshape(C,p^2,1);
    
    int_S = (K_hat'*K_hat - 2*coef'*(N'*K_hat)+ (coef'*N')*N*coef)*int^2;
    int_D = (K_hat(diag_ind)'*K_hat(diag_ind) - 2*coef'*(N(diag_ind,:)'*K_hat(diag_ind))+ (coef'*N(diag_ind,:)')*N(diag_ind,:)*coef)*int; 
    M_D = (K_hat(diag_ind)'*K_hat(diag_ind) - 2*coef'*(N(diag_ind,:)'*K_hat(diag_ind))+ (coef'*N(diag_ind,:)')*N(diag_ind,:)*coef)*int^2; 
    
    f = int_S-M_D+int_D;
    
   if nargout > 1 % supply gradient
    gfc = (-2*(N'*K_hat)+2*(N'*N)*coef)*int^2 -(-2*(N(diag_ind,:)'*K_hat(diag_ind))+2*(N(diag_ind,:)'*N(diag_ind,:))*coef)*int^2+(-2*(N(diag_ind,:)'*K_hat(diag_ind))+2*(N(diag_ind,:)'*N(diag_ind,:))*coef)*int;
    gcx = zeros(p^2,p^2);
    for i = 1:p
        for j=1:p
            
              for s =1:p  
                    gcx((i-1)*p+s,(j-1)*p+i)=par((j-1)*p+s);
                    gcx((j-1)*p+s,(j-1)*p+i)=par((i-1)*p+s);
                    gcx((i-1)*p+s,(i-1)*p+i)=2*par((i-1)*p+s);
                   
                
             end
        end   
    end
    g = gcx*gfc;
    
    end
end

%Constraint estimator
function fhat=f_hat(x,y,K,Cv)

  Bs1 = Gram_Schmidt(x,@(x)fourierx(x,K,1),0,1);
  Bs2 = Gram_Schmidt(y,@(x)fourierx(x,K,1),0,1);
    
  N=zeros(length(x),size(Bs1,2)*size(Bs2,2));
  
  for i = 1:length(x)
  
    N_temp=kron(Bs2(i,:),Bs1(i,:));
    N(i,:) = N_temp;
  end
  
  fhat = N*Cv;
end

%Fourier basis function (augmented by x)
function f = fourierx(x,M,period)
 scale = sqrt(2/period);
 f = zeros(length(x),2*M+2);
 f(:,1)=1;
 
 for i=1:M
     f(:,2*i)=cos(2*i*pi*x/period)*scale;
     f(:,2*i+1) = sin(2*i*pi*x/period)*scale;
 end

 f(:,2*M+2)=x;
end

%Gram_Schmidt orthogonolisation
function f_n = Gram_Schmidt(x,f,ax,bx)
   
 M = length(f(ax));
 f_n = zeros(length(x),M);
 fx = f(x);
 
 int = 0.001;
 x_seq = linspace(ax,bx,(bx-ax)/int);
 
 f_int = f(x_seq);
 v=f_int(:,1);
 normv=sqrt((v')*v*int);
 f_n(:,1)=fx(:,1)/normv;
 f_n_int = zeros(length(x_seq),M);
 f_n_int(:,1) = f_int(:,1)/normv;
 
 if M>1
 for j=2:M
     proj = (f_int(:,j)')*f_n_int(:,1:j-1).*int;
     proj = repmat(proj,length(x),1);
     projfx = proj.*fx(:,1:j-1);
     f_n(:,j)=fx(:,j)-sum(projfx,2);
     
     proj = (f_int(:,j)')*f_n_int(:,1:j-1).*int;
     proj = repmat(proj,length(x_seq),1);
     projfint = proj.*f_int(:,1:j-1);
     f_n_int(:,j)=f_int(:,j)-sum(projfint,2);
     v= f_n_int(:,j);
     normv = sqrt((v')*v*int);
     
     f_n(:,j)=f_n(:,j)/normv;
     f_n_int(:,j)=f_n_int(:,j)/normv;
     
 end
 end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%