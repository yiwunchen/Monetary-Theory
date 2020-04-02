beta=0.988;
theta=2; % suppose theta=2
alpha=0.42;
delta=0.025;
rho=0.95;

nc=2 ;  % the number of control variables(c,n)
ns=1 ;  % the number of state variables(k)
ncs=1 ; % the number of co-state variable (lamda)
nex=1 ;  % the number of exogenous variable (A)
nf=2 ;  % the number of additional variables (y and i)

kyratio=alpha/((1/beta)-1+delta) ;
cyratio=1-(kyratio)*(delta) ;
iyratio=1-cyratio ;
ciratio=cyratio*(1/iyratio);
ckratio=(1/kyratio)-delta ;

%calculate the steady-state value for n
ssn=1/(cyratio*(theta/(1-alpha))+1); 

%calculate the steady-state value for k
ssk=ssn/((ckratio+delta)^(1/(1-alpha))) ;

%calculate the steady-state value for y
ssy=(ssk^(alpha))*(ssn^(1-alpha)) ;

mu=alpha*beta*(1/kyratio);

mcc = [-1 0; 0 (ssn+alpha*(1-ssn))/(1-ssn)] ;

mcs = [ 0 1; alpha 1 ] ;

mce = [ 0; 1 ] ;

mss0 = [ (1-alpha)*mu -1; 1 0 ]  ;

mss1 = [ 0 1; -((mu/beta)+1-delta) 0 ] ;

msc0 = [ 0 (1-alpha)*mu; 0 0 ] ;

msc1 = [  0 0; -(ckratio) (1-alpha)*(1/kyratio)] ;

mse0 = [ mu; 0] ;

mse1 = [ 0; (1/kyratio)] ;

fc = [ 0 (1-alpha); -ciratio (1/iyratio)*(1-alpha)];
    
fx = [ alpha 0; (1/iyratio)*alpha 0]; 

fe = [ 1; (1/iyratio)]; 

w = -inv(mss0 - msc0*inv(mcc)*mcs)*(mss1 - msc1*inv(mcc)*mcs);
r =  inv(mss0 - msc0*inv(mcc)*mcs)*(mse0 + msc0*inv(mcc)*mce);
q =  inv(mss0 - msc0*inv(mcc)*mcs)*(mse1 + msc1*inv(mcc)*mce);

[pr,lambr]=eig(w);
alamb=abs(diag(lambr)) ; %the abs value of eigenvalues
[lambs,lambz]=sort(alamb) ; % order the abs value of eigenvalues

lambda=lambr(lambz,lambz) ; %diagnoal matrix of eigenvalues
p=pr(:,lambz) ; %eigenvectors

lamb1=lambda(1:ns,1:ns) ;
lamb2=lambda(ns+1:ns+ncs,ns+1:ns+ncs) ;

p11=p(1:ns,1:ns) ;
p12=p(1:ns,ns+1:ns+ncs) ;
p21=p(ns+1:ns+ncs,1:ns) ;
p22=p(ns+1:ns+ncs,ns+1:ns+ncs) ;

ps=inv(p) ;
ps11=ps(1:ns,1:ns) ;
ps12=ps(1:ns,ns+1:ns+ncs) ;
ps21=ps(ns+1:ns+ncs,1:ns) ;
ps22=ps(ns+1:ns+ncs,ns+1:ns+ncs) ;

rxe=r(1:ns,1:nex) ;
rle=r(ns+1:ns+ncs,1:nex) ;
qxe=q(1:ns,1:nex) ;
qle=q(ns+1:ns+ncs,1:nex) ;

phi0=ps21*rxe+ps22*rle ;
phi1=ps21*qxe+ps22*qle ;

psi=zeros(ncs,nex) ;

for i=1:ncs
 psi(i,:)=-(phi0(i,:)*rho+phi1(i,:))*inv(eye(nex)-rho/lamb2(i,i))/lamb2(i,i) ;
end


xx=p11*lamb1*inv(p11) ; %% x(t)=state variable, z(t):exogenous variable
xe=(p11*lamb1*ps12+p12*lamb2*ps22)*inv(ps22)*psi+rxe*rho+qxe ;
solx=[ xx xe ] ; % express x(t+1) as fun. of s(t)
%s(t) is vector of x(t) and z(t) 

lx=-inv(ps22)*ps21 ;
lex=inv(ps22)*psi ;
soll=[lx lex] ;  %express co-state(t) as fun. of s(t)

cxl=inv(mcc)*mcs ;
ce=inv(mcc)*mce ;
solc=[ cxl*[ eye(ns) ; lx ] cxl*[ zeros(ns,nex) ; lex ]+ce ] ;
% calculate the solutions for control variables

fx(2,:)=solc(1,1:ns)*xx ;
fx=diag(fx);
fe(2,:)=solc(1,1:ns)*xe+solc(1,ns+1:ns+nex)*rho ;

solf=[ fx fe ]+fc*solc ; %the solution for other interested variables(y&i)

m = [ solx ; [ zeros(nex,ns) rho ] ] ;
h = [ soll ; solc ; solf ] ;
mh=[ m ; h ] ;

sigma=0.007;
nimp=50; % variables in response to shock up to nimp periods

drules=mh;
[nex,nex2]=size(sigma) ;
[ny,ns]=size(drules) ;
ns=ns-nex ;
m=drules(1:nex+ns,:) ;
h=drules(nex+ns+1:ny,:) ;
irf=zeros(nimp,6) ;

for k=1:nimp
  if k==1
    irfs=eye(ns+nex) ;
  else
    irfs=m*irfs ;
  end
  irff=h*irfs ;
  irf(k,1)=irfs(ns+1,ns+1) ;
  irf(k,2)=irfs(1,ns+1) ;
  irf(k,3)=irff(2,ns+1) ;
  irf(k,4)=irff(3,ns+1) ;
  irf(k,5)=irff(4,ns+1) ;
  irf(k,6)=irff(5,ns+1) ;
end

t=(1:1:nimp)' ;
zers=zeros(nimp,1) ;

figure(1)

subplot(121)
plot(t,[ irf(:,1) zers],'k' )
title('A')
set(gca,'XLim',[0 nimp+1]) ;

subplot(122)
plot(t, irf(:,4),'b',t,irf(:,5),'k.')
title('N and Y')
set(gca,'XLim',[0 nimp+1]) ;

figure(2)
subplot(121)
plot(t,irf(:,3),'b',t,irf(:,5),'k.')
title('C and Y')
set(gca,'XLim',[0 nimp+1]) ;

subplot(122)
plot(t,irf(:,2),'b',t,irf(:,5),'k.')
title('K and Y')
set(gca,'XLim',[0 nimp+1]) ;

figure(3)
subplot(121)
plot(t,irf(:,6),'b',t,irf(:,5),'k.')
title('I and Y')
set(gca,'XLim',[0 nimp+1]) ;