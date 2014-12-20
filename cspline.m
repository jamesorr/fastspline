function[x,score]=cspline(y,r,T,J,lam)
n=length(y); nc=ceil(n/2); Tcu=T^3; rn=r*n;
x=zeros(1,rn+r-1); c=zeros(1,n); w=zeros(1,n-2); z=zeros(1,n);
for j=1:n-2
   w(j)=y(j)-2*y(j+1)+y(j+2);
end
lamr=lam*Tcu;
a0=6+lamr*2/3; a1=4-lamr/6; a2=sqrt(a1^2-4*(a0-2));
x1=(a1+a2)/2; x2=-(a2-a1)/2;
if lamr>24
   alpha=0.5*(x1+sqrt(x1^2-4)); beta=0.5*(x2+sqrt(x2^2-4));
elseif lamr<24
   alpha=0.5*(x1-sqrt(x1^2-4)); beta=0.5*(x2-sqrt(x2^2-4));
else
   alpha=0.5*(x1-sqrt(x1^2-4)); beta=0.5*(x2+sqrt(x2^2-4));
end

if J>log10(alpha*beta)-(nc-1)*2*log10(abs(alpha))

   %Use untruncated algorithm since N > nc-2
   %Factor coefficient matrix, solve triangular systems, find trace

   e=zeros(1,n-2); f=zeros(1,n-2);
   d=a0; f(1)=1/d; c(2)=f(1)*w(1); mu=a1; e(1)=mu*f(1);
   d=a0-mu*e(1); f(2)=1/d; c(3)=f(2)*(w(2)+mu*c(2));
   mu=a1-e(1); e(2)=mu*f(2);
   for j=3:n-2
       d=a0-mu*e(j-1)-f(j-2); f(j)=1/d;
       c(j+1)=f(j)*(w(j)+mu*c(j)-c(j-1));
       mu=a1-e(j-1); e(j)=mu*f(j);
   end
   c(n-2)=c(n-2)+e(n-3)*c(n-1);
   for j=n-4:-1:1
       c(j+1)=c(j+1)+e(j)*c(j+2)-f(j)*c(j+3);
   end
   g2=f(n-2); tr1=g2; h=e(n-3)*g2; tr2=h;
   g1=f(n-3)+e(n-3)*h; tr1=tr1+g1; tr3=0;
   for k=n-4:-1:n-nc
       q=e(k)*h-f(k)*g2; tr3=tr3+q;
       h=e(k)*g1-f(k)*h; tr2=tr2+h; g2=g1;
       g1=f(k)*(1-q)+e(k)*h; tr1=tr1+g1;
   end
   q=e(n-nc-1)*h-f(n-nc-1)*g2; tr3=tr3+q;
   h=e(n-nc-1)*g1-f(n-nc-1)*h; tr2=tr2+h;

   tr1=6*(2*tr1-rem(n,2)*g1);
   tr2=-8*(2*tr2-(1+rem(n,2))*h);
   tr3=2*(2*tr3-rem(n,2)*q);
   tr=(tr1+tr2+tr3)/n;
else %Use truncated algorithm since N < nc-1
   %Factor coefficient matrix, solve triangular systems, find trace
   flim=alpha*beta; elim=alpha+beta;
   glim=flim*(1+flim)/((1-flim)*((1+flim)^2-elim^2));
   hlim=elim*glim/(1+flim); qlim=elim*hlim-flim*glim;
   N=ceil((log10(flim)-J)/(2*log10(abs(alpha))));
   e=zeros(1,N); f=zeros(1,N);
   d=a0; f(1)=1/d; c(2)=f(1)*w(1); mu=a1; e(1)=mu*f(1);
   d=a0-mu*e(1); f(2)=1/d; c(3)=f(2)*(w(2)+mu*c(2));
   mu=a1-e(1); e(2)=mu*f(2);
   g2=flim; tr1=g2; h=elim*g2; tr2=h;
   g1=flim+elim*h; tr1=tr1+g1; tr3=0;
   for j=3:N
       d=a0-mu*e(j-1)-f(j-2); f(j)=1/d;
       c(j+1)=f(j)*(w(j)+mu*c(j)-c(j-1));
       mu=a1-e(j-1); e(j)=mu*f(j);
       q=elim*h-flim*g2; tr3=tr3+q;
       h=elim*g1-flim*h; tr2=tr2+h; g2=g1;
       g1=flim*(1-q)+elim*h; tr1=tr1+g1;
   end
   tr1=tr1+(nc-N-1)*glim; tr2=tr2+(nc-N)*hlim;
   tr3=tr3+(nc-N)*qlim; tr1=6*(2*tr1-rem(n,2)*glim);
   tr2=-8*(2*tr2-(1+rem(n,2))*hlim); tr3=2*(2*tr3-rem(n,2)*qlim);
   tr=(tr1+tr2+tr3)/n; mu=a1-elim;
   for j=N+1:n-2
       c(j+1)=flim*(w(j)+mu*c(j)-c(j-1));
   end
   c(n-2)=c(n-2)+elim*c(n-1);
   for j=n-3:-1:N+2
       c(j)=c(j)+elim*c(j+1)-flim*c(j+2);
   end
   for j=N:-1:1
       c(j+1)=c(j+1)+e(j)*c(j+2)-f(j)*c(j+3);
   end
end

%Compute GCV score
z(1)=c(2); z(2)=c(3)-2*c(2);
for j=3:n-2
    z(j)=c(j-1)-2*c(j)+c(j+1);
end
z(n-1)=c(n-2)-2*c(n-1); z(n)=c(n-1);
sq=(z*z')/n; score=sq/tr^2;

%Compute estimates
x(r:r:rn)=y-z;
if r<8
   fac1=x(2*r)-x(r)-lamr*c(2)/6; fac2=x(rn)-x(rn-r)+lamr*c(n-1)/6;
   for j=1:r-1
       j1=j/r; j2=1-j1; v=lamr*j1*j2/6; j3=v*(1+j1); j4=v*(2-j1);
       for i=1:n-1
           ir=i*r;
           x(ir+j)=j2*x(ir)+j1*x(ir+r)-j3*c(i+1)-j4*c(i);
       end
       x(j)=x(r)-j2*fac1; x(rn+j)=x(rn)+j1*fac2;
   end
else
   lamc=lamr/(6*r^3); r2=2*r; rsq=r^2;
   for i=1:n-1
       ir=i*r; tmp1=x(ir)/r; tmp2=x(ir+r)/r; tmp3=lamc*c(i+1); tmp4=lamc*c(i);
       for j=1:r-1
           x(ir+j)=(r-j)*tmp1+j*tmp2-j*(rsq-j*j)*tmp3-j*(r-j)*(r2-j)*tmp4;
       end
   end
   tmp1=x(r)-x(r+1); tmp2=x(rn)-x(rn-1);
   for j=1:r-1
       x(j)=x(r)+(r-j)*tmp1; x(rn+j)=x(rn)+j*tmp2;
   end
end
