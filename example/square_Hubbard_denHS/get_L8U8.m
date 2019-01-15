beta=1:10;
m=0;
cd 1dL8beta1U8/
load jxf_k.dat; load jxh_k.dat; Nt=length(jxh_k)/3; tau=jxh_k(1:Nt,1);
y=jxh_k(1:Nt,5)*4+jxf_k(1:Nt,5)*2;
dy=sqrt(jxh_k(1:Nt,7).^2*16+jxf_k(1:Nt,7).^2*4);
n=1;
s(1)=y(n);ds(1)=dy(n);
dtau=tau(2)-tau(1);fy(1)=sum(y(1:Nt-1).*exp(2i*m*pi*tau(1:Nt-1)/1))*dtau;
cd ../1dL8beta2U8/
load jxf_k.dat; load jxh_k.dat; Nt=length(jxh_k)/3; tau=jxh_k(1:Nt,1);
y=jxh_k(1:Nt,5)*4+jxf_k(1:Nt,5)*2;
dy=sqrt(jxh_k(1:Nt,7).^2*16+jxf_k(1:Nt,7).^2*4);
n=1;
s(2)=y(n);ds(2)=dy(n);
dtau=tau(2)-tau(1);fy(2)=sum(y(1:Nt-1).*exp(2i*m*pi*tau(1:Nt-1)/2))*dtau;
cd ../1dL8beta3U8/
load jxf_k.dat; load jxh_k.dat; Nt=length(jxh_k)/3; tau=jxh_k(1:Nt,1);
y=jxh_k(1:Nt,5)*4+jxf_k(1:Nt,5)*2;
dy=sqrt(jxh_k(1:Nt,7).^2*16+jxf_k(1:Nt,7).^2*4);
n=1;
s(3)=y(n);ds(3)=dy(n);
dtau=tau(2)-tau(1);fy(3)=sum(y(1:Nt-1).*exp(2i*m*pi*tau(1:Nt-1)/3))*dtau;
cd ../1dL8beta4U8/
load jxf_k.dat; load jxh_k.dat; Nt=length(jxh_k)/3; tau=jxh_k(1:Nt,1);
y=jxh_k(1:Nt,5)*4+jxf_k(1:Nt,5)*2;
dy=sqrt(jxh_k(1:Nt,7).^2*16+jxf_k(1:Nt,7).^2*4);
n=1;
s(4)=y(n);ds(4)=dy(n);
dtau=tau(2)-tau(1);fy(4)=sum(y(1:Nt-1).*exp(2i*m*pi*tau(1:Nt-1)/4))*dtau;
cd ../1dL8beta5U8/
load jxf_k.dat; load jxh_k.dat; Nt=length(jxh_k)/3; tau=jxh_k(1:Nt,1);
y=jxh_k(1:Nt,5)*4+jxf_k(1:Nt,5)*2;
dy=sqrt(jxh_k(1:Nt,7).^2*16+jxf_k(1:Nt,7).^2*4);
n=1;
s(5)=y(n);ds(5)=dy(n);
dtau=tau(2)-tau(1);fy(5)=sum(y(1:Nt-1).*exp(2i*m*pi*tau(1:Nt-1)/5))*dtau;
cd ../1dL8beta6U8/
load jxf_k.dat; load jxh_k.dat; Nt=length(jxh_k)/3; tau=jxh_k(1:Nt,1);
y=jxh_k(1:Nt,5)*4+jxf_k(1:Nt,5)*2;
dy=sqrt(jxh_k(1:Nt,7).^2*16+jxf_k(1:Nt,7).^2*4);
n=1;
s(6)=y(n);ds(6)=dy(n);
dtau=tau(2)-tau(1);fy(6)=sum(y(1:Nt-1).*exp(2i*m*pi*tau(1:Nt-1)/6))*dtau;
cd ../1dL8beta7U8/
load jxf_k.dat; load jxh_k.dat; Nt=length(jxh_k)/3; tau=jxh_k(1:Nt,1);
y=jxh_k(1:Nt,5)*4+jxf_k(1:Nt,5)*2;
dy=sqrt(jxh_k(1:Nt,7).^2*16+jxf_k(1:Nt,7).^2*4);
n=1;
s(7)=y(n);ds(7)=dy(n);
dtau=tau(2)-tau(1);fy(7)=sum(y(1:Nt-1).*exp(2i*m*pi*tau(1:Nt-1)/7))*dtau;
cd ../1dL8beta8U8/
load jxf_k.dat; load jxh_k.dat; Nt=length(jxh_k)/3; tau=jxh_k(1:Nt,1);
y=jxh_k(1:Nt,5)*4+jxf_k(1:Nt,5)*2;
dy=sqrt(jxh_k(1:Nt,7).^2*16+jxf_k(1:Nt,7).^2*4);
n=1;
s(8)=y(n);ds(8)=dy(n);
dtau=tau(2)-tau(1);fy(8)=sum(y(1:Nt-1).*exp(2i*m*pi*tau(1:Nt-1)/8))*dtau;
cd ../1dL8beta9U8/
load jxf_k.dat; load jxh_k.dat; Nt=length(jxh_k)/3; tau=jxh_k(1:Nt,1);
y=jxh_k(1:Nt,5)*4+jxf_k(1:Nt,5)*2;
dy=sqrt(jxh_k(1:Nt,7).^2*16+jxf_k(1:Nt,7).^2*4);
n=1;
s(9)=y(n);ds(9)=dy(n);
dtau=tau(2)-tau(1);fy(9)=sum(y(1:Nt-1).*exp(2i*m*pi*tau(1:Nt-1)/9))*dtau;
cd ../1dL8beta10U8/
load jxf_k.dat; load jxh_k.dat; Nt=length(jxh_k)/3; tau=jxh_k(1:Nt,1);
y=jxh_k(1:Nt,5)*4+jxf_k(1:Nt,5)*2;
dy=sqrt(jxh_k(1:Nt,7).^2*16+jxf_k(1:Nt,7).^2*4);
n=1;
s(10)=y(n);ds(10)=dy(n);
dtau=tau(2)-tau(1);fy(10)=sum(y(1:Nt-1).*exp(2i*m*pi*tau(1:Nt-1)/10))*dtau;
cd ..
fy0=fy;
%-----------------------------------------------------------------------------
m=1;
cd 1dL8beta1U8/
load jxf_k.dat; load jxh_k.dat; Nt=length(jxh_k)/3; tau=jxh_k(1:Nt,1);
y=jxh_k(1:Nt,5)*4+jxf_k(1:Nt,5)*2;
dy=sqrt(jxh_k(1:Nt,7).^2*16+jxf_k(1:Nt,7).^2*4);
dtau=tau(2)-tau(1);fy(1)=sum(y(1:Nt-1).*exp(2i*m*pi*tau(1:Nt-1)/1))*dtau;
cd ../1dL8beta2U8/
load jxf_k.dat; load jxh_k.dat; Nt=length(jxh_k)/3; tau=jxh_k(1:Nt,1);
y=jxh_k(1:Nt,5)*4+jxf_k(1:Nt,5)*2;
dy=sqrt(jxh_k(1:Nt,7).^2*16+jxf_k(1:Nt,7).^2*4);
dtau=tau(2)-tau(1);fy(2)=sum(y(1:Nt-1).*exp(2i*m*pi*tau(1:Nt-1)/2))*dtau;
cd ../1dL8beta3U8/
load jxf_k.dat; load jxh_k.dat; Nt=length(jxh_k)/3; tau=jxh_k(1:Nt,1);
y=jxh_k(1:Nt,5)*4+jxf_k(1:Nt,5)*2;
dy=sqrt(jxh_k(1:Nt,7).^2*16+jxf_k(1:Nt,7).^2*4);
dtau=tau(2)-tau(1);fy(3)=sum(y(1:Nt-1).*exp(2i*m*pi*tau(1:Nt-1)/3))*dtau;
cd ../1dL8beta4U8/
load jxf_k.dat; load jxh_k.dat; Nt=length(jxh_k)/3; tau=jxh_k(1:Nt,1);
y=jxh_k(1:Nt,5)*4+jxf_k(1:Nt,5)*2;
dy=sqrt(jxh_k(1:Nt,7).^2*16+jxf_k(1:Nt,7).^2*4);
dtau=tau(2)-tau(1);fy(4)=sum(y(1:Nt-1).*exp(2i*m*pi*tau(1:Nt-1)/4))*dtau;
cd ../1dL8beta5U8/
load jxf_k.dat; load jxh_k.dat; Nt=length(jxh_k)/3; tau=jxh_k(1:Nt,1);
y=jxh_k(1:Nt,5)*4+jxf_k(1:Nt,5)*2;
dy=sqrt(jxh_k(1:Nt,7).^2*16+jxf_k(1:Nt,7).^2*4);
dtau=tau(2)-tau(1);fy(5)=sum(y(1:Nt-1).*exp(2i*m*pi*tau(1:Nt-1)/5))*dtau;
cd ../1dL8beta6U8/
load jxf_k.dat; load jxh_k.dat; Nt=length(jxh_k)/3; tau=jxh_k(1:Nt,1);
y=jxh_k(1:Nt,5)*4+jxf_k(1:Nt,5)*2;
dy=sqrt(jxh_k(1:Nt,7).^2*16+jxf_k(1:Nt,7).^2*4);
dtau=tau(2)-tau(1);fy(6)=sum(y(1:Nt-1).*exp(2i*m*pi*tau(1:Nt-1)/6))*dtau;
cd ../1dL8beta7U8/
load jxf_k.dat; load jxh_k.dat; Nt=length(jxh_k)/3; tau=jxh_k(1:Nt,1);
y=jxh_k(1:Nt,5)*4+jxf_k(1:Nt,5)*2;
dy=sqrt(jxh_k(1:Nt,7).^2*16+jxf_k(1:Nt,7).^2*4);
dtau=tau(2)-tau(1);fy(7)=sum(y(1:Nt-1).*exp(2i*m*pi*tau(1:Nt-1)/7))*dtau;
cd ../1dL8beta8U8/
load jxf_k.dat; load jxh_k.dat; Nt=length(jxh_k)/3; tau=jxh_k(1:Nt,1);
y=jxh_k(1:Nt,5)*4+jxf_k(1:Nt,5)*2;
dy=sqrt(jxh_k(1:Nt,7).^2*16+jxf_k(1:Nt,7).^2*4);
dtau=tau(2)-tau(1);fy(8)=sum(y(1:Nt-1).*exp(2i*m*pi*tau(1:Nt-1)/8))*dtau;
cd ../1dL8beta9U8/
load jxf_k.dat; load jxh_k.dat; Nt=length(jxh_k)/3; tau=jxh_k(1:Nt,1);
y=jxh_k(1:Nt,5)*4+jxf_k(1:Nt,5)*2;
dy=sqrt(jxh_k(1:Nt,7).^2*16+jxf_k(1:Nt,7).^2*4);
dtau=tau(2)-tau(1);fy(9)=sum(y(1:Nt-1).*exp(2i*m*pi*tau(1:Nt-1)/9))*dtau;
cd ../1dL8beta10U8/
load jxf_k.dat; load jxh_k.dat; Nt=length(jxh_k)/3; tau=jxh_k(1:Nt,1);
y=jxh_k(1:Nt,5)*4+jxf_k(1:Nt,5)*2;
dy=sqrt(jxh_k(1:Nt,7).^2*16+jxf_k(1:Nt,7).^2*4);
dtau=tau(2)-tau(1);fy(10)=sum(y(1:Nt-1).*exp(2i*m*pi*tau(1:Nt-1)/10))*dtau;
cd ..