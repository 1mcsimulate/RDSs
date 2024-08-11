clear;
%a=[0.2,0.5,0.6,0.3,0.5,0.5,0.5,0.5,0.5,0.5];
%b=[5.6,6.2,3.3,4.1,3.8,5.6,6.2,2.3,3.1,5.9];
%c=[0.93,0.81,0.91,0.78,0.55,0.71,0.82,0.67,0.85,0.85];
d=[1,1,1,1,1,1];
%e=[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5];
D=[4,4,4,4,4,4];
w=0.01;cc=0.1;
A=ones(6);
p=0.6;
q=1.4;
rho_1=0.2;
rho_2=0.6;
rho_3=2.6;
rho_4=1;
rho_5=1;    


B=[0 0.04 0.03 0 0 0
   0.03 0 0.06 0.04 0 0 
   0 0 0 0 0.03 0 
   0 0 0 0 0.04 0.06 
   0 0 0 0 0 0.01 
   0 0 0 0 0.03 0 ];



tau=10;tao=1;xi=1;
%N是tau/h(h为时间上划分的步长)
h=0.01; 
N=tau/h;
N1=tao/h;
%P=10;%每个周期步数
%K=500;%时间方向计算步数
%T=h*K;%时间最值为T
%p=K+N+1;%时间点总数
M=25;%Delta{x}为空间上的划分,M=pi/Delta{x}
dx=10/M;
%初始化，求u，v在每个网格的节点上的值

%%%%%%%%%初始化   先 驱动  在响应  在误差
 y1=zeros(M+1,N);
 y2=zeros(M+1,N);
 y3=zeros(M+1,N);
 y4=zeros(M+1,N);
 y5=zeros(M+1,N);
 y6=zeros(M+1,N);
%y7=zeros(M+1,N);
%y8=zeros(M+1,N);
%y9=zeros(M+1,N);
%y10=zeros(M+1,N);

s1=zeros(M+1,N);
s2=zeros(M+1,N);
s3=zeros(M+1,N);
s4=zeros(M+1,N);
s5=zeros(M+1,N);
s6=zeros(M+1,N);
%s7=zeros(M+1,N);
%s8=zeros(M+1,N);
%s9=zeros(M+1,N);
%s10=zeros(M+1,N);




%%%%%%%%%%%%%%%赋边值
 for j=1:N
     y1(1,j)=0;     y1(M+1,j)=0;
     y2(1,j)=0;     y2(M+1,j)=0;
     y3(1,j)=0;     y3(M+1,j)=0;
     y4(1,j)=0;     y4(M+1,j)=0;
     y5(1,j)=0;     y5(M+1,j)=0;
     y6(1,j)=0;     y6(M+1,j)=0;
     %y7(1,j)=0;     y7(M+1,j)=0;
     %y8(1,j)=0;     y8(M+1,j)=0;
     %y9(1,j)=0;     y9(M+1,j)=0;
     %y10(1,j)=0;    y10(M+1,j)=0;
 
     s1(1,j)=0;     s1(M+1,j)=0;
     s2(1,j)=0;     s2(M+1,j)=0;
     s3(1,j)=0;     s3(M+1,j)=0;
     s4(1,j)=0;     s4(M+1,j)=0;
     s5(1,j)=0;     s5(M+1,j)=0;
     s6(1,j)=0;     s6(M+1,j)=0;
     %s7(1,j)=0;     s7(M+1,j)=0;
     %s8(1,j)=0;     s8(M+1,j)=0;
     %s9(1,j)=0;     s9(M+1,j)=0;
     %s10(1,j)=0;    s10(M+1,j)=0;
     
 end
%%%%%%%%%%%%赋初值
for i=2:M
   for j=1:N1
    y1(i,j)=11.2*sin((i*h*pi)^6);
     y2(i,j)=7*sin((i*h*pi)^4);
     y3(i,j)=1*sin((i*h*pi)^3);
      y4(i,j)=12*sin((i*h*pi)^6);
     y5(i,j)=6*sin((i*h*pi)^4);
     y6(i,j)=1*sin((i*h*pi)^3);
     %y7(i,1)=11.8*sin((i*h*pi)^6);
     %y8(i,1)=7.2*sin((i*h*pi)^4);
     %y9(i,1)=1*sin((i*h*pi)^3);
     %y10(i,1)=13.2*sin((i*h*pi)^6);
      
     s1(i,j)=4*sin((i*h*pi)^4);
     s2(i,j)=14*sin((i*h*pi)^5);
     s3(i,j)=3.8*sin((i*h*pi)^3);
     s4(i,j)=4.2*sin((i*h*pi)^5);
     s5(i,j)=4.8*sin((i*h*pi)^3);
     s6(i,j)=3.4*sin((i*h*pi)^4);
     %s7(i,1)=4.5*sin((i*h*pi)^5);
     %s8(i,1)=3.8*sin((i*h*pi)^3);
     %s9(i,1)=4*sin((i*h*pi)^5);
     %s10(i,1)=3*sin((i*h*pi)^3);
     
                            
end
   
end
%%%%%%%%%%%%积分初始化
yita=[2.8,3.1,4.5,2.4,3.5,3.2];

JD=1;%节点数
JXKZ=zeros(N,JD);%放间歇控制，如果=1就是有控制，=0没有控制
time=0:h:tau;
%接下来是第一个节点的间歇控制
Kk=0;
for t11=1:N 
    if time(t11)>=Kk && time(t11)<Kk+0.4 %事先定好的
        JXKZ(t11,1)=1; %1表示第1、2个节点；=1表示这个时间段内有控制
    end
    if time(t11)>=Kk+0.5 && time(t11)<Kk+0.95 %事先定好的
        JXKZ(t11,1)=1;%1表示第1、2个节点；=1表示这个时间段内有控制
    end
    if time(t11)>Kk+0.95 %事先定好的，看上数三行最后
        Kk=Kk+1;
    end
    if Kk>=tau-h
        break;
    end
end




%原来的f 算控制项
% x1=zeros(M+1,N);
% x2=zeros(M+1,N);
% x3=zeros(M+1,N);
% x4=zeros(M+1,N);
% x5=zeros(M+1,N);
% x6=zeros(M+1,N);
%x7=zeros(M+1,N);
%x8=zeros(M+1,N);
%x9=zeros(M+1,N);
%x10=zeros(M+1,N);

%原来的ff  算耦合项和hanshu f项
% xx1=zeros(M+1,N);
% xx2=zeros(M+1,N);
% xx3=zeros(M+1,N);
% xx4=zeros(M+1,N);
% xx5=zeros(M+1,N);
% xx6=zeros(M+1,N);
%xx7=zeros(M+1,N);
%xx8=zeros(M+1,N);
%xx9=zeros(M+1,N);
%xx10=zeros(M+1,N);
u1=zeros(M+1,N);
u2=zeros(M+1,N);
u3=zeros(M+1,N);
u4=zeros(M+1,N);
u5=zeros(M+1,N);
u6=zeros(M+1,N);
%u7=zeros(M+1,N);
%u8=zeros(M+1,N);
%u9=zeros(M+1,N);
%u10=zeros(M+1,N);


tao()




for j=N1:N-1  
    for i=2:M
       %%%%%%%%%%%%%%%%%%%%%%%%drive system驱动系统  
        xx1(i,j)=atan(y1(i,j ))...
            +(2*A(1,1)*B(1,1)+A(1,2)*B(1,2)+A(1,3)*B(1,3)+A(1,4)*B(1,4)+A(1,5)*B(1,5)+A(1,6)*B(1,6))*(y1(i,j))...
            +A(1,2)*B(1,2)*(y2(i,j))+A(1,3)*B(1,3)*(y3(i,j))+A(1,4)*B(1,4)*(y4(i,j))+A(1,5)*B(1,5)*(y5(i,j))...
            +A(1,6)*B(1,6)*(y6(i,j));
        
        xx2(i,j)=atan(y2(i,j ))...
            +(2*A(2,2)*B(2,2)+A(2,1)*B(2,1)+A(2,3)*B(2,3)+A(2,4)*B(2,4)+A(2,5)*B(2,5)+A(2,6)*B(2,6))*(y2(i,j))...
            +A(2,1)*B(2,1)*(y1(i,j))+A(2,3)*B(2,3)*(y3(i,j))+A(2,4)*B(2,4)*(y4(i,j))+A(2,5)*B(2,5)*(y5(i,j))...
            +A(2,6)*B(2,6)*(y6(i,j));
        
        xx3(i,j)=atan(y3(i,j-tao))...
            +(2*A(3,3)*B(3,3)+A(3,1)*B(3,1)+A(3,2)*B(3,2)+A(3,4)*B(3,4)+A(3,5)*B(3,5)+A(3,6)*B(3,6))*(y3(i,j))...
            +A(3,1)*B(3,1)*(y1(i,j))+A(3,2)*B(3,2)*(y2(i,j))+A(3,4)*B(3,4)*(y4(i,j))+A(3,5)*B(3,5)*(y5(i,j))...
            +A(3,6)*B(3,6)*(y6(i,j));
    
        xx4(i,j)=atan(y4(i,j-tao))...
            +(2*A(4,4)*B(4,4)+A(4,1)*B(4,1)+A(4,2)*B(4,2)+A(4,3)*B(4,3)+A(4,5)*B(4,5)+A(4,6)*B(4,6))*(y4(i,j))...
            +A(4,1)*B(4,1)*(y1(i,j))+A(4,2)*B(4,2)*(y2(i,j))+A(4,3)*B(4,3)*(y3(i,j))+A(4,5)*B(4,5)*(y5(i,j))...
            +A(4,6)*B(4,6)*(y6(i,j));
        
        xx5(i,j)=atan(y5(i,j-tao))...
            +(2*A(5,5)*B(5,5)+A(5,1)*B(5,1)+A(5,2)*B(5,2)+A(5,3)*B(5,3)+A(5,4)*B(5,4)+A(5,6)*B(5,6))*(y5(i,j))...
            +A(5,1)*B(5,1)*(y1(i,j))+A(5,2)*B(5,2)*(y2(i,j))+A(5,3)*B(5,3)*(y3(i,j))+A(5,4)*B(5,4)*(y4(i,j))...
            +A(5,6)*B(5,6)*(y6(i,j));
        
        xx6(i,j)=atan(y6(i,j-tao))...
            +(2*A(6,6)*B(6,6)+A(6,1)*B(6,1)+A(6,2)*B(6,2)+A(6,3)*B(6,3)+A(6,4)*B(6,4)+A(6,5)*B(6,5)  )*(y6(i,j))...
            +A(6,1)*B(6,1)*(y1(i,j))+A(6,2)*B(6,2)*(y2(i,j))+A(6,3)*B(6,3)*(y3(i,j))+A(6,4)*B(6,4)*(y4(i,j))...
            +A(6,5)*B(6,5)*(y5(i,j));
        
        
     %wei gai quan
        y1(i,j+1)=y1(i,j)+h*D(1)*(y1(i+1,j)-2*y1(i,j)+y1(i-1,j))/dx/dx+h*xx1(i,j);
        y2(i,j+1)=y2(i,j)+h*D(2)*(y2(i+1,j)-2*y2(i,j)+y2(i-1,j))/dx/dx+h*xx2(i,j);
        y3(i,j+1)=y3(i,j)+h*D(3)*(y3(i+1,j)-2*y3(i,j)+y3(i-1,j))/dx/dx+h*xx3(i,j);
        y4(i,j+1)=y4(i,j)+h*D(4)*(y4(i+1,j)-2*y4(i,j)+y4(i-1,j))/dx/dx+h*xx4(i,j);
        y5(i,j+1)=y5(i,j)+h*D(5)*(y5(i+1,j)-2*y5(i,j)+y5(i-1,j))/dx/dx+h*xx5(i,j);
        y6(i,j+1)=y6(i,j)+h*D(6)*(y6(i+1,j)-2*y6(i,j)+y6(i-1,j))/dx/dx+h*xx6(i,j);
       
 
        
        
        %%%%%%%%%%%response system 控制器出现
       
          x1(i,j)=+atan(s1(i,j-tao))...
            +(2*A(1,1)*B(1,1)+A(1,2)*B(1,2)+A(1,3)*B(1,3)+A(1,4)*B(1,4)+A(1,5)*B(1,5)+A(1,6)*B(1,6))*(s1(i,j))...
            +A(1,2)*B(1,2)*(s2(i,j))+A(1,3)*B(1,3)*(s3(i,j))+A(1,4)*B(1,4)*(s4(i,j))+A(1,5)*B(1,5)*(s5(i,j))...
            +A(1,6)*B(1,6)*(s6(i,j));
      
      
      
          x2(i,j)=+atan(s2(i,j-tao))...
            +(2*A(2,2)*B(2,2)+A(2,1)*B(2,1)+A(2,3)*B(2,3)+A(2,4)*B(2,4)+A(2,5)*B(2,5)+A(2,6)*B(2,6))*(s2(i,j))...
            +A(2,1)*B(2,1)*(s1(i,j))+A(2,3)*B(2,3)*(s3(i,j))+A(2,4)*B(2,4)*(s4(i,j))+A(2,5)*B(2,5)*(s5(i,j))...
            +A(2,6)*B(2,6)*(s6(i,j));
      
      
          x3(i,j)=+atan(s3(i,j-tao))...
            +(2*A(3,3)*B(3,3)+A(3,1)*B(3,1)+A(3,2)*B(3,2)+A(3,4)*B(3,4)+A(3,5)*B(3,5)+A(3,6)*B(3,6))*(s3(i,j))...
            +A(3,1)*B(3,1)*(s1(i,j))+A(3,2)*B(3,2)*(s2(i,j))+A(3,4)*B(3,4)*(s4(i,j))+A(3,5)*B(3,5)*(s5(i,j))...
            +A(3,6)*B(3,6)*(s6(i,j));
      
      
          x4(i,j)=+atan(s4(i,j-tao))...
            +(2*A(4,4)*B(4,4)+A(4,1)*B(4,1)+A(4,2)*B(4,2)+A(4,3)*B(4,3)+A(4,5)*B(4,5)+A(4,6)*B(4,6))*(s4(i,j))...
            +A(4,1)*B(4,1)*(s1(i,j))+A(4,2)*B(4,2)*(s2(i,j))+A(4,3)*B(4,3)*(s3(i,j))+A(4,5)*B(4,5)*(s5(i,j))...
            +A(4,6)*B(4,6)*(s6(i,j));
          
          
          x5(i,j)=+atan(s5(i,j-tao))...
             +(2*A(5,5)*B(5,5)+A(5,1)*B(5,1)+A(5,2)*B(5,2)+A(5,3)*B(5,3)+A(5,4)*B(5,4)+A(5,6)*B(5,6))*(s5(i,j))...
            +A(5,1)*B(5,1)*(s1(i,j))+A(5,2)*B(5,2)*(s2(i,j))+A(5,3)*B(5,3)*(s3(i,j))+A(5,4)*B(5,4)*(s4(i,j))...
            +A(5,6)*B(5,6)*(s6(i,j));
          
          
         
          x6(i,j)=+atan(s6(i,j-tao))...
            +(2*A(6,6)*B(6,6)+A(6,1)*B(6,1)+A(6,2)*B(6,2)+A(6,3)*B(6,3)+A(6,4)*B(6,4)+A(6,5)*B(6,5))*(s6(i,j))...
            +A(6,1)*B(6,1)*(s1(i,j))+A(6,2)*B(6,2)*(s2(i,j))+A(6,3)*B(6,3)*(s3(i,j))+A(6,4)*B(6,4)*(s4(i,j))...
            +A(6,5)*B(6,5)*(s5(i,j));
        
       
   
          
          s1(i,j+1)=s1(i,j)+h*D(1)*(s1(i+1,j)-2*s1(i,j)+s1(i-1,j))/dx/dx+h*x1(i,j)-yita(1)*JXKZ(j,1)*(s1(i-1,j)-y1(i-1,j))*h-rho_4*JXKZ(j,1)*sign(s1(i-1,j)-y1(i-1,j))*(abs(s1(i-1,j)-y1(i-1,j)))^p*h-rho_5*JXKZ(j,1)*sign(s1(i-1,j)-y1(i-1,j))*abs((s1(i-1,j)-y1(i-1,j)))^q*h-JXKZ(j,1)*(rho_1*((s1(i-1,j)-y1(i-1,j))^2+(s1(i-1,j-tao)-y1(i-1,j-tao))^2)*h/2)^((1+p)/2)*sign(s1(i-1,j)-y1(i-1,j))*h-JXKZ(j,1)*(rho_2*((s1(i-1,j)-y1(i-1,j))^2+(s1(i-1,j-tao)-y1(i-1,j-tao))^2)*h/2)^((1+q)/2)*sign(s1(i-1,j)-y1(i-1,j))*h-JXKZ(j,1)*(rho_3*((s1(i-1,j)-y1(i-1,j))^2+(s1(i-1,j-tao)-y1(i-1,j-tao))^2)*h/2)^1*sign(s1(i-1,j)-y1(i-1,j))*h;
          s2(i,j+1)=s2(i,j)+h*D(2)*(s2(i+1,j)-2*s2(i,j)+s2(i-1,j))/dx/dx+h*x2(i,j)-yita(2)*JXKZ(j,1)*(s2(i-1,j)-y2(i-1,j))*h-rho_4*JXKZ(j,1)*sign(s2(i-1,j)-y2(i-1,j))*(abs(s2(i-1,j)-y2(i-1,j)))^p*h-rho_5*JXKZ(j,1)*sign(s2(i-1,j)-y2(i-1,j))*abs((s2(i-1,j)-y2(i-1,j)))^q*h-JXKZ(j,1)*(rho_1*((s2(i-1,j)-y2(i-1,j))^2+(s2(i-1,j-tao)-y2(i-1,j-tao))^2)*h/2)^((1+p)/2)*sign(s2(i-1,j)-y2(i-1,j))*h-JXKZ(j,1)*(rho_2*((s2(i-1,j)-y2(i-1,j))^2+(s2(i-1,j-tao)-y2(i-1,j-tao))^2)*h/2)^((1+q)/2)*sign(s2(i-1,j)-y2(i-1,j))*h-JXKZ(j,1)*(rho_3*((s2(i-1,j)-y2(i-1,j))^2+(s2(i-1,j-tao)-y2(i-1,j-tao))^2)*h/2)^1*sign(s2(i-1,j)-y2(i-1,j))*h;
          s3(i,j+1)=s3(i,j)+h*D(3)*(s3(i+1,j)-2*s3(i,j)+s3(i-1,j))/dx/dx+h*x3(i,j)-yita(3)*JXKZ(j,1)*(s3(i-1,j)-y3(i-1,j))*h-rho_4*JXKZ(j,1)*sign(s3(i-1,j)-y3(i-1,j))*(abs(s3(i-1,j)-y3(i-1,j)))^p*h-rho_5*JXKZ(j,1)*sign(s3(i-1,j)-y3(i-1,j))*abs((s3(i-1,j)-y3(i-1,j)))^q*h-JXKZ(j,1)*(rho_1*((s3(i-1,j)-y3(i-1,j))^2+(s3(i-1,j-tao)-y3(i-1,j-tao))^2)*h/2)^((1+p)/2)*sign(s3(i-1,j)-y3(i-1,j))*h-JXKZ(j,1)*(rho_2*((s3(i-1,j)-y3(i-1,j))^2+(s3(i-1,j-tao)-y3(i-1,j-tao))^2)*h/2)^((1+q)/2)*sign(s3(i-1,j)-y3(i-1,j))*h-JXKZ(j,1)*(rho_3*((s3(i-1,j)-y3(i-1,j))^2+(s3(i-1,j-tao)-y3(i-1,j-tao))^2)*h/2)^1*sign(s3(i-1,j)-y3(i-1,j))*h;
          s4(i,j+1)=s4(i,j)+h*D(4)*(s4(i+1,j)-2*s4(i,j)+s4(i-1,j))/dx/dx+h*x4(i,j)-yita(4)*JXKZ(j,1)*(s4(i-1,j)-y4(i-1,j))*h-rho_4*JXKZ(j,1)*sign(s4(i-1,j)-y4(i-1,j))*(abs(s4(i-1,j)-y4(i-1,j)))^p*h-rho_5*JXKZ(j,1)*sign(s4(i-1,j)-y4(i-1,j))*abs((s4(i-1,j)-y4(i-1,j)))^q*h-JXKZ(j,1)*(rho_1*((s4(i-1,j)-y4(i-1,j))^2+(s4(i-1,j-tao)-y4(i-1,j-tao))^2)*h/2)^((1+p)/2)*sign(s4(i-1,j)-y4(i-1,j))*h-JXKZ(j,1)*(rho_2*((s4(i-1,j)-y4(i-1,j))^2+(s4(i-1,j-tao)-y4(i-1,j-tao))^2)*h/2)^((1+q)/2)*sign(s4(i-1,j)-y4(i-1,j))*h-JXKZ(j,1)*(rho_3*((s4(i-1,j)-y4(i-1,j))^2+(s4(i-1,j-tao)-y4(i-1,j-tao))^2)*h/2)^1*sign(s4(i-1,j)-y4(i-1,j))*h;
          s5(i,j+1)=s5(i,j)+h*D(5)*(s5(i+1,j)-2*s5(i,j)+s5(i-1,j))/dx/dx+h*x5(i,j)-yita(5)*JXKZ(j,1)*(s5(i-1,j)-y5(i-1,j))*h-rho_4*JXKZ(j,1)*sign(s5(i-1,j)-y5(i-1,j))*(abs(s5(i-1,j)-y5(i-1,j)))^p*h-rho_5*JXKZ(j,1)*sign(s5(i-1,j)-y5(i-1,j))*abs((s5(i-1,j)-y5(i-1,j)))^q*h-JXKZ(j,1)*(rho_1*((s5(i-1,j)-y5(i-1,j))^2+(s5(i-1,j-tao)-y5(i-1,j-tao))^2)*h/2)^((1+p)/2)*sign(s5(i-1,j)-y5(i-1,j))*h-JXKZ(j,1)*(rho_2*((s5(i-1,j)-y5(i-1,j))^2+(s5(i-1,j-tao)-y5(i-1,j-tao))^2)*h/2)^((1+q)/2)*sign(s5(i-1,j)-y5(i-1,j))*h-JXKZ(j,1)*(rho_3*((s5(i-1,j)-y5(i-1,j))^2+(s5(i-1,j-tao)-y5(i-1,j-tao))^2)*h/2)^1*sign(s5(i-1,j)-y5(i-1,j))*h;
          s6(i,j+1)=s6(i,j)+h*D(6)*(s6(i+1,j)-2*s6(i,j)+s6(i-1,j))/dx/dx+h*x6(i,j)-yita(6)*JXKZ(j,1)*(s6(i-1,j)-y6(i-1,j))*h-rho_4*JXKZ(j,1)*sign(s6(i-1,j)-y6(i-1,j))*(abs(s6(i-1,j)-y6(i-1,j)))^p*h-rho_5*JXKZ(j,1)*sign(s6(i-1,j)-y6(i-1,j))*abs((s6(i-1,j)-y6(i-1,j)))^q*h-JXKZ(j,1)*(rho_1*((s6(i-1,j)-y6(i-1,j))^2+(s6(i-1,j-tao)-y6(i-1,j-tao))^2)*h/2)^((1+p)/2)*sign(s6(i-1,j)-y6(i-1,j))*h-JXKZ(j,1)*(rho_2*((s6(i-1,j)-y6(i-1,j))^2+(s6(i-1,j-tao)-y6(i-1,j-tao))^2)*h/2)^((1+q)/2)*sign(s6(i-1,j)-y6(i-1,j))*h-JXKZ(j,1)*(rho_3*((s6(i-1,j)-y6(i-1,j))^2+(s6(i-1,j-tao)-y6(i-1,j-tao))^2)*h/2)^1*sign(s6(i-1,j)-y6(i-1,j))*h;
          %s7(i,j+1)=s7(i,j)+h*D(7)*(s7(i+1,j)-2*s7(i,j)+s7(i-1,j))/dx/dx+h*x7(i,j)-yita(7)*JXKZ(j,1)*(s7(i-1,j)-y7(i-1,j))*h;
          %s8(i,j+1)=s8(i,j)+h*D(8)*(s8(i+1,j)-2*s8(i,j)+s8(i-1,j))/dx/dx+h*x8(i,j)-yita(8)*JXKZ(j,1)*(s8(i-1,j)-y8(i-1,j))*h;
          %s9(i,j+1)=s9(i,j)+h*D(9)*(s9(i+1,j)-2*s9(i,j)+s9(i-1,j))/dx/dx+h*x9(i,j)-yita(9)*JXKZ(j,1)*(s9(i-1,j)-y9(i-1,j))*h;
          %s10(i,j+1)=s10(i,j)+h*D(10)*(s10(i+1,j)-2*s10(i,j)+s10(i-1,j))/dx/dx+h*x10(i,j)-yita(10)*JXKZ(j,1)*(s10(i-1,j)-y10(i-1,j))*h;
          
          u1(i,j+1)=-yita(1)*JXKZ(j,1)*(s1(i-1,j)-y1(i-1,j))*h-rho_4*JXKZ(j,1)*sign(s1(i-1,j)-y1(i-1,j))*(abs(s1(i-1,j)-y1(i-1,j)))^p*h-rho_5*JXKZ(j,1)*sign(s1(i-1,j)-y1(i-1,j))*abs((s1(i-1,j)-y1(i-1,j)))^q*h-JXKZ(j,1)*(rho_1*((s1(i-1,j)-y1(i-1,j))^2+(s1(i-1,j-tao)-y1(i-1,j-tao))^2)*h/2)^((1+p)/2)*sign(s1(i-1,j)-y1(i-1,j))*h-JXKZ(j,1)*(rho_2*((s1(i-1,j)-y1(i-1,j))^2+(s1(i-1,j-tao)-y1(i-1,j-tao))^2)*h/2)^((1+q)/2)*sign(s1(i-1,j)-y1(i-1,j))*h-JXKZ(j,1)*(rho_3*((s1(i-1,j)-y1(i-1,j))^2+(s1(i-1,j-tao)-y1(i-1,j-tao))^2)*h/2)^1*sign(s1(i-1,j)-y1(i-1,j))*h;
          u2(i,j+1)=-yita(2)*JXKZ(j,1)*(s2(i-1,j)-y2(i-1,j))*h-rho_4*JXKZ(j,1)*sign(s2(i-1,j)-y2(i-1,j))*(abs(s2(i-1,j)-y2(i-1,j)))^p*h-rho_5*JXKZ(j,1)*sign(s2(i-1,j)-y2(i-1,j))*abs((s2(i-1,j)-y2(i-1,j)))^q*h-JXKZ(j,1)*(rho_1*((s2(i-1,j)-y2(i-1,j))^2+(s2(i-1,j-tao)-y2(i-1,j-tao))^2)*h/2)^((1+p)/2)*sign(s2(i-1,j)-y2(i-1,j))*h-JXKZ(j,1)*(rho_2*((s2(i-1,j)-y2(i-1,j))^2+(s2(i-1,j-tao)-y2(i-1,j-tao))^2)*h/2)^((1+q)/2)*sign(s2(i-1,j)-y2(i-1,j))*h-JXKZ(j,1)*(rho_3*((s2(i-1,j)-y2(i-1,j))^2+(s2(i-1,j-tao)-y2(i-1,j-tao))^2)*h/2)^1*sign(s2(i-1,j)-y2(i-1,j))*h;
          u3(i,j+1)=-yita(3)*JXKZ(j,1)*(s3(i-1,j)-y3(i-1,j))*h-rho_4*JXKZ(j,1)*sign(s3(i-1,j)-y3(i-1,j))*(abs(s3(i-1,j)-y3(i-1,j)))^p*h-rho_5*JXKZ(j,1)*sign(s3(i-1,j)-y3(i-1,j))*abs((s3(i-1,j)-y3(i-1,j)))^q*h-JXKZ(j,1)*(rho_1*((s3(i-1,j)-y3(i-1,j))^2+(s3(i-1,j-tao)-y3(i-1,j-tao))^2)*h/2)^((1+p)/2)*sign(s3(i-1,j)-y3(i-1,j))*h-JXKZ(j,1)*(rho_2*((s3(i-1,j)-y3(i-1,j))^2+(s3(i-1,j-tao)-y3(i-1,j-tao))^2)*h/2)^((1+q)/2)*sign(s3(i-1,j)-y3(i-1,j))*h-JXKZ(j,1)*(rho_3*((s3(i-1,j)-y3(i-1,j))^2+(s3(i-1,j-tao)-y3(i-1,j-tao))^2)*h/2)^1*sign(s3(i-1,j)-y3(i-1,j))*h;
          u4(i,j+1)=-yita(4)*JXKZ(j,1)*(s4(i-1,j)-y4(i-1,j))*h-rho_4*JXKZ(j,1)*sign(s4(i-1,j)-y4(i-1,j))*(abs(s4(i-1,j)-y4(i-1,j)))^p*h-rho_5*JXKZ(j,1)*sign(s4(i-1,j)-y4(i-1,j))*abs((s4(i-1,j)-y4(i-1,j)))^q*h-JXKZ(j,1)*(rho_1*((s4(i-1,j)-y4(i-1,j))^2+(s4(i-1,j-tao)-y4(i-1,j-tao))^2)*h/2)^((1+p)/2)*sign(s4(i-1,j)-y4(i-1,j))*h-JXKZ(j,1)*(rho_2*((s4(i-1,j)-y4(i-1,j))^2+(s4(i-1,j-tao)-y4(i-1,j-tao))^2)*h/2)^((1+q)/2)*sign(s4(i-1,j)-y4(i-1,j))*h-JXKZ(j,1)*(rho_3*((s4(i-1,j)-y4(i-1,j))^2+(s4(i-1,j-tao)-y4(i-1,j-tao))^2)*h/2)^1*sign(s4(i-1,j)-y4(i-1,j))*h;
          u5(i,j+1)=-yita(5)*JXKZ(j,1)*(s5(i-1,j)-y5(i-1,j))*h-rho_4*JXKZ(j,1)*sign(s5(i-1,j)-y5(i-1,j))*(abs(s5(i-1,j)-y5(i-1,j)))^p*h-rho_5*JXKZ(j,1)*sign(s5(i-1,j)-y5(i-1,j))*abs((s5(i-1,j)-y5(i-1,j)))^q*h-JXKZ(j,1)*(rho_1*((s5(i-1,j)-y5(i-1,j))^2+(s5(i-1,j-tao)-y5(i-1,j-tao))^2)*h/2)^((1+p)/2)*sign(s5(i-1,j)-y5(i-1,j))*h-JXKZ(j,1)*(rho_2*((s5(i-1,j)-y5(i-1,j))^2+(s5(i-1,j-tao)-y5(i-1,j-tao))^2)*h/2)^((1+q)/2)*sign(s5(i-1,j)-y5(i-1,j))*h-JXKZ(j,1)*(rho_3*((s5(i-1,j)-y5(i-1,j))^2+(s5(i-1,j-tao)-y5(i-1,j-tao))^2)*h/2)^1*sign(s5(i-1,j)-y5(i-1,j))*h;
          u6(i,j+1)=-yita(6)*JXKZ(j,1)*(s6(i-1,j)-y6(i-1,j))*h-rho_4*JXKZ(j,1)*sign(s6(i-1,j)-y6(i-1,j))*(abs(s6(i-1,j)-y6(i-1,j)))^p*h-rho_5*JXKZ(j,1)*sign(s6(i-1,j)-y6(i-1,j))*abs((s6(i-1,j)-y6(i-1,j)))^q*h-JXKZ(j,1)*(rho_1*((s6(i-1,j)-y6(i-1,j))^2+(s6(i-1,j-tao)-y6(i-1,j-tao))^2)*h/2)^((1+p)/2)*sign(s6(i-1,j)-y6(i-1,j))*h-JXKZ(j,1)*(rho_2*((s6(i-1,j)-y6(i-1,j))^2+(s6(i-1,j-tao)-y6(i-1,j-tao))^2)*h/2)^((1+q)/2)*sign(s6(i-1,j)-y6(i-1,j))*h-JXKZ(j,1)*(rho_3*((s6(i-1,j)-y6(i-1,j))^2+(s6(i-1,j-tao)-y6(i-1,j-tao))^2)*h/2)^1*sign(s6(i-1,j)-y6(i-1,j))*h;
          %u7(i,j+1)=-yita(7)*JXKZ(j,1)*(s7(i-1,j)-y7(i-1,j))-rho_4*sign(JXKZ(j,1)*(s7(i-1,j)-y7(i-1,j)))*abs(JXKZ(j,1)*(s7(i-1,j)-y7(i-1,j)))^p-rho_5*sign(JXKZ(j,1)*(s7(i-1,j)-y7(i-1,j)))*abs(JXKZ(j,1)*(s7(i-1,j)-y7(i-1,j)))^q;
          %u8(i,j+1)=-yita(8)*JXKZ(j,1)*(s8(i-1,j)-y8(i-1,j))-rho_4*sign(JXKZ(j,1)*(s8(i-1,j)-y8(i-1,j)))*abs(JXKZ(j,1)*(s8(i-1,j)-y8(i-1,j)))^p-rho_5*sign(JXKZ(j,1)*(s8(i-1,j)-y8(i-1,j)))*abs(JXKZ(j,1)*(s8(i-1,j)-y8(i-1,j)))^q;
          %u9(i,j+1)=-yita(9)*JXKZ(j,1)*(s9(i-1,j)-y9(i-1,j))-rho_4*sign(JXKZ(j,1)*(s9(i-1,j)-y9(i-1,j)))*abs(JXKZ(j,1)*(s9(i-1,j)-y9(i-1,j)))^p-rho_5*sign(JXKZ(j,1)*(s9(i-1,j)-y9(i-1,j)))*abs(JXKZ(j,1)*(s9(i-1,j)-y9(i-1,j)))^q;
          %u10(i,j+1)=-yita(10)*JXKZ(j,1)*(s10(i-1,j)-y10(i-1,j))-rho_4*sign(JXKZ(j,1)*(s10(i-1,j)-y10(i-1,j)))*abs(JXKZ(j,1)*(s10(i-1,j)-y10(i-1,j)))^p-rho_5*sign(JXKZ(j,1)*(s10(i-1,j)-y10(i-1,j)))*abs(JXKZ(j,1)*(s10(i-1,j)-y10(i-1,j)))^q;
      
        end    
      end
       
    
%end





x=-5:dx:5;
t=-1:h:tau-1-h;





%figure
%mesh(t,x,err3),xlabel('t'),ylabel('x'),zlabel('e_3');
%axis([-3.5 5 -5 5 -0.8 0.8])
%subplot(2,3,1),mesh(t,x,err1)
%subplot(2,3,3),mesh(t,x,err2)
%subplot(2,3,5),mesh(t,x,err3)




figure
subplot(2,2,1),mesh(t,x,y1),xlabel('t'),ylabel('x'),zlabel('u_1')
subplot(2,2,2),mesh(t,x,s1),xlabel('t'),ylabel('x'),zlabel('v_1')
subplot(2,2,3),mesh(t,x,y2),xlabel('t'),ylabel('x'),zlabel('u_2')
subplot(2,2,4),mesh(t,x,s2),xlabel('t'),ylabel('x'),zlabel('v_2')
%xlim([-1 20])

figure
subplot(2,2,1),mesh(t,x,y3),xlabel('t'),ylabel('x'),zlabel('u_3')
subplot(2,2,2),mesh(t,x,s3),xlabel('t'),ylabel('x'),zlabel('v_3')
subplot(2,2,3),mesh(t,x,y4),xlabel('t'),ylabel('x'),zlabel('u_4')
subplot(2,2,4),mesh(t,x,s4),xlabel('t'),ylabel('x'),zlabel('v_4')

figure
subplot(2,2,1),mesh(t,x,y5),xlabel('t'),ylabel('x'),zlabel('u_5')
subplot(2,2,2),mesh(t,x,s5),xlabel('t'),ylabel('x'),zlabel('v_5')
subplot(2,2,3),mesh(t,x,y6),xlabel('t'),ylabel('x'),zlabel('u_6')
subplot(2,2,4),mesh(t,x,s6),xlabel('t'),ylabel('x'),zlabel('v_6')


%figure
%subplot(2,2,1),mesh(t,x,y7),xlabel('t'),ylabel('s'),zlabel('z_7')
%subplot(2,2,2),mesh(t,x,s7),xlabel('t'),ylabel('s'),zlabel('w_7')
%subplot(2,2,3),mesh(t,x,y8),xlabel('t'),ylabel('s'),zlabel('z_8')
%subplot(2,2,4),mesh(t,x,s8),xlabel('t'),ylabel('s'),zlabel('w_8')

%figure
%subplot(2,2,1),mesh(t,x,y9),xlabel('t'),ylabel('s'),zlabel('z_9')
%subplot(2,2,2),mesh(t,x,s9),xlabel('t'),ylabel('s'),zlabel('w_9')
%subplot(2,2,3),mesh(t,x,y10),xlabel('t'),ylabel('s'),zlabel('z_{10}')
%subplot(2,2,4),mesh(t,x,s10),xlabel('t'),ylabel('s'),zlabel('w_{10}')




figure
subplot(2,2,1),mesh(t,x,s1-y1),xlabel('t'),ylabel('x'),zlabel('z_1')
subplot(2,2,2),mesh(t,x,s2-y2),xlabel('t'),ylabel('x'),zlabel('z_2')
subplot(2,2,3),mesh(t,x,s3-y3),xlabel('t'),ylabel('x'),zlabel('z_3')
subplot(2,2,4),mesh(t,x,s4-y4),xlabel('t'),ylabel('x'),zlabel('z_4')
figure
subplot(2,2,1),mesh(t,x,s5-y5),xlabel('t'),ylabel('x'),zlabel('z_5')
subplot(2,2,2),mesh(t,x,s6-y6),xlabel('t'),ylabel('x'),zlabel('z_6')
%subplot(2,2,3),mesh(t,x,s7-y7),xlabel('t'),ylabel('s'),zlabel('v_7')
%subplot(2,2,4),mesh(t,x,s8-y8),xlabel('t'),ylabel('s'),zlabel('v_8')
%figure
%subplot(2,2,1),mesh(t,x,s9-y9),xlabel('t'),ylabel('s'),zlabel('v_9')
%subplot(2,2,2),mesh(t,x,s10-y10),xlabel('t'),ylabel('s'),zlabel('v_{10}')
subplot(2,2,3),mesh(t,x,u1),xlabel('t'),ylabel('x'),zlabel('w_1')
subplot(2,2,4),mesh(t,x,u2),xlabel('t'),ylabel('x'),zlabel('w_2')
figure
subplot(2,2,1),mesh(t,x,u3),xlabel('t'),ylabel('x'),zlabel('w_3')
subplot(2,2,2),mesh(t,x,u4),xlabel('t'),ylabel('x'),zlabel('w_4')
subplot(2,2,3),mesh(t,x,u5),xlabel('t'),ylabel('x'),zlabel('w_5')
subplot(2,2,4),mesh(t,x,u6),xlabel('t'),ylabel('x'),zlabel('w_6')
%figure
%subplot(2,2,1),mesh(t,x,u7),xlabel('t'),ylabel('s'),zlabel('e_7')
%subplot(2,2,2),mesh(t,x,u8),xlabel('t'),ylabel('s'),zlabel('e_8')
%subplot(2,2,3),mesh(t,x,u9),xlabel('t'),ylabel('s'),zlabel('e_9')
%subplot(2,2,4),mesh(t,x,u10),xlabel('t'),ylabel('s'),zlabel('e_{10}')



