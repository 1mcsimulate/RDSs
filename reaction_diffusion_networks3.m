clear;
%a=[0.2,0.5,0.6,0.3,0.5,0.5,0.5,0.5,0.5,0.5];
%b=[5.6,6.2,3.3,4.1,3.8,5.6,6.2,2.3,3.1,5.9];
%c=[0.93,0.81,0.91,0.78,0.55,0.71,0.82,0.67,0.85,0.85];
d=[1,1];
%e=[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5];
D=[0.3,0.3];
w=0.01;cc=0.1;
A=ones(2);
p=0.6;
q=1.4;
rho_1=0.2;
rho_2=0.6;
rho_3=2.6;
rho_4=1;
rho_5=1;    




tau=10;tao=100;xi=1;
%N是tau/h(h为时间上划分的步长)
h=0.01; 
N=tau/h;
N1=tao;
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

s1=zeros(M+1,N);
s2=zeros(M+1,N);




%%%%%%%%%%%%%%%赋边值
 for j=1:N
     y1(1,j)=0;     y1(M+1,j)=0;
     y2(1,j)=0;     y2(M+1,j)=0;
     s1(1,j)=0;     s1(M+1,j)=0;
     s2(1,j)=0;     s2(M+1,j)=0;
     
 end
%%%%%%%%%%%%赋初值
for i=2:M
   for j=1:N1
    y1(i,j)=1.2*sin((i*h*pi)^2);
     y2(i,j)=3.2*sin((i*h*pi)^4);      
     s1(i,j)=4*sin((i*h*pi)^4);
     s2(i,j)=14*sin((i*h*pi)^5);                            
   end   
end
%%%%%%%%%%%%积分初始化
yita=[2.8,3.1];

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
u1=zeros(M+1,N);
u2=zeros(M+1,N);
%u7=zeros(M+1,N);
%u8=zeros(M+1,N);
%u9=zeros(M+1,N);
%u10=zeros(M+1,N);


for j=N1:N-1  
    for i=2:M
        tao=round((exp(j*0.01)/(1+exp(j*0.01)))/0.01);
       %%%%%%%%%%%%%%%%%%%%%%%%drive system驱动系统  
        xx1(i,j)=tanh(y1(i,j-tao ));   
        xx2(i,j)=tanh(y2(i,j-tao ));     
     %wei gai quan
        y1(i,j+1)=y1(i,j)+h*D(1)*(y1(i+1,j)-2*y1(i,j)+y1(i-1,j))/dx/dx+h*xx1(i,j);
        y2(i,j+1)=y2(i,j)+h*D(2)*(y2(i+1,j)-2*y2(i,j)+y2(i-1,j))/dx/dx+h*xx2(i,j);   
        %%%%%%%%%%%response system 控制器出现
          x1(i,j)=+tanh(s1(i,j-tao)); 
          x2(i,j)=+tanh(s2(i,j-tao));  
          s1(i,j+1)=s1(i,j)+h*D(1)*(s1(i+1,j)-2*s1(i,j)+s1(i-1,j))/dx/dx+h*x1(i,j)-yita(1)*JXKZ(j,1)*(s1(i-1,j)-y1(i-1,j))*h-rho_4*JXKZ(j,1)*sign(s1(i-1,j)-y1(i-1,j))*(abs(s1(i-1,j)-y1(i-1,j)))^p*h-rho_5*JXKZ(j,1)*sign(s1(i-1,j)-y1(i-1,j))*abs((s1(i-1,j)-y1(i-1,j)))^q*h-JXKZ(j,1)*(rho_1*((s1(i-1,j)-y1(i-1,j))^2+(s1(i-1,j-tao)-y1(i-1,j-tao))^2)*h/2)^((1+p)/2)*sign(s1(i-1,j)-y1(i-1,j))*h-JXKZ(j,1)*(rho_2*((s1(i-1,j)-y1(i-1,j))^2+(s1(i-1,j-tao)-y1(i-1,j-tao))^2)*h/2)^((1+q)/2)*sign(s1(i-1,j)-y1(i-1,j))*h-JXKZ(j,1)*(rho_3*((s1(i-1,j)-y1(i-1,j))^2+(s1(i-1,j-tao)-y1(i-1,j-tao))^2)*h/2)^1*sign(s1(i-1,j)-y1(i-1,j))*h;
          s2(i,j+1)=s2(i,j)+h*D(2)*(s2(i+1,j)-2*s2(i,j)+s2(i-1,j))/dx/dx+h*x2(i,j)-yita(2)*JXKZ(j,1)*(s2(i-1,j)-y2(i-1,j))*h-rho_4*JXKZ(j,1)*sign(s2(i-1,j)-y2(i-1,j))*(abs(s2(i-1,j)-y2(i-1,j)))^p*h-rho_5*JXKZ(j,1)*sign(s2(i-1,j)-y2(i-1,j))*abs((s2(i-1,j)-y2(i-1,j)))^q*h-JXKZ(j,1)*(rho_1*((s2(i-1,j)-y2(i-1,j))^2+(s2(i-1,j-tao)-y2(i-1,j-tao))^2)*h/2)^((1+p)/2)*sign(s2(i-1,j)-y2(i-1,j))*h-JXKZ(j,1)*(rho_2*((s2(i-1,j)-y2(i-1,j))^2+(s2(i-1,j-tao)-y2(i-1,j-tao))^2)*h/2)^((1+q)/2)*sign(s2(i-1,j)-y2(i-1,j))*h-JXKZ(j,1)*(rho_3*((s2(i-1,j)-y2(i-1,j))^2+(s2(i-1,j-tao)-y2(i-1,j-tao))^2)*h/2)^1*sign(s2(i-1,j)-y2(i-1,j))*h;
          
          u1(i,j+1)=-yita(1)*JXKZ(j,1)*(s1(i-1,j)-y1(i-1,j))*h-rho_4*JXKZ(j,1)*sign(s1(i-1,j)-y1(i-1,j))*(abs(s1(i-1,j)-y1(i-1,j)))^p*h-rho_5*JXKZ(j,1)*sign(s1(i-1,j)-y1(i-1,j))*abs((s1(i-1,j)-y1(i-1,j)))^q*h-JXKZ(j,1)*(rho_1*((s1(i-1,j)-y1(i-1,j))^2+(s1(i-1,j-tao)-y1(i-1,j-tao))^2)*h/2)^((1+p)/2)*sign(s1(i-1,j)-y1(i-1,j))*h-JXKZ(j,1)*(rho_2*((s1(i-1,j)-y1(i-1,j))^2+(s1(i-1,j-tao)-y1(i-1,j-tao))^2)*h/2)^((1+q)/2)*sign(s1(i-1,j)-y1(i-1,j))*h-JXKZ(j,1)*(rho_3*((s1(i-1,j)-y1(i-1,j))^2+(s1(i-1,j-tao)-y1(i-1,j-tao))^2)*h/2)^1*sign(s1(i-1,j)-y1(i-1,j))*h;
          u2(i,j+1)=-yita(2)*JXKZ(j,1)*(s2(i-1,j)-y2(i-1,j))*h-rho_4*JXKZ(j,1)*sign(s2(i-1,j)-y2(i-1,j))*(abs(s2(i-1,j)-y2(i-1,j)))^p*h-rho_5*JXKZ(j,1)*sign(s2(i-1,j)-y2(i-1,j))*abs((s2(i-1,j)-y2(i-1,j)))^q*h-JXKZ(j,1)*(rho_1*((s2(i-1,j)-y2(i-1,j))^2+(s2(i-1,j-tao)-y2(i-1,j-tao))^2)*h/2)^((1+p)/2)*sign(s2(i-1,j)-y2(i-1,j))*h-JXKZ(j,1)*(rho_2*((s2(i-1,j)-y2(i-1,j))^2+(s2(i-1,j-tao)-y2(i-1,j-tao))^2)*h/2)^((1+q)/2)*sign(s2(i-1,j)-y2(i-1,j))*h-JXKZ(j,1)*(rho_3*((s2(i-1,j)-y2(i-1,j))^2+(s2(i-1,j-tao)-y2(i-1,j-tao))^2)*h/2)^1*sign(s2(i-1,j)-y2(i-1,j))*h;
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
%mesh(t,x,y3),xlabel('t'),ylabel('x'),zlabel('y_3');
%axis([-3.5 5 -5 5 -0.8 0.8])
%subplot(2,3,1),mesh(t,x,y1)
%subplot(2,3,3),mesh(t,x,y2)
%subplot(2,3,5),mesh(t,x,y3)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % figure
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % mesh(t,x,y3),xlabel('t'),ylabel('x'),zlabel('y_3');
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % axis([-3.5 5 -5 5 -0.8 0.8])
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % subplot(1,3,1),mesh(t,x,y1)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % subplot(1,3,2),mesh(t,x,y2)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % subplot(1,3,3),mesh(t,x,y3)
% surf(x,t,z1);



%figure
%mesh(t,x,s3),xlabel('t'),ylabel('x'),zlabel('S_3');
%axis([-3.5 5 -5 5 -0.8 0.8])
%subplot(2,3,1),mesh(t,x,s1)
%subplot(2,3,3),mesh(t,x,s2)
%subplot(2,3,5),mesh(t,x,s3)






% % % % % % % figure
% % % % % % % mesh(t,x,s3),xlabel('t'),ylabel('x'),zlabel('S_3');
% % % % % % % axis([-3.5 5 -5 5 -0.8 0.8])
% % % % % % % subplot(1,3,1),mesh(t,x,s1)
% % % % % % % subplot(1,3,2),mesh(t,x,s2)
% % % % % % % subplot(1,3,3),mesh(t,x,s3)
% surf(x,t,z1);



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
%subplot(2,2,3),mesh(t,x,s7-y7),xlabel('t'),ylabel('s'),zlabel('v_7')
%subplot(2,2,4),mesh(t,x,s8-y8),xlabel('t'),ylabel('s'),zlabel('v_8')
%figure
%subplot(2,2,1),mesh(t,x,s9-y9),xlabel('t'),ylabel('s'),zlabel('v_9')
%subplot(2,2,2),mesh(t,x,s10-y10),xlabel('t'),ylabel('s'),zlabel('v_{10}')
subplot(2,2,3),mesh(t,x,u1),xlabel('t'),ylabel('s'),zlabel('w_1')
subplot(2,2,4),mesh(t,x,u2),xlabel('t'),ylabel('s'),zlabel('w_2')
%figure
%subplot(2,2,1),mesh(t,x,u7),xlabel('t'),ylabel('s'),zlabel('e_7')
%subplot(2,2,2),mesh(t,x,u8),xlabel('t'),ylabel('s'),zlabel('e_8')
%subplot(2,2,3),mesh(t,x,u9),xlabel('t'),ylabel('s'),zlabel('e_9')
%subplot(2,2,4),mesh(t,x,u10),xlabel('t'),ylabel('s'),zlabel('e_{10}')


