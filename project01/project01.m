%划分网格
Nr=50;%在r方向打51个点
Ntheta=100;%在theta方向打101个点
r=linspace(0,1,Nr+1);%网格点的r坐标
r=r(2:Nr);%内部网格点的r坐标
hr=1/Nr;%r方向的步长
htheta=2*pi/Ntheta;%theta方向的步长
theta=linspace(0,2*pi,Ntheta+1);%网格点的theta坐标
theta=theta(1:Ntheta);%除去theta=2*pi的点后的网格坐标

%计算方程右端项F
F=[];
for k = 1:Ntheta
    F=[F,0,4*exp(-8^2*((r.*cos(theta(k))).^2+(r.*sin(theta(k))-0.6).^2))];
end

%构建矩阵A，其中记号与上机报告中保持一致
a=-Nr^2.+1./(2.*r.*hr);%差分格式中U_{i-1,j}前的系数
b=-Nr^2.-1./(2.*r(1:Nr-2).*hr);%差分格式中U_{i+1,j}前的系数
c=2*Nr^2.+2./(r.^2.*htheta^2);%差分格式中U_{i,j}前的系数
T0=diag([0,c])+diag(a,-1)+diag([0,b],1);
T1=diag([1,c])+diag(a,-1)+diag([0,b],1);
d=-1./((r.^2).*(htheta^2));%差分格式中U_{i,j+1}和U_{i,j-1}前的系数
K1=diag([1,d]);
K_1=diag([-1,d]);
K0=diag([0,d]);
M=diag([-1,d]);
M(1,2)=1;
A=blkdiag(T1,kron(diag(linspace(1,1,Ntheta-1)),T0));%A的对角块
A=A+kron(diag((linspace(1,1,Ntheta-1)),1),K_1);
A=A+kron(diag((linspace(1,1,Ntheta-1)),-1),K1);
A(Nr*(Ntheta-1)+1,Nr*(Ntheta-2)+1)=0;
A(Nr*(Ntheta-1)+1:Nr*Ntheta,1:Nr)=M;
A(1:Nr,Nr*(Ntheta-1)+1:Nr*Ntheta)=K0;

%求解并画图
U=A\F';
U=U';
%画等高线图
W=[];
for k = 1:Ntheta
    W=[W;U(Nr*(k-1)+2:Nr*k),0];
end
W=W';
r=[r,1];
[TH,R]=meshgrid(theta,r);
[X,Y]=pol2cart(TH,R);
p=4*exp(-8^2*((R.*cos(TH)).^2+(R.*sin(TH)-0.6).^2));
contour(X,Y,p);

saveas(gcf,'load_contour.jpg');
contour(X,Y,W);

saveas(gcf,'deflection_contour.jpg');
%画y轴上的图
y=[fliplr(-r),0,r];
u_y=50.*[0,fliplr(U(Nr*(3*Ntheta/4)+2:Nr*(3*Ntheta/4+1))),U(1),U(Nr*(Ntheta/4)+2:Nr*(Ntheta/4+1)),0];
load_y=[fliplr(4*exp(-8^2*(-r-0.6).^2)),4*exp(-8^2*0.6.^2),4*exp(-8^2*(r-0.6).^2)];
plot(y,u_y);
hold on;
plot(y,load_y);
legend('Deflection','Load');
saveas(gcf,'yaxsis.jpg')












