for theta=[0,1/2-1/(12*6/5),1/2,1]
    iter_box=[];%用于存储运算次数
    e_theta_box=[];%用于存储误差
    for n=3:7
        N=2^n;%计算N
        x=linspace(0,1,N+1);%打均匀网格
        U=[linspace(0,0,floor(1/3*N+1)),linspace(1,1,floor(2/3*N+1)-floor(1/3*N+1)),linspace(0,0,N+1-floor(2/3*N+1))];%初值条件
        h=1/N;%空间步长

        k=1:1:30;
        kp=k.*pi;
        m1=2*(cos(kp/3)-cos(2*kp/3))./(kp);
        m2=sin(kp'*x);
        m31=exp(-(kp).^2*0.1);
        m32=exp(-(kp).^2*1);
        m33=exp(-(kp).^2*10);
        u1=m1.*m31*m2;
        u2=m1.*m32*m2;
        u3=m1.*m33*m2;
        sol=[u1;u2;u3];

        %确定时间步长
        if theta==1/2-1/(12*6/5)
            tau=5/6*h^2;
        elseif theta==1/2
            tau=h^2;
        elseif theta==1
            tau=h;
        else
            tau=1/2*h^2;
        end

        iter_t=0;

        iter_box_temp=[];%用于存储N下运算次数
        e_theta_box_temp=[];%用于存储N下运算误差

        N_tem=0;%运算到第几个时间层
        j=0;
        for T=[0.1,1,10]
            j=j+1;

            %构建方程组
            diag_ele=1/tau+2*theta/h^2;%系数矩阵的对角元
            iter_t=iter_t+5;%更新iter_t
            diag_ele_vec=linspace(diag_ele,diag_ele,N-1);%系数矩阵的对角元（排成向量）
            while tau*N_tem<T

                F=(1-theta)*(U(1:N-1)-2*U(2:N)+U(3:N+1))/h^2+(1/tau)*U(2:N);%线性方程组右端项

                if theta==0
                    % F=F./diag_ele_vec;
                    A=diag(diag_ele_vec);
                    F=A\F';
                    F=F';
                    iter_t=iter_t+N-1;
                else
                    sub_diag_ele=-theta/h^2;%系数矩阵的次对角元
                    iter_t=iter_t+2;%更新iter_t
                    sub_diag_ele_vec=linspace(sub_diag_ele,sub_diag_ele,N-2);%系数矩阵的次对角元（排成向量）
                    A=diag(diag_ele_vec)+diag(sub_diag_ele_vec,1)+diag(sub_diag_ele_vec,-1);%系数矩阵
                    
                    %追赶法求解方程组
                    %带状LU分解
                    for i = 1:N-2
                        A(i+1,i)=A(i+1,i)/A(i,i);
                        A(i+1,i+1)=A(i+1,i+1)-A(i,i+1)*A(i,i+1)/A(i,i);
                        iter_t=iter_t+4;%更新iter_t
                    end
                    %解下三角方程组
                    for k = 2:N-1
                        F(k:N-1)=F(k:N-1)-A(k:N-1,k-1)'*F(k-1);
                        iter_t=iter_t+2*(N-k);%更新iter_t
                    end
                    %解上三角方程组
                    for k = N-2:-1:1
                        F(k+1)=F(k+1)/A(k+1,k+1);
                        F(1:k)=F(1:k)-A(1:k,k+1)'*F(k+1);
                        iter_t=iter_t+2*k+1;%更新iter_t
                    end
                    F(1)=F(1)/A(1,1);
                    iter_t=iter_t+1;%更新iter_t
                end
                
                U=[0,F,0];%加边界条件
                N_tem=N_tem+1;%更新迭代步数
            end
            e_theta=max(abs(U-sol(j,:)));
            iter_box_temp=[iter_box_temp;iter_t];
            e_theta_box_temp=[e_theta_box_temp;e_theta];

        end

        iter_box=[iter_box,iter_box_temp];
        e_theta_box=[e_theta_box,e_theta_box_temp];

    end

    %处理数据并作图
    XN=[8,16,32,64,128];
    loglog(XN,e_theta_box(1,:));
    hold on;
    loglog(XN,e_theta_box(2,:));
    hold on;
    loglog(XN,e_theta_box(3,:));
    legend('T=0.1','T=1','T=10','Location','east');
    title(strcat('e_{\theta}-N'));
    saveas(gcf,strcat('theta=',num2str(theta),',N_e,no.jpg'));
    clf;
    loglog(iter_box(1,:),e_theta_box(1,:));
    hold on;
    loglog(iter_box(2,:),e_theta_box(2,:));
    hold on;
    loglog(iter_box(3,:),e_theta_box(3,:));
    legend('T=0.1','T=1','T=10','Location','east');
    title(strcat('e_{\theta}-运算次数'));
    saveas(gcf,strcat('theta=',num2str(theta),',iter_e,no.jpg'));
    clf;
end
