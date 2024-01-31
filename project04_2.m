%守恒型Lax-Wendroff格式
for h = [0.1,0.01,0.001,0.0001]
    x=-1:h:1;
    tau=h;
    for T=[0.2,0.4,0.6,0.8,1.0]
        u=[linspace(1,1,(length(x)-1)/2),linspace(0,0,(length(x)+1)/2)];%初值
        for m = 1:T/tau
            a=u(1:length(u)-1)+u(2:length(u));
            u(2:length(u)-1)=u(2:length(u)-1)-1/4*(u(3:length(u)).^2-u(1:length(u)-2).^2)+1/4*(a(2:length(a)).*(u(3:length(u)).^2-u(2:length(u)-1).^2)-a(1:length(a)-1).*(u(2:length(u)-1).^2-u(1:length(u)-2).^2));
            u(1)=1;%加左边值条件
            u(length(u))=0;%加右边值条件
        end
        plot(x,u);
        saveas(gcf,strcat("h=",num2str(h),",T=",num2str(T),"2.jpg"));
        clf;
    end
end