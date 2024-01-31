for h = [0.1,0.01,0.001,0.0001]
    x=-1:h:1;
    tau=h;
    for T=[0.2,0.4,0.6,0.8,1.0]
        u=[linspace(1,1,(length(x)-1)/2),linspace(0,0,(length(x)+1)/2)];
        for m = 1:T/tau
            %求sgn(a)，由于U>=0总成立，只需判断a是否等于0
            for j = 1:length(u)-1
                if u(j)==u(j+1)
                    a(j)=0;
                else
                    a(j)=1;
                end
            end
            u(2:length(u)-1)=u(2:length(u)-1)...
                -1/4*((1.+a(2:length(a))).*(u(2:length(u)-1).^2) ...
                +(1.-a(2:length(a))).*(u(3:length(u)).^2) ...
                -(1.+a(1:length(a)-1)).*(u(1:length(u)-2).^2) ...
                -(1-a(1:length(a)-1)).*(u(2:length(u)-1).^2));
            u(1)=1;%加左边值条件
            u(length(u))=0;%加右边值条件
        end
        plot(x,u);
        saveas(gcf,strcat("h=",num2str(h),",T=",num2str(T),".jpg"));
        clf;
    end
end
