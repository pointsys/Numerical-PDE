n=0;%用于记录是第几次循环的数据，仅为方便区分图片文件名
for ep = [1/10,1/20,1/40,1/80,1/160] %epsilon的不同取值
    n=n+1;
    h = linspace(0.1,0.1,10);%首先对网格进行均分
    while length(h)<20
	    eta =calcu_eta(h,ep);%计算eta
	    h = thining(h,eta,0.5);%对网格进行加细
    end
    u = line_eq(h,ep);%用有限差分法求得数值解u
    plt(u,h);
    hold on;
    x=0:0.001:1;
    y=(1-exp(-x./ep))./(1-exp(-1/ep));%计算真解
    plot(x,y);%做出真解图像
    legend('数值解','网格','真解','Location','southeast');
    title('epsilon=',num2str(ep));
    file_name=['answer',num2str(n)];
    saveas(gcf,file_name,'jpg');
    clf;
end
