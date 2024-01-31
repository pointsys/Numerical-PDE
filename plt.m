function x = plt(u,h)
	%u:除边值外节点处有限元解值
	%h:每两节点间间距
	y = [0,u,1];
	x(1)=0;
	for m = 1:length(h)
		x(m+1)=x(m)+h(m);
	end
	plot(x,y,'-*',x,linspace(0,0,length(x)),'-*');
end