function u = line_eq(h,ep)
	%h:该划分每两节点间的间距
	
	%a_u/a_m/a_d:矩阵A上对角线/对角线/下对角线的值
	a_u=-2./h(2:length(h)-1);
	a_m=2./h(1:length(h)-1)+2./h(2:length(h));
	a_d=-2./h(2:length(h)-1);
	
	%求节点的横坐标
	n =length(h);
	x(1)=0;
	for m = 1:n
		x(m+1)=x(m)+h(m);
	end
	x = x(2:n);
	%右端项
	u = 1/ep^2.*(exp(-x./ep)./(1-exp(-1/ep))^2).*(h(1:length(h)-1)+h(2:length(h)));
    u(length(u))=u(length(u))+2/h(length(h));
	
	%追赶法解三对角方程
	a_u(1) = a_u(1)/a_m(1); 
	for k = 2:n-2
		a_m(k) = a_m(k)-a_d(k-1)*a_u(k-1);
		a_u(k) = a_u(k)/a_m(k);
	end
	a_m(n-1) = a_m(n-1)-a_d(n-2)*a_u(n-2);
	u(1) = u(1)/a_m(1);
	for i = 2:n-1
		u(i) = (u(i)-a_d(i-1)*u(i-1))/a_m(i);
	end
	for j = n-2:-1:1
		u(j) = u(j)-a_u(j)*u(j+1);
	end
end
