function eta = calcu_eta(h,ep)
	%用于已知划分求解每一小区间上的eta值
	n = length(h);
	x(1)=0;
	for m = 1:n	
		x(m+1)=x(m)+h(m);
	end	
	eta = [];
	for m =1:n
		%数值积分采用梯形公式
		f1 = 1/ep^2.*(exp(-x(m)./ep)./(1-exp(-1/ep))^2);
		f2 = 1/ep^2.*(exp(-x(m+1)./ep)./(1-exp(-1/ep))^2);
		eta(m)=h(m)*(f1+f2)/2;
		eta(m)=h(m)*(eta(m)^(1/2));
	end
end
