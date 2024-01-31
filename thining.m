function h = thining(h,eta,alpha)
	%用于对符合条件的区间进行二分加密
	n = length(h);
	em = alpha*max(eta);
	k=0;
	for m = 1:n
		if eta(m)>em
			h = [h(1:m-1+k),h(m+k)/2,h(m+k)/2,h(m+1+k:n+k)];
			k=k+1;
		end
	end
end