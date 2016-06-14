function [] = plotxy(x,y,a,desc)
	%figure(a)
	plot(x,y,'Linewidth',3)
	xlabel('x')
	ylabel('y')
	legend(desc)
end