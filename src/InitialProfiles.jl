#Hermite-Gaussian Modes
HG(n,m,x,y) = hermiteh(n,√2*x)*hermiteh(m,√2*y)*exp(-x^2-y^2);

function θ(x,y)
	#θ ∈ [-π,π]
	if x ≥ 0
		atan(y/x)
	else
		if y≥ 0
			atan(y/x) + π
		else
			atan(y/x) - π
		end
	end
end;

function LG(p,l,x,y)
	#Laguerre-Gaussian Modes
    if x==0 && y==0
        0.0
    else
        (2*(x^2+y^2))^(0.5*abs(l))*laguerrel(p,abs(l),2*(x^2+y^2))*exp(-(x^2+y^2))*exp(im*l*θ(x,y))
    end
end

function airy_beam(x,α)
	#Airy beam (in a single direction)
	airyai(x)*exp(α*x)
end

#Series of obstacles to simulate difraction
function vertical_obstacle(x,y,l)
	if abs(x)<l
		0.0
	else
		1.0
	end
end

function horizontal_obstacle(x,y,l)
	if abs(y)<l
		0.0
	else
		1.0
	end
end

function circular_obstacle(x,y,r)
	if x^2+y^2<r
		0.0
	else
		1.0
	end
end

function rectangular_obstacle(x,y,a,b,c,d)
	if (a<x<b) & (c<y<d)
		0.0
	else
		1.0
	end
end