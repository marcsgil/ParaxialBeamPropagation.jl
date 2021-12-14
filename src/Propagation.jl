function get_dft_interval(max_val::Float64,N::Int64)
    #Returns the interval over which we obtain the Fourier transform
	(-N÷2:N-1-N÷2)*(2*max_val/N)
end

function get_reciprocal_dft_interval(max_val::Float64,N::Int64)
    #Returns the interval over which we have to sample our function to perform the DFT
	(-N÷2:N-1-N÷2)*(π/max_val)
end

function get_fourier_transform(k_max::Float64,N::Int64,ψ0)
	#Calculates the Fourier transform of the given initial profile ψ0

    ks = get_dft_interval(k_max,N)
    rs = get_reciprocal_dft_interval(k_max,N)

	#This is a matrix containing the ψ0s
	ψ0s = convert(Matrix{ComplexF64},[ ψ0(x,y) for x in rs, y in rs])

	#Then, we calculate the Fourier transform of ψ, which returns a matrix
	ϕs =  (0.5*π/k_max^2)*ifftshift( fft( fftshift(  ifftshift(fft(fftshift(ψ0s,1),1),1)  ,2),2),2)
    
	#Now, we interpolate this matrix, to be able to calculate the function at an arbitrary point.
	CubicSplineInterpolation( (ks , ks),ϕs,extrapolation_bc=0.0)
end

function display_fourier_transform(plot_range::Float64,plot_N::Int64,clim::Float64,ϕ)
	#Displays the fourier transform ϕ. It is useful to see if the chosen parameters are adequate
	ks = LinRange(-plot_range,plot_range,plot_N)

	I = [ real(ϕ(kx,ky)*conj(ϕ(kx,ky))) for ky in ks, kx in ks ]

	heatmap(ks,ks,I,clims=(0,clim),aspect_ratio=:equal,
	xlabel=L"k_x",ylabel=L"k_y",xlims=(-plot_range,plot_range),ylims=(-plot_range,plot_range),
	title="Fourier transform of initial profile")
end

function propagate_beam(r_max::Float64,N::Int64,z::Float64,ϕ)
	#Propagates a beam to a distance z given the fourier transform of its initial profile

	rs = get_dft_interval(r_max,N)
    ks = get_reciprocal_dft_interval(r_max,N)
			
	s = zeros(ComplexF64,N,N)

	Threads.@threads for n in 1:N
		Threads.@threads for m in 1:N
			s[n,m] = ϕ(ks[m],ks[n])*exp(-im*z*(ks[n]^2+ks[m]^2)/2)
		end
	end

	( 0.5*π/r_max^2 )*ifftshift( bfft( fftshift(  ifftshift(bfft(fftshift(s,1),1),1)  ,2),2),2)
end

function propagate_beam(r_max::Float64,N::Int64,z::Float64,ϕ,bplan1,bplan2)
	#Same as before, but accepts plans for the FFT.

	rs = get_dft_interval(r_max,N)
    ks = get_reciprocal_dft_interval(r_max,N)
			
	s = zeros(ComplexF64,N,N)
	
	Threads.@threads for n in 1:N
		Threads.@threads for m in 1:N
			s[n,m] = ϕ(ks[m],ks[n])*exp(-im*z*(ks[n]^2+ks[m]^2)/2)
		end
	end

	( 0.5*π/r_max^2 )*ifftshift(bplan2*fftshift(ifftshift(bplan1*fftshift(s,1),1),2),2)
end

function decrease_size(a,rs,I)
	#This function decreses the size of a range and corresponding intensities. 
	#It is analogous to setting xlims=(-a,a), ylims=(-a,a) but somehow is faster
	K=1
	for k in eachindex(rs)
		if rs[k]>-a
			K = k
			break
		end
	end
	rs[ K:(length(rs)+1-K) ],I[K:(length(rs)+1-K),K:(length(rs)+1-K)]
end;

function display_beam(r_max::Float64,N::Int64,z::Float64,plot_range::Float64,clim::Float64,ϕ)
	#Displays a propageted beam
	rs = get_dft_interval(r_max,N)
	ψs = propagate_beam(r_max,N,ϕ,z)

	I = real(ψs.*conj(ψs))

	drs,dI = decrease_size(plot_range,rs,I)

	heatmap(drs,drs,dI,clims=(0,clim),aspect_ratio=:equal,title="z="*(@sprintf "%.2f" z),xlabel=L"x",ylabel=L"y",xlims=(-plot_range,plot_range),ylims=(-plot_range,plot_range))
end

function animate_beam(r_max::Float64,N::Int64,zmax::Float64,
	plot_range::Float64,clim::Float64,ϕ,fps::Int64,nframes::Int64,file_name::String)
	#Animates a beam propagation
	#The file_name should end in .gif
	
	rs = get_dft_interval(r_max,N)
	ks = get_reciprocal_dft_interval(r_max,N)

	plan1 = plan_bfft( zeros(ComplexF64, N,N ) , 1)
	plan2 = plan_bfft( zeros(ComplexF64, N,N ) , 2)

	anim = @animate for z in LinRange(0,zmax,nframes)
		ψs = propagate_beam(r_max,N,ϕ,z,plan1,plan2)

		I = real(ψs.*conj(ψs))

		drs,dI = decrease_size(plot_range,rs,I)

		heatmap(drs,drs,dI,clims=(0,clim),aspect_ratio=:equal,
		title=L"z="*(@sprintf "%.2f" z),xlabel=L"x",ylabel=L"y",
		xlims=(-plot_range,plot_range),size=(1200,900),labelfontsize=20,titlefontsize=20,tickfontsize=14)
	end

	gif(anim,file_name,fps=fps)
end;