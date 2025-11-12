module Math
	using ForwardDiff #package used for chain rules of large derivatives
    using LinearAlgebra
	using Memoization
    using Base.Threads

    # Constants
    const c::Float64 = 2.998e5    # km/s
    const G::Float64 = 1.327e11   # km^3 / M⊙ / s^2
    const period_factor::Float64 = sqrt(4.0 * π^2 / G)
    const dJ_rr_factor::Float64 = -2.0 * G / (5.0 * c^5)
    const Mpc_to_km_factor = 2.0*G*3.2e-20/c^4. #for conversion from Mpc to km
    const M_ch=1.44 #Chandresekhar mass in M⊙
    
    """
        Returns the natural radius of a white dwarf.
        ```julia
        R0_WD(M_WD,mu_e=2.0)
        ```
        M_WD = mass of the white dwarf.
        mu_e = mean molecular mass per electron (default 2 for white dwarfs of all but the most exotic chemical compositions (e.g., ONeMG)).
    """
    @memoize function R0_WD(M_WD;mu_e=2.0,equation=1) 
    	#@assert M_WD > 0 # In production ready code, consider commenting out assertion
    	#@assert mu_e > 0 # In production ready code, consider commenting out assertion
    	#white dwarf radius in km as ftn of mass in M_odot and mean molecular weight/electron=2
        
        x=M_WD/M_ch
        y=M_WD/0.00057
        #see Even (2009; https://arxiv.org/pdf/0908.2116) for summary
        if equation==1 #Nauenberg equation
        	return(15584.0/mu_e * x^(-1.0/3.0) * (1-x^(4.0/3.0))^(1.0/2.0))
        else #Eggleton equation
            return(7931.0*(x^(-2.0/3.0)-x^(2.0/3.0))^(1.0/2.0) * (1.0+3.5*y^(-2.0/3.0)+y^(-1.0))^(-2.0/3.0)) 
        end
    end
    
    """
        Returns the dynamical timescale of the white dwarf.
        ```julia
        dynamical_timescale(M_WD,R_WD)
        ```
        M_WD = mass of the white dwarf.
        R_WD = radius of the white dwarf.
    """
    @memoize function dynamical_timescale(M_WD,R_WD)
    	#@assert M_WD > 0 # In production ready code, consider commenting out assertion
    	#@assert R_WD > 0 # In production ready code, consider commenting out assertion
    	#returns dynamical timescale for perturbed white dwarf
    	#with radius R and mass M
    	return((R_WD^3.0/(G*M_WD))^(1.0/2.0))
    end
    
    """
        Returns the Roche limit of the binary, using the Eggleton approximation.
        ```julia
        Roche_Limit(a,M_WD,M_NS)
        ```
        a = total orbital separation of the binary.
        M_WD = mass of the white dwarf.
        M_NS = mass of the neutron star.
    """
    @memoize function Roche_Limit(a,M_WD,M_NS)
    	#@assert M_WD < M_NS # In production ready code, consider commenting out assertion
    	@assert a > 0 # In production ready code, consider commenting out assertion
    	#returns Roche limit of M_WD for companion of mass M_NS 
    	#at a separation of distance r 
    	q=M_WD/M_NS
        q_2_3rds = q^(2.0/3.0)
    	return(a * 0.49*q_2_3rds/(0.6*q_2_3rds+log(1.0+q^(1.0/3.0))))
    end
    
    
    """
        Returns the orbital period of the binary.
        ```julia
        period(a, M_tot)
        ```
        a = total orbital separation of the binary.
        M_tot = total mass of the binary (white dwarf + neutron star).
    """
    @memoize function period(a, M_tot)
        @assert a > 0 # In production ready code, consider commenting out assertion
    	@assert M_tot > 0 # In production ready code, consider commenting out assertion
        #return sqrt(4.0 * π^2 * a^3 / (G * Mtot))
        return period_factor * sqrt(a^3.0 / M_tot)
    end
    
    
    """
        Returns the decretion (mass loss) rate of the white dwarf.
        ```julia
        decretion_rate(a, R_WD, M_WD, M_NS; A=10.0)
        ```
        a = total orbital separation of the binary.
        R_WD = radius of the white dwarf.
        M_WD = mass of the white dwarf.
        M_NS = mass of the neutron star.
        A = numerical constant determined from fluid mechanics (default 10).
        If Roche lobe is greater than radius of white dwarf, returns 0.
    """
    @memoize function decretion_rate(a, R_WD,
                            M_WD, M_NS; A=10.0)
    	#@assert R_WD > 0 # In production ready code, consider commenting out assertion
        RL = Roche_Limit(a, M_WD, M_NS)
        if RL > R_WD
            return 0.0
        end
        return -A * M_WD / period(a, M_WD + M_NS) * ((R_WD - RL) / R_WD)^3.0
    end
    
    """
        Returns the specific angular momentum of material ejected from binary.
        ```julia
        gamma(M_WD, M_NS; a_ring_frac=1.0, mode::String="Isotropic Re-emission") 
        ```
        M_WD = mass of the white dwarf.
        M_NS = mass of the neutron star.
        a_ring_frac = fractional radius (in units of the binary’s total orbital separation) at which material is forms a circumbinary ring in the “Circumbinary Ring” mass loss mode. Default: 1.0.
        mode = mass loss mode. One of "Jeans mode" ('j'), "Isotropic Re-emission" ('i'), or "Circumbinary Ring" ('c'). Default: "Isotropic Re-emission")
    
        Jeans mode - material leaves the donor star directly, carrying that star’s specific angular momentum.
        Isotropic re-emission - material first accretes onto the neutron star, then is ejected isotropically, carrying the accretor’s angular momentum.
        Circumbinary ring - material forms a ring around both stars at some radius a_ring = a x a_ring_frac, where a is the total binary separation.
    """
    @memoize function gamma(M_WD, M_NS; 
                   a_ring_frac=1.0, mode::String="Isotropic Re-emission") 
    	#@assert M_WD > 0 # In production ready code, consider commenting out assertion
    	@assert M_NS > 0 # In production ready code, consider commenting out assertion
    	@assert a_ring_frac > 0 # In production ready code, consider commenting out assertion
        first_char = lowercase(first(mode))
        if first_char == 'j'
            return M_NS / M_WD
        elseif first_char == 'i'
            return M_WD / M_NS
        elseif first_char == 'c'
            return (M_WD + M_NS)^2.0 / (M_WD * M_NS) * sqrt(a_ring_frac)
        else
            error("$mode not defined. Use 'Jeans Mode', 'Isotropic Re-emission', or 'Circumbinary Ring'")
        end
    end
    
    
    """
        Returns the Eddington accretion rate onto an object.
        ```julia
        eddington_rate(M; eta_acc=0.1)
        ```
        M = mass of the accretor.
        eta_acc = accretion rate efficiency (default 0.1).
        Default values apply to Eddington accretion onto the neutron star from the white dwarf.
    """
    @memoize function eddington_rate(M; eta_acc=0.1)
        @assert M > 0 # In production ready code, consider commenting out assertion
    	@assert 0 < eta_acc <= 1 # In production ready code, consider commenting out assertion
        # 4πG*m_p / (σ_T * eta_acc * c) * M
        # = 7.05e-16 * M * (0.1 / eta_acc)
        return 7.05e-17 * M / eta_acc  # M in solar masses
    end
    
    
    """
        Returns the (fraction of) mass accreted onto the neutron star over the mass lost from the white dwarf per unit time.
    """
    @memoize function beta(a, R_WD, M_WD, M_NS; 
                  A=10.0, eta_acc=0.1)
        M_WD_dot = -decretion_rate(a, R_WD, M_WD, M_NS; A=A)
        lam_edd = eddington_rate(M_NS; eta_acc=eta_acc)
        if M_WD_dot < lam_edd
            return 1.0
        end
        return lam_edd / M_WD_dot
    end
    
    
    """
        Returns the rate of total angular momentum loss of the binary due to non-conservative mass transfer.
    """
    @memoize function J_dot(J, a, R_WD, M_WD, M_NS; 
                   A=10.0, eta_acc=0.1, 
                   a_ring_frac=1.0, mode::String="Isotropic Re-emission")
        g = gamma(M_WD, M_NS; a_ring_frac=a_ring_frac, mode=mode)
        b = beta(a, R_WD, M_WD, M_NS; A=A, eta_acc=eta_acc)
    	@assert b <= 1.0 # In production-ready code, consider commenting this out
        M_WD_dot = decretion_rate(a, R_WD, M_WD, M_NS; A=A)
        return J * g * (1.0 - b) * M_WD_dot / (M_WD + M_NS)
    end

    """
        Returns the rate of change for total orbital separation.
    """
    @memoize function a_dot(a, R_WD, M_WD, M_NS; 
                   A=10.0, eta_acc=0.1, 
                   a_ring_frac=1.0, mode::String="Isotropic") 
        # Promote keyword values to T so arithmetic mixes consistently
        g = gamma(M_WD, M_NS; a_ring_frac=a_ring_frac, mode=mode)
        b = beta(a, R_WD, M_WD, M_NS; A=A, eta_acc=eta_acc)
        @assert b <= 1.0 # In production-ready code, consider commenting this out
        M_WD_dot = decretion_rate(a, R_WD, M_WD, M_NS; A=A)
        one_minus_b = 1.0 - b
        return -2.0 * a * M_WD_dot / M_WD * (one_minus_b * M_WD / M_NS - one_minus_b * (g + 0.5) * M_WD / (M_WD + M_NS))
    end
    
    
    """
        Returns rate of change for white dwarf radius.
        ```julia
        dot_R_WD(M_WD, R_WD)
        ```
        M_WD = mass of white dwarf.
        R_WD = radius of white dwarf.
    """
    @memoize function dot_R_WD(M_WD, R_WD; mu_e=2.0)
        natural_R = R0_WD(M_WD, mu_e=mu_e)
        tau = dynamical_timescale(M_WD, R_WD)
        return (natural_R - R_WD) / tau
    end
    
    
    """
        Returns reduced mass.
        ```julia
        mu(M_WD, M_NS)
        ```
        M_WD = mass of white dwarf.
        M_NS = mass of neutron star.
    """
    @memoize function mu(M_WD, M_NS)
    	#@assert M_WD > 0 # In production-ready code, consider commenting this out
    	@assert M_NS > 0 # In production-ready code, consider commenting this out
        return M_WD * M_NS / (M_WD + M_NS)
    end
    
    """
        Returns the rate of total angular momentum loss of the binary due to the radiation reaction (see math notes).
        ```julia
        dJ_rr(I2::Matrix{Float64}, I3::Matrix{Float64})::Float64
        ```
        I2 = second order time derivative of quadrupole moment (quadrupole acceleration tensor).
        I3 = third order time derivative of quadrupole moment (quadrupole jerk tensor).
    """
    @memoize function dJ_rr(I2::Matrix{Float64}, I3::Matrix{Float64})::Float64
        # Corresponds to: -2G/(5c^5) * (dot(I2[0],I3[1]) - dot(I2[1],I3[0]))
        term::Float64 = dot(I2[:, 1], I3[:, 2]) - dot(I2[:, 2], I3[:, 1])
        #return -2.0 * G / (5.0 * c^5) * term
        return dJ_rr_factor * term
    end
    
    """
        Converts polar coordinates to Cartesian coordinates.
    """
    @memoize function x(r::T, phase::T)::T where {T<:Real}
        @assert r > 0 # In production-ready code, consider commenting this out
        return r * cos(phase)
    end
    
    """
        Converts polar coordinates to Cartesian coordinates.
    """
    @memoize function y(r::T, phase::T)::T where {T<:Real}
        @assert r > 0 # In production-ready code, consider commenting this out
        return r * sin(phase)
    end
    
    """ 
        Returns rate of change of neutron star mass.
    """
    @memoize function dM_NS(a::T, R_WD::T, M_WD::T, M_NS::T;
            A=10., eta_acc=0.1)::T where T<:Real
        t_dM_WD=decretion_rate(a,R_WD,M_WD,M_NS,A=A)
        
        t_edd=eddington_rate(M_NS, eta_acc=eta_acc)
        return(min(-t_dM_WD,t_edd))
    end
    
    """ 
        Returns rate of change of orbital phase of the white dwarf.
    """
    @memoize function dphase(J::T,a::T,M_WD::T, M_NS::T; mode::String="Debug")::T where T<:Real
        if mode=="Debug"
            @assert M_WD > 0 # In production-ready code, consider commenting this out
            @assert M_NS > 0 # In production-ready code, consider commenting this out
            @assert a > 0 # In production-ready code, consider commenting this out
        end
        return(J/(a^2.0*mu(M_WD,M_NS)))
    end   
    
    
    """
        Returns I_{xx}, where I is the quadrupole moment tensor.
        ```julia
        I_xx(M_WD::T, M_NS::T, r::T, phase::T; mode::String="Debug")
        ```
    """
    @memoize function I_xx(M_WD::T, M_NS::T, r::T, phase::T; mode::String="Debug")::T where T<:Real
        if mode=="Debug"
            @assert M_WD > 0 # In production-ready code, consider commenting this out
            @assert M_NS > 0 # In production-ready code, consider commenting this out
            @assert r > 0 # In production-ready code, consider commenting this out
        end
        return (mu(M_WD, M_NS) * (x(r,phase))^2.0)
    end           
    
    """Given polar orbital parameters, calculates I_{xy}, where I is the quadrupole moment tensor."""
    @memoize function I_xy(M_WD::T, M_NS::T, r::T, phase::T; mode::String="Debug")::T where T<:Real
        if mode=="Debug"
            @assert M_WD > 0 # In production-ready code, consider commenting this out
            @assert M_NS > 0 # In production-ready code, consider commenting this out
            @assert r > 0 # In production-ready code, consider commenting this out
        end
        return mu(M_WD, M_NS) * x(r,phase) * y(r,phase)
    end                
    
    """Fixes the higher derivatives of the function so that it is formatted in a way that ForwardDiff is comfortable with."""
    function fix_higher_derivatives(ftn::Function,higher_derivatives::AbstractVector{T})::Function where {T<:Real} 
    	function inner_ftn(gen_theta::AbstractVector{T}) where {T<:Real}
    		return(ftn(gen_theta,higher_derivatives))
    	end
    	return ftn 
    end;
    
    """Given a function, its vectorized parameters, and the time derivatives of those parameters, returns the time derivative of the function as calculated via the chain rule."""
    function time_derivative(ftn::Function, ind_vars_and_derivatives::AbstractVector{T}; mode::String="Debug") where {T<:Real}
        #ftn is the function to be taken a time derivative of
        #ind_vars is an array with the values of the arguments of ftn
        #derivatives is an array with the known time derivatives of ind_vars
    	#print(length(ind_vars_and_derivatives))
        if mode=="Debug"
            @assert mod(length(ind_vars_and_derivatives),6)==0
        end
    
    	ind_vars=@view ind_vars_and_derivatives[1:end-6]
    	derivatives=@view ind_vars_and_derivatives[7:length(ind_vars_and_derivatives)]
    
    	# ftn_at_ind_vars is the function evaluated at the parameters in ind_vars
        ftn_at_ind_vars = DiffResults.GradientResult(ind_vars)
    	ForwardDiff.gradient!(ftn_at_ind_vars,ftn,ind_vars)

        # DiffResults.grad(ftn_at_ind_vars) gives the gradient of the function, evaluated at the parameters in ind_vars
        return(sum(DiffResults.gradient(ftn_at_ind_vars) .* derivatives))#chain rule
    end
    
    """Given a function, defines its next derivative based on the function's parameters and their time derivatives. `n` determines how many time derivatives of \vec{phase} (i.e., the parameters to be integrated) are parameters of the function itself. """
    function next_derivative(in_ftn::Function; n::Integer=1, mode::String="Debug")::Function
        #returns a function which is the next higher derivative of the in_ftn
        #n is the total derivatives past the given analytic form
        if mode=="Debug"
        	@assert n>0
        end
        function ftn(gen_theta_and_derivatives::AbstractVector{T}) where {T<:Real}
            if mode=="Debug"
                @assert length(gen_theta_and_derivatives)==6*(n+1) 
            end
    		
            return(time_derivative(in_ftn, gen_theta_and_derivatives))
        end
        return ftn
    end
    
    #The following functions wrap the original derivatives in vectors to make life easier when taking derivatives due to mass transfer only. In particular, `theta` is of the form `[a, M_WD, M_NS, R_WD, J, phase]`.
    """Wrapper function for time derivative of total orbital separation."""
    function da_vector(theta::AbstractVector{T}; 
            A=10.0, eta_acc=0.1, 
            a_ring_frac=1.0, mode::String="Isotropic")::T where {T<:Real}
    
        # promote the numeric keyword args to match the type of the inputs
        A = T(A)
        eta_acc = T(eta_acc)
        a_ring_frac = T(a_ring_frac)
    
        #unpacking theta
        a, M_WD, M_NS, R_WD, J, phase=theta
    
        
        return(a_dot(a,R_WD, M_WD, M_NS;
                A=A,eta_acc=eta_acc,
                a_ring_frac=a_ring_frac,mode=mode))
    end
    
    """Wrapper function for time derivative of white dwarf mass."""
    function dM_WD_vector(theta::AbstractVector{T}; A=10.)::T where {T<:Real}
        A = T(A)
        
        a, M_WD, M_NS, R_WD, J, phase=theta
        return(decretion_rate(a, R_WD, M_WD, M_NS; A=A))
    end
    
    """Wrapper function for time derivative of neutron star mass."""
    function dM_NS_vector(theta::AbstractVector{T};
            A=10, eta_acc=0.1)::T where {T<:Real}
    
        A = T(A)
        eta_acc = T(eta_acc)
        
        a, M_WD, M_NS, R_WD, J, phase=theta
        return(dM_NS(a, R_WD, M_WD, M_NS;
            A=a, eta_acc=eta_acc))
    end
    
    """Wrapper function for time derivative of white dwarf radius."""
    function dR_WD_vector(theta::AbstractVector{T}; mu_e=2.)::T where {T<:Real} 
        mu_e = T(mu_e)
        
        a, M_WD, M_NS, R_WD, J, phase=theta
        return( dot_R_WD(M_WD, R_WD, mu_e=mu_e))
    end
    
    """Wrapper function for time derivative of total orbital momentum of the binary."""
    function dJ_vector(theta::AbstractVector{T}; 
                   A=10.0, eta_acc=0.1, 
                   a_ring_frac=1.0, mode::String="Isotropic")::T where {T<:Real}
        A = T(A)
        eta_acc = T(eta_acc)
        a_ring_frac = T(a_ring_frac)
        
        a, M_WD, M_NS, R_WD, J, phase=theta
        return(J_dot(J, a, R_WD, M_WD, M_NS, 
                   A=A, eta_acc=eta_acc, 
                   a_ring_frac=a_ring_frac, mode=mode))
    end
    
    """Wrapper function for time derivative of orbital phase of white dwarf."""
    function dphase_vector(theta::AbstractVector{T})::T where {T<:Real} 
        a, M_WD, M_NS, R_WD, J, phase=theta
        return(dphase(J,a,M_WD, M_NS))
    end
    
    """Wrapper function for quadrupole moment tensor's xx component."""
    function Ixx_vector(theta::AbstractVector{T})::T where {T<:Real} 
        a, M_WD, M_NS, R_WD, J, phase=theta
        return(I_xx(M_WD, M_NS, a, phase))
    end
    
    """Wrapper function for quadrupole moment tensor's xy component."""
    function Ixy_vector(theta::AbstractVector{T})::T  where {T<:Real} 
        a, M_WD, M_NS, R_WD, J, phase=theta
        return(I_xy(M_WD, M_NS, a, phase))
    end
    
    #"The below block of code uses the vectorized derivatives above to output functions which take higher time derivatives of $\vec{\theta}$ and $I_{ij}$ as a function of phase. "
    dda_vector=next_derivative(da_vector)
    ddM_WD_vector=next_derivative(dM_WD_vector)
    ddM_NS_vector=next_derivative(dM_NS_vector)
    ddR_WD_vector=next_derivative(dR_WD_vector)
    ddJ_vector=next_derivative(dJ_vector)
    ddphase_vector=next_derivative(dphase_vector)
    
    ddda_vector=next_derivative(dda_vector,n=2)
    dddM_WD_vector=next_derivative(ddM_WD_vector,n=2)
    dddM_NS_vector=next_derivative(ddM_NS_vector,n=2)
    dddR_WD_vector=next_derivative(ddR_WD_vector,n=2)
    dddJ_vector=next_derivative(ddJ_vector,n=2)
    dddphase_vector=next_derivative(ddphase_vector,n=2)
    
    dI_xx_vector=next_derivative(Ixx_vector)
    ddI_xx_vector=next_derivative(Ixx_vector,n=2)
    dddI_xx_vector=next_derivative(Ixx_vector,n=3)
    
    dI_xy_vector=next_derivative(Ixy_vector)
    ddI_xy_vector=next_derivative(Ixy_vector,n=2)
    dddI_xy_vector=next_derivative(Ixy_vector,n=3)
    
    #=bundling derivatives together which will prove convenient later=#
    first_derivatives=[da_vector, dM_WD_vector, dM_NS_vector, dR_WD_vector, dJ_vector, dphase_vector]
    second_derivatives=[dda_vector, ddM_WD_vector, ddM_NS_vector, ddR_WD_vector, ddJ_vector, ddphase_vector]
    third_derivatives=[ddda_vector, dddM_WD_vector, dddM_NS_vector, dddR_WD_vector, dddJ_vector, dddphase_vector]

    """Constructs a symmetric 2x2 gravitational-wave strain tensor corresponding to given “plus” (+) and “cross” (×) polarization amplitudes."""
    function tensor_plus_cross(plus::Float64, cross::Float64)::Matrix{Float64}
        return([plus cross; cross -plus        
        ])
    end
    
    #"The below block of code is central to the interface between the mass transfer and radiation reaction effects. It also sits in the innermost loop of our integrator. This is a good place to peer review."
    #=Central Engine of the integrator=#
    """Calculates the second derivative of the quadrupole moment `t_ddI` due to kinematics of keplerian motion and mass transfer `dtheta`. Also returns higher derivatives of theta useful for calculating the third derivative of the quadrupole moment."""
    function ddI_helper(theta::AbstractVector{T})::Tuple{Vector{T},Vector{T},Vector{T},Array{T, 2}} where {T}
        #calulates a list of values which will be helpful in the next two functions
        dtheta=zeros(6)
        #calculating first time derivatives due to MT only
        
        Threads.@threads for i in 1:6
           @inbounds dtheta[i]=first_derivatives[i](theta)
        end
    
        theta_dtheta=vcat(theta,dtheta)
        ddtheta=zeros(6)
        #calculating second time derivatives due to MT only
        Threads.@threads for i in 1:6
            @inbounds ddtheta[i]=second_derivatives[i](theta_dtheta)
        end
    
        theta_ddtheta=vcat(theta_dtheta,ddtheta)
    
        #using orbital elements and first and second time derivatives 
        #to calculate second derivative of quadrupole moment
        #(for use in calculating RR acceleration)
        t_ddI=tensor_plus_cross(ddI_xx_vector(theta_ddtheta),
            ddI_xy_vector(theta_ddtheta)
        )
        return(dtheta, ddtheta, theta_ddtheta, t_ddI)
    end


    """The first derivative is much costlier than the following derivatives (about as costly as all the others put together), so this function parallelizes using only two workers by having the first worker calculate the first derivative while the second worker calculates the second through sixth derivatives."""
    function uneven_parallelization_helper(theta::AbstractVector{T}, derivatives)::Vector{T} where {T}
        dtheta = zeros(6)

        # Task 1: expensive separation derivative (separation--many chain rules) 
        task_1 = Threads.@spawn begin
            dtheta[1] = derivatives[1](theta)
        end
    
        # Task 2: lighter subsequent derivatives
        task_2 = Threads.@spawn begin
            @inbounds for j in 2:6
                dtheta[j] = derivatives[j](theta)
            end
        end
    
        # Wait for completion
        fetch(task_1)
        fetch(task_2)
    
        return dtheta
end
        

    """Same as above function, but attempting to parallelize by dedicating one worker to more intensive calculation of derivatives of separation (a) and one worker to calculating all other derivatives """
    function ddI_helper_uneven(theta::AbstractVector{T})::Tuple{Vector{T},Vector{T},Vector{T},Array{T, 2}} where {T}
        #calulates a list of values which will be helpful in the next two functions
        dtheta=uneven_parallelization_helper(theta, first_derivatives)
        #calculating first time derivatives due to MT only
    
    
        theta_dtheta=vcat(theta,dtheta)
        ddtheta=uneven_parallelization_helper(theta_dtheta, second_derivatives)
        #calculating second time derivatives due to MT only
        
    
        theta_ddtheta=vcat(theta_dtheta,ddtheta)
    
        #using orbital elements and first and second time derivatives 
        #to calculate second derivative of quadrupole moment
        #(for use in calculating RR acceleration)
        t_ddI=tensor_plus_cross(ddI_xx_vector(theta_ddtheta),
            ddI_xy_vector(theta_ddtheta)
        )
        return(dtheta, ddtheta, theta_ddtheta, t_ddI)
    end


    """Returns the second derivative of the quadrupole moment as calculated in `ddI_Helper`. """
    function ddI_from_theta(theta::Vector{Float64})::Array{Float64,2}
        return(ddI_helper(theta)[end])
    end
    	    
    """Calculates the third derivative of the quadrupole moment and returns it along with the time derivative of theta and the second derivative of the quadrupole moment as calculated in `ddI_Helper`."""
    function ddI_dddI_from_theta(theta::Vector{Float64};ddI_only::Bool=false)::Tuple{Vector{Float64},Array{Float64,2},Array{Float64,2}}
        #returns time derivatives of theta and second and third derivatives of 
        #the quadrupole moment due to MT effects
        
        dtheta, ddtheta, theta_ddtheta, t_ddI = ddI_helper(theta)
        #using orbital elements and first through third time derivatives  
        #to calculate third derivative of quadrupole moment 
        #(for use in calculating RR acceleration)
        
        dddtheta=zeros(6)
        #calculating second time derivatives due to MT only
        Threads.@threads for i in 1:6
            @inbounds dddtheta[i]=third_derivatives[i](theta_ddtheta)
        end
    
        theta_dddtheta=vcat(theta_ddtheta,dddtheta)
    
        t_dddI=tensor_plus_cross(dddI_xx_vector(theta_dddtheta),
            dddI_xy_vector(theta_dddtheta)
        )
        return(dtheta,t_ddI,t_dddI)
    end


"""Uneven parallelization to calculate third derivatives """
    function ddI_dddI_from_theta_uneven(theta::Vector{Float64};ddI_only::Bool=false)::Tuple{Vector{Float64},Array{Float64,2},Array{Float64,2}}
        #returns time derivatives of theta and second and third derivatives of 
        #the quadrupole moment due to MT effects
        
        dtheta, ddtheta, theta_ddtheta, t_ddI = ddI_helper_uneven(theta)
        #using orbital elements and first through third time derivatives  
        #to calculate third derivative of quadrupole moment 
        #(for use in calculating RR acceleration)
        
        dddtheta=uneven_parallelization_helper(theta_ddtheta, third_derivatives)
        #calculating second time derivatives due to MT only
    
        theta_dddtheta=vcat(theta_ddtheta,dddtheta)
    
        t_dddI=tensor_plus_cross(dddI_xx_vector(theta_dddtheta),
            dddI_xy_vector(theta_dddtheta)
        )
        return(dtheta,t_ddI,t_dddI)
    end

    
    """ An "in place" implementation of the full time derivative of theta including effects from both mass transfer and radiation reaction.""" 
    function acceleration!(du::Vector{Float64}, theta::Vector{Float64},p::Vector{Any},t::Float64; uneven::Bool=false)
        if uneven 
            #keyword argument dictates the style of parallelization to be used 
            #(mostly for benchmarking purposes)
            dtheta,t_ddI,t_dddI=ddI_dddI_from_theta_uneven(theta)
        else
            dtheta,t_ddI,t_dddI=ddI_dddI_from_theta(theta)
        end
    
        #forces due to RR
        t_dJ_rr=dJ_rr(t_ddI,t_dddI)
        @inbounds begin
            a = theta[1]
            J = theta[5]
            t_da_rr= 2.0 * a * t_dJ_rr / J
        
            #adding RR forces to MT forces
            dtheta[1]+=t_da_rr
            dtheta[5]+=t_da_rr
        end
    
        #mutate du
        du .= dtheta
    end

    """ Normal parallelization implementation of in-place acceleration! function"""    
    function acceleration_even!(du::Vector{Float64}, theta::Vector{Float64},p::Vector{Any},t::Float64)
        acceleration!(du, theta, p, t)
    end

    """ Uneven parallelization implementation of in-place acceleration! function"""    
    function acceleration_uneven!(du::Vector{Float64}, theta::Vector{Float64},p::Vector{Any},t::Float64)
        acceleration!(du, theta, p, t, uneven=true)
    end
        

    
    #"The next several blocks define functions intended to initiate the orbit beginning at mass transfer for initial masses, take an integration step forward calculating waveform, and integrate until the White Dwarf dissipates completely"
    """Calculates in km the orbital separation for the white dwarf to overflow its Roche Lobe and begin mass transfer"""
    function RL_contact(M_WD::Float64,M_NS::Float64;mu_e::Float64=2.)::Float64
        #returns separation at which a donor WD overflows its
        #Roche Lobe
        R_WD=R0_WD(M_WD,mu_e=mu_e)
        RL_a1=Roche_Limit(1.,M_WD,M_NS) #roche limit at separation of 1 km
        return(R_WD / RL_a1)#initial separation in km
    end
    
    """Calculates the angular momentum of the system at initial Roche Lobe overflow """	
    function circular_J(a::Float64,M_WD::Float64,M_NS::Float64)::Float64
        #returns angular momentum for a circular orbit
        return(( a * G * M_WD ^ 2.0 * M_NS ^ 2.0 / (M_WD + M_NS))^(1.0/2.0))
    end
end