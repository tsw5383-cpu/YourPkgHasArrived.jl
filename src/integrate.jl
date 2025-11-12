module Integrate
    using OrdinaryDiffEq
    include("math.jl")
    include("save.jl")
    
    """
        Returns initial conditions for the beginning of Roche Lobe overflow as a vector.
        Vector output is [a, M_WD, M_NS, R_WD, J, θ], where:
        a = total orbital separation of the binary.
        M_WD = mass of the white dwarf.
        M_NS = mass of the neutron star.
        R_WD = radius of the white dwarf.
        J = total angular momentum of the binary.
        θ = angular parameter
        ```julia
        initial_circular_orbit(M_WD::Float64;M_NS::Float64=1.35,mu_e::Float64=2.)
        ```
        mu_e = mean molecular mass per electron (default 2 for white dwarfs of all but the most exotic chemical compositions (e.g., ONeMG)).
    """
    function initial_circular_orbit(M_WD::Float64;M_NS::Float64=1.35,mu_e::Float64=2., mode::String="Debug")::Vector{Float64}
        #returns keplerian orbital parameters, theta, at initial Roche Lobe overflow
        #M_NS=1.35 M_\odot is typical NS mass
        if mode=="Debug"
            @assert M_WD<1.4 #Chandrasekar mass, beyond which we cannot calculate R_WD(see initial function)
        end
        a=Math.RL_contact(M_WD,M_WD,mu_e=mu_e)
        J=Math.circular_J(a,M_WD,M_WD)
        R_WD=Math.R0_WD(M_WD)
        return([a, M_WD, M_NS, R_WD, J, 0.])
    end
    
    #"Below is the heart of the integration process"
    """Integrates theta forward one step with time interval either given explicitly (dt!=-1) or as a fraction of the orbital period determined by dt_period. Returns the absolute time step taken and the updated theta parameters."""
    function one_step!(theta::Vector{Float64}; dt_period::Float64=1e-2, dt::Float64=-1., mode::String="Debug", uneven::Bool=false)::Tuple{Float64, Vector{Float64}}
        if mode=="Debug"
            @assert length(theta)==6
        end
        #integrates forward one step with time interval either given explicitly (dt!=-1)
        #or as a fraction of the orbital period determined by dt_period
        if dt==-1
            @inbounds begin
                a = theta[1]
                M_WD = theta[2]
                M_NS = theta[3]
                dt=dt_period*Math.period(a, M_WD + M_NS)
            end
        end
        if uneven
            prob = ODEProblem(Math.acceleration_uneven!,theta,(0.,dt),[])
        else
            prob = ODEProblem(Math.acceleration_even!,theta,(0.,dt),[])
        end
        theta .= solve(prob, alg=Vern9()).u[end]
        dt,theta
    end
    
    """Given an initial set of masses, `integrate` solves for the equations of motion from initial Roche Lobe overflow until the eventual complete dissipation of the white dwarf or the integration finishes according to user defined cutoff (either a limit on the number of orbits completed or on the integration steps taken). `integrate` stores and returns the second derivative of the quadrupole moment as a function of time for the orbit. The function also saves the run to a file named by the date and initial conditions if `write_file`. """
    function integrate(M_WD::Float64; M_NS::Float64=1.35,mu_e::Float64=2., 
            orbit_limit::Float64=1e5, step_limit::Int=10^7, mass_limit::Float64=1e-3,
                   write_file::Bool=true, filename::Union{String,Int}=-1,uneven::Bool=false, dt_period=1e-3)::AbstractMatrix{<:Real}
        #returns time series second derivatives of the quadrupole moment
        t=0.
        theta=initial_circular_orbit(M_WD,M_NS=M_NS,mu_e=mu_e)
        ddI=Math.ddI_from_theta(theta)
        datumses=zeros(step_limit+1,9)
        datumses[1,:].=[t, ddI[1,1], ddI[1,2], theta... ]
        i=1

        @inbounds begin
            phase = theta[6]
            M_WD = theta[2]
            while phase<2*π*orbit_limit && i<step_limit && theta[2]>mass_limit
                #integrate until time or step limit are elapsed or WD completely dissipates
                
                dt,theta=one_step!(theta; dt_period=dt_period,uneven )
                ddI=Math.ddI_from_theta(theta)
                i+=1
                t+=dt
            
                datumses[i,:].=[t, ddI[1,1], ddI[1,2], theta... ]
            end
        end
    
        datumses=datumses[1:i,:]
    
        if write_file #saves datumses to a csv in directory "Output"
            Save.save_as_csv(datumses,M_WD,M_NS,filename=filename)
        end
        
        return(datumses)
    end
end