module Save
    using CSV, Dates
    
    """Formats the filename as YYYYMMDD_M_WD=*M_WD*_M_NS=*M_NS*.csv
		if a filename is not provided"""
	function format_filename(M_WD::Float64,M_NS::Float64;filename::Union{String,Int}=-1)::String
		#
		if filename==-1 
		        today_str = Dates.format(Dates.now(), "yyyymmdd")
		        base_filename = "Outputs/$(today_str)_M_WD=$(M_WD)_M_NS=$(M_NS).csv"
		        filename = base_filename
		        suffix = 2
		        while isfile(filename)
		            filename = "Outputs/$(today_str)_M_WD=$(M_WD)_M_NS=$(M_NS)_$(suffix).csv"
		            suffix += 1
				end
		end
		return filename
	end
		
	"""Writes array of format given in Integrate.integrate() to a file."""
	function save_as_csv(datumses::AbstractMatrix{<:Real},
					  M_WD::Float64,M_NS::Float64; filename::Union{String,Int}=-1)::Nothing
		
		
		filename=format_filename(M_WD,M_NS;filename=filename)
		# Writes the file to a CSV titled filename

		isdir("Outputs") || mkdir("Outputs")
		#creates output directory if not already extant
		
			CSV.write(filename, 
						(t = datumses[:,1],
     					ddI_xx = datumses[:,2],
     					ddI_xy = datumses[:,3],
						a = datumses[:,4],
     					M_WD = datumses[:,5],
     					M_NS = datumses[:,6],
						R_WD = datumses[:,7],
     					J = datumses[:,8],
     					θ = datumses[:,9]))
		return nothing
	end
end