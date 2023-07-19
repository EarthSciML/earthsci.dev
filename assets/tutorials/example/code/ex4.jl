# This file was generated, do not modify it. # hide
function ModelingToolkit.check_units(eqs...) # Skip unit enforcement for now
    nothing
    #ModelingToolkit.validate(eqs...) || @info "Some equations had invalid units. See warnings for details."
end

# Skip returning the observed variables (i.e. variables that are not state variables)
# because it is currently extremely slow to do so for this demo. 
SciMLBase.observed(sol::SciMLBase.AbstractTimeseriesSolution, sym, i::Colon) = zeros(Float64, length(sol.t))