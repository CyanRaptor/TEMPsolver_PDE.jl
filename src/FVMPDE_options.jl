abstract type Schemes end
struct UpWind <: Schemes end
struct Second <: Schemes end
struct Third <: Schemes end
struct Fourth <: Schemes end

mutable struct FVMPDEOption

    Δt
    Δx
    Δs
    function HPDEOption(;kwargs...)
        Δt, Δx, Δs = nothing,nothing,nothing
        if haskey(kwargs,:dt)
            Δt = kwargs[:dt]
        end
        if haskey(kwargs,:dx)
            Δx = kwargs[:dx]
        end
        if haskey(kwargs,:ds)
            Δs = kwargs[:ds]
        end
        new(Δt, Δx, Δs)
    end
end



HPDEOption()
HPDEOption(dt=0.1)
HPDEOption(dx=0.1)
HPDEOption(ds=0.1)
HPDEOption(ds=0.1,dx=0.2)
HPDEOption(ds=0.1,dx=0.2,dt=0.3)
