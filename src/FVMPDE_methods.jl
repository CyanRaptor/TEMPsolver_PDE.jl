function HPDE_∂F∂x(func_F,u,Δx,dim)
    Uₓᴸ, Uₓᴿ = GeneralFluxLimiter(u,Δx)
    n = size(u,dim)
    ∂f∂x = copy(u)

    for i=1:n
        A⁺ , A⁻ = LRJacobian(func_F,u[:,i])
        ∂f∂x[:,i] = (A⁻ * Uₓᴸ[:,i] .+ A⁺ * Uₓᴿ[:,i])
    end
    return ∂f∂x
end

function HPDE_∂u∂t(prob::HPDEProblem,u)

    ∂F∂x = HPDE_∂F∂x(prob.func_F,u,prob.Grid.Δx,2)

    # future extentions
    #∂F∂x = HPDE_∂F∂x(prob.func_F,u,prob.Grid.Δx,1)
    #∂G∂y = HPDE_∂F∂x(prob.func_G,u,prob.Grid.Δy,2)
    #∂H∂z = HPDE_∂F∂x(prob.func_H,u,prob.Grid.Δz,3)
    #S = prob.func_S(u)

    ∂G∂y = 0
    ∂H∂z = 0
    #S = func_S(u)
    S = 0

    ∂u∂t = S .- (∂F∂x .+ ∂G∂y .+ ∂H∂z)

    return ∂u∂t
end


function LRJacobian(func_F,u)
    A = ForwardDiff.jacobian(func_F, u)
    #A = (abs.(A) .> 1e-10) .* A
    #println(A)
    Λ = diagm(eigvals(A))
    Λ⁺ = 0.5 .* (Λ .+ abs.(Λ))
    Λ⁻ = 0.5 .* (Λ .- abs.(Λ))
    P = eigvecs(A)
    P⁻¹ = inv(P)
    A⁺ = P * Λ⁺ * P⁻¹
    A⁻ = P * Λ⁻ * P⁻¹
    return A⁺, A⁻
end

function GeneralFluxLimiter(u,Δx)
    ϵ = 1.0e-8
    ∂u∂x(U⁺,U⁻) = (U⁺ .- U⁻) ./ Δx
    func_Ψ(r) = 0
    #func_Ψ(r) = κ_Scheme_new(r,1)
    func_Ψ(r) = SuperBee(r)
    #func_Ψ(r) = Koren(r)
    #func_Ψ(r) = 0
    func_r(Uₓ₁,Uₓ₂) = (Uₓ₁ .+ ϵ) ./ (Uₓ₂ .+ ϵ)
    func_U(Δxᵢ,r,uₓ) = (Δxᵢ / 2) .* func_Ψ(r) .* (uₓ .+ ϵ)
    Uᵢ₋₂ = [u[:,1] u[:,1] u[:,1:end-2]]
    Uᵢ₋₁ = [u[:,1] u[:,1:end-1]]
    Uᵢ   = copy(u)
    Uᵢ₊₁ = [u[:,2:end] u[:,end]]
    Uᵢ₊₂ = [u[:,3:end] u[:,end] u[:,end]]

    Uₓi₊¾ = ∂u∂x(Uᵢ₊₂,Uᵢ₊₁)
    Uₓi₊½ = ∂u∂x(Uᵢ₊₁,Uᵢ)
    Uₓi₋½ = ∂u∂x(Uᵢ,Uᵢ₋₁)
    Uₓi₋¾ = ∂u∂x(Uᵢ₋₁,Uᵢ₋₂)

    Uᴸi₊½ = Uᵢ₊₁ .- func_U(Δx,func_r(Uₓi₊½,Uₓi₊¾),Uₓi₊¾)
    Uᴸi₋½ = Uᵢ   .- func_U(Δx,func_r(Uₓi₋½,Uₓi₊½),Uₓi₊½)
    Uᴿi₊½ = Uᵢ   .+ func_U(Δx,func_r(Uₓi₊½,Uₓi₋½),Uₓi₋½)
    Uᴿi₋½ = Uᵢ₋₁ .+ func_U(Δx,func_r(Uₓi₋½,Uₓi₋¾),Uₓi₋¾)

    Uₓᴸ = (Uᴸi₊½ .- Uᴸi₋½) ./ Δx
    Uₓᴿ = (Uᴿi₊½ .- Uᴿi₋½) ./ Δx
    #return Uᴿi₊½, Uᴿi₋½, Uᴸi₊½ , Uᴸi₋½
    #return Uₓi₊¾, Uₓi₊½, Uₓi₋½, Uₓi₋¾
    return Uₓᴸ, Uₓᴿ
end





################## METHODS ####################

function Koren(r)
    return max.(0,min.(2 .* r,min.(1/3 .+ 2 /3 .* r,2)))
end

function UpWind_FO(r)
    return 0
end

function UpWind_SO(r)
    return κ_Scheme_new(r,-1)
end

function SCD(r)
    return κ_Scheme_new(r,1)
end

function UpWind_Q(r)
    return κ_Scheme_new(r,0.5)
end

function UpWind_C(r)
    return κ_Scheme_new(r,1/3.0)
end

function SuperBee(r)
    return max.(0,min.(r,1))
end
