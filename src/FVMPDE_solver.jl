
(prob::DiffEqBase.AbstractODEProblem{uType,tType,isinplace},
                          Alg::FwdEulerAlg;
                          dt=(prob.tspan[2]-prob.tspan[1])/100,
                          tstops=tType[],
                          kwargs... # ignored kwargs
                          ) where {uType,tType,isinplace}


function solve(prob::FVMPDEProblem,
               U0::Matrix{Float64})

    ODE_Function(u,p,t) = HPDE_∂u∂t(prob,u)
    #ODE_Function(u,p,t) = κ_Scheme(prob,u,1,1)
    #ODE_Function(u,p,t) = FTCS(prob,u)
    ODE_Problem = ODEProblem(ODE_Function,U0,(0,prob.Grid.S[end]));
    temp_u = DifferentialEquations.solve(ODE_Problem,dt=prob.Grid.Δs,adaptive=false)
    return temp_u
end



include("HPDE_options.jl")
include("HPDE_problem.jl")

#initialize


#function solve(prob::HPDEProblem;kwargs...)
    #function ∂u∂t(u,p,t)
#q = prob.func_Q(V0)
function HPDEsolve(prob::HPDEProblem,U0::Matrix{Float64})
    U = copy(U0)

    ODE_Function(u,p,t) = HPDE_∂u∂t(prob,u,0,-1)
    #ODE_Function(u,p,t) = FTCS(prob,u)

    for i=2:length(prob.Grid.S)
        println(i)
        prob.Grid.Δy =0
        ODE_Problem = ODEProblem(ODE_Function,U,(0,prob.Grid.Δs));
        condition(u,t,integrator) = true
        affect!(integrator) = println(integrator.t)
        cb = DiscreteCallback(condition,affect!)
        temp_u = DifferentialEquations.solve(ODE_Problem)#,dt=0.001)#,callback=cb,dt=0.1,adaptive=false)
        #println(size(temp_u.t))
        U = copy(temp_u[:,:,end])
        #U = max.(U,zeros(size(U)))
        prob.Grid.U[:,:,i] = copy(U)
    end
    return prob
end




sol = HPDEsolve(prob, U0)


plot(sol.Grid.X,sol.Grid.U[:,:,end]')

    #end
using PlotlyJS
PP = copy(sol.Grid.U)
PP[2,:,:] = sol.Grid.U[2,:,:] ./ sol.Grid.U[1,:,:]
PP[3,:,:] = sol.Grid.U[3,:,:] ./ sol.Grid.U[1,:,:]

plot(PP[:,:,end]')
#end
using Plots
gr()


function animPlot(PP,i)
    animTime = 10
    fps = 10
    _plotsInterval = convert(Int64, sol.Grid.S[end] / animTime / sol.Grid.Δs / fps)
    anims = @animate for j =1:_plotsInterval:size(PP,3)
        Plots.plot(sol.Grid.X,PP[i,:,j],label="t = $((j-1)*sol.Grid.Δs)")
    end
    gif(anims, "anim_U3.gif", fps = fps)
end
animPlot(PP,1)
