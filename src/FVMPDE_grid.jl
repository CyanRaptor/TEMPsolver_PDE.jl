abstract type AbstractGridComponent end

mutable struct ExtendableGrid{Tc,Ti}
    components::Dict{Type{<:AbstractGridComponent},Any}
    ExtendableGrid{Tc,Ti}() where{Tc,Ti} =new(Dict{Type{<:AbstractGridComponent},Any}())
end


function simplexgrid(_X::AbstractVector)
    X=collect_or_assign(_X)
    #    is_monotone(X) || error("X not monotone")
    coord=reshape(X,1,length(X))
    cellnodes=zeros(Cint,2,length(X)-1)
    cellregions=zeros(Cint,length(X)-1)
    for i=1:length(X)-1
        cellnodes[1,i]=i
        cellnodes[2,i]=i+1
        cellregions[i]=1
    end
    bfacenodes=Array{Cint}(undef,1,2)
    bfaceregions=zeros(Cint,2)
    bfacenodes[1,1]=1
    bfacenodes[1,2]=length(X)
    bfaceregions[1]=1
    bfaceregions[2]=2
    grid=simplexgrid(coord,
                     cellnodes,
                     cellregions,
                     bfacenodes,
                     bfaceregions)
    grid[XCoordinates]=X
    grid
end


function FVMPDEGrid(cellCenters::Array{Tc,2},
                    cellnodes::Array{Ti,2},
                    cellregions::Array{Ti,1},
              bfacenodes::Array{Ti,2},
              bfaceregions::Array{Ti,1}
              ) where {Tc,Ti}
    @assert size(coord,2)>0
    dim=size(coord,1)
    # for 2D it should be Array{Float64,4}
    # for 3D it should be Array{Float64,5}
    # this will be fixed later
    U::Array{Float64,3};

    Δs::Float64;
    S::Array{Float64,1};

    Δx::Float64;
    X::Array{Float64,1};
    F::Array{Float64,3};

    Δy::Float64;
    Y::Array{Float64,1};
    G::Array{Float64,3};

    Δz::Float64;
    Z::Array{Float64,1};
    H::Array{Float64,3};

    _grid=ExtendableGrid{Tc,Ti}()
    _grid[Coordinates]=coord
    _grid[CellNodes]=cellnodes
    _grid[CellRegions]=cellregions
    _grid[CellGeometries]=VectorOfConstants(eltype,length(cellregions))
    _grid[BFaceNodes]=bfacenodes
    _grid[BFaceRegions]=bfaceregions
    _grid[BFaceGeometries]=VectorOfConstants(btype,length(bfaceregions))
    _grid[CoordinateSystem]=csys
    return _grid

    # this function only works for 1D
    # for 2D and 3D, two other functions should be added later
    function HPDEGrid(tspan,xspan,U0,func_F,option::HPDEOption)
        println(size(U0))
        nvars = size(U0,1)
        nx = size(U0,2)
        Δx = option.Δx = (xspan[2]-xspan[1])/nx
        X = xspan[1] + 0.5*Δx : Δx : xspan[2] - 0.5*Δx

        if option.Δs === nothing
            option.Δs = Δx/10.0 # according to Zahir equation 1.32
        end
        Δs = option.Δs
        S = tspan[1] : Δs : tspan[2]

        U = zeros(nvars,length(X),length(S));
        F = zeros(nvars,length(X),length(S));
        U[:,:,1] = U0
        F[:,:,1] = func_F(U[:,:,1])
        new(nvars,U,Δs,S,Δx,X,F)
    end
end

HPDEGrid((0,0.2),(0,1.0),V0,func_F,HPDEOption())


X=collect(0:0.5:1.5)
Y=collect(0:0.2:0.2)
Z=collect(0:0.2:0.2)

A = simplexgrid(X)


A = simplexgrid(X,Y)

A = simplexgrid(X,Y,Z)


A[Coordinates]


A[CellNodes]


A[CellRegions]


A[XCoordinates]


A[BFaceNodes]


A[BFaceRegions]


A[BFaceGeometries]
