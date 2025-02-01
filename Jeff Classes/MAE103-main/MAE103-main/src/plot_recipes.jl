using RecipesBase
using ColorTypes
import PlotUtils: cgrad, palette, color_list
using LaTeXStrings

@userplot Staticpressure

@recipe function f(h::Staticpressure;ambient=Pressure(1u"atm"),numarrows=20)
    if length(h.args) != 2
      error("`staticpressure` should be given two arguments.  Got: $(typeof(h.args))")
     end
    depths, fluids = h.args

    _depths = typeof(depths) <: AbstractArray ? value.(depths) : [value(depths)]
    _fluids = typeof(fluids) <: AbstractArray ? fluids : [fluids]
    p = Pressure[]
    d = Depth[]
    fs = Depth(0u"m")
    push!(p,ambient)
    push!(d,fs)
    prepend!(_depths,value(fs))
    dd = ustrip(_depths[end])/numarrows
    for i in eachindex(_fluids)
        di = Depth.(range(0,ustrip(_depths[i+1]-_depths[i]),step=dd)*unit(_depths[i+1]))
        pi = Pressure.(p[end] .+ SpecificWeight(_fluids[i])*di)
        append!(d,Depth.(d[end] .+ di))
        append!(p,pi)
    end

    flip := true
    xguide := "Pressure"
    yguide := "Depth"
    widen := false
    ymirror := true
    @series begin
        p, d
    end
    @series begin
        seriestype := Plots.quiver
        quiver := (-ustrip.(p),zeros(length(p)))
        p,d
    end
end
