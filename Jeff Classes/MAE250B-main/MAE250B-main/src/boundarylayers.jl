### Boundary layer routines

export u, v, vorticity, dudy, d2udy2, eta, η, d99, dstar, theta, Cf, falknerskan

struct BLSolution{ETAT,UT}
  β :: Float64
  Vw :: Float64
  ηmax :: Float64
  η :: ETAT
  f :: UT
  fp :: UT
  fpp :: UT
  h0 :: Float64
  d99coeff :: Float64
  dstarcoeff :: Float64
  θcoeff :: Float64
  cfcoeff :: Float64
end

"""
    u(sol::BLSolution) -> Vector{Float64}

Return the streamwise velocity array, scaled by the external velocity, ``u/U_e``
"""
u(sol::MAE250B.BLSolution) = sol.fp;

"""
    η(sol::BLSolution) -> Vector{Float64}
    eta(sol::BLSolution) -> Vector{Float64}

Return the array of wall-normal coordinates, scaled by the
boundary layer thickness factor ``\\eta = y/\\delta(x) = y/(\\nu x/U_e)^{1/2}``
"""
η(sol::MAE250B.BLSolution) = sol.η
eta(sol::MAE250B.BLSolution) = η(sol)

"""
    v(sol::BLSolution) -> Vector{Float64}

Return the array of wall-normal velocities, ``v``, scaled by the
boundary layer factor: ``v Re_x^{1/2}/U_e``
"""
v(sol::MAE250B.BLSolution) = 0.5*(sol.η.*sol.fp - sol.f)

"""
    vorticity(sol::BLSolution) -> Vector{Float64}

Return the array of vorticity, ``\\omega``, scaled by the
boundary layer factor: ``\\omega \\delta(x)/U_e``
"""
vorticity(sol::MAE250B.BLSolution) = -sol.fpp

"""
    dudy(sol::BLSolution) -> Vector{Float64}

Return the array of ``du/dy``, scaled by the
boundary layer factor: ``du/dy (\\delta(x)/U_e)``
"""
dudy(sol::MAE250B.BLSolution) = sol.fpp

"""
    d2udy2(sol::BLSolution) -> Vector{Float64}

Return the array of ``d^2u/dy^2``, scaled by the
boundary layer factor: ``d^2u/dy^2 (\\delta^2(x)/U_e)``
"""
function d2udy2(sol::MAE250B.BLSolution)
  uyy = zero(sol.fpp)
  uyy[1] = (sol.fpp[2] - sol.fpp[1])/(sol.η[2] - sol.η[1])
  for i in 2:lastindex(sol.fpp)-1
    uyy[i] = (sol.fpp[i+1] - sol.fpp[i-1])/(sol.η[i+1] - sol.η[i-1])
  end
  uyy[end] = (sol.fpp[end] - sol.fpp[end-1])/(sol.η[end] - sol.η[end-1])
  return uyy
end


"""
    d99(sol::BLSolution) -> Vector{Float64}

Return the boundary layer 99 percent thickness divided by
the boundary layer thickness factor ``\\delta_{99}(x)/\\delta(x)``
"""
d99(sol::MAE250B.BLSolution) = sol.d99coeff

"""
    dstar(sol::BLSolution) -> Vector{Float64}

Return the boundary layer displacement thickness divided by
the boundary layer thickness factor ``\\delta^{*}(x)/\\delta(x)``
"""
dstar(sol::MAE250B.BLSolution) = sol.dstarcoeff

"""
    theta(sol::BLSolution) -> Vector{Float64}

Return the boundary layer displacement thickness divided by
the boundary layer thickness factor ``\\theta(x)/\\delta(x)``
"""
theta(sol::MAE250B.BLSolution) = sol.θcoeff

"""
    Cf(sol::BLSolution) -> Vector{Float64}

Return the boundary layer skin friction coefficient scaled by
the boundary layer factor ``C_f(x)Re_x^{1/2}``
"""
Cf(sol::MAE250B.BLSolution) = sol.cfcoeff

### Solver routines ###

function fsrhs!(du,u,p,t)
  _, m, _ = p

  du[1] = u[2]
  du[2] = u[3]
  du[3] = -0.5*(m+1)*u[1]*u[3] - m*(1-u[2]^2)
end

function fs_integrate(h0,p)
  Vw, m, ηmax = p

  f0 = -2Vw
  g0 = 0.0
  gL = 1.0

  u0 = [f0;g0;h0]
  ηspan = (0.0,ηmax)


  prob = ODEProblem(fsrhs!,u0,ηspan,p)
  sol = OrdinaryDiffEq.solve(prob,Tsit5(),reltol=1e-12,abstol=1e-15)

  resid = sol[2,end] - gL

  return resid, sol

end

"""
    falknerskan(β[,Vw=0][,ηmax=10][,h0init=1.3])

Compute the Falkner-Skan boundary layer solution for parameter `β`,
where `β` can be a value larger than -0.1999 (the separation case).
It returns the data in a solution structure which can be used
to return individual results, e.g.,

`sol = falknerskan(0)`

for ``\\beta=0`` and then `u(sol)` to return the streamwise velocity profile.

You can use the optional `Vw=` argument to set a wall-normal velocity
at the wall. (This value represents ``V_w Re_x^{1/2}/U_e``.)
"""
function falknerskan(β;Vw = 0.0,ηmax = 10.0,h0init=1.3,kwargs...)

  m = β/(2-β)
  p = [Vw,m,ηmax]

  h0 = find_zero(x -> fs_integrate(x,p)[1],h0init;kwargs...)
  resid,sol = fs_integrate(h0,p)

  sol = BLSolution(convert(Float64,β),convert(Float64,Vw),convert(Float64,ηmax),_fs_eta(sol),_fs_streamfunction(sol),_fs_velocity(sol),
                    _fs_vorticity(sol),
                    convert(Float64,h0),blthickness_99(sol), blthickness_displacement(sol),
                    blthickness_momentum(sol), blskinfriction(sol))

  return sol

end

_fs_streamfunction(sol::OrdinaryDiffEq.ODESolution) = sol[1,:]
_fs_velocity(sol::OrdinaryDiffEq.ODESolution) = sol[2,:]
_fs_vorticity(sol::OrdinaryDiffEq.ODESolution) = sol[3,:]
_fs_eta(sol::OrdinaryDiffEq.ODESolution) = sol.t

function blthickness_99(sol::OrdinaryDiffEq.ODESolution)

    u, η = _fs_velocity(sol), _fs_eta(sol)

    sign_diff = sign.(u .- 0.99)  # +1 where u > 0.99, -1 otherwise
    i0 = findfirst(diff(sign_diff) .== 2.0)
    return sol.t[i0]+ (0.99-u[i0])/(u[i0+1]-u[i0])*(η[i0+1]-η[i0])
end

function blthickness_displacement(sol::OrdinaryDiffEq.ODESolution)
    u, η = _fs_velocity(sol), _fs_eta(sol)

    uhalf = 0.5*(u[2:end].+u[1:end-1])
    return sum((1 .- uhalf).*diff(η))
end

function blthickness_momentum(sol::OrdinaryDiffEq.ODESolution)
    u, η = _fs_velocity(sol), _fs_eta(sol)

    uhalf = 0.5*(u[2:end].+u[1:end-1])
    return sum(uhalf.*(1 .- uhalf).*diff(η))
end

blskinfriction(sol::OrdinaryDiffEq.ODESolution) = 2*sol[3,1]
