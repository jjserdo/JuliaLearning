using OrdinaryDiffEq
using Interpolations
using LinearAlgebra

export orrsomm, grad_optimization

"""
    orrsomm(Rstar,αstar,cinit) -> ComplexF64

Given the Reynolds number based on displacement thickness `Rstar`,
and the wavenumber scaled by displacement thickness `αstar`, find
the corresponding complex wavespeed `c`. Use the initial guess `cinit`
to start reasonably close.
"""
function orrsomm(Rstar,αstar,cinit)

  # Numerical parameters
  toler = 1e-6  # Tolerance for shooting, integration, etc
  ymaxos = 10  # The upper limit of vertical range for eigenvalue calculation

  # Solve the Falkner-Skan problem for Blasius case
  fssol = falknerskan(0.0,ηmax=13.0)

  # Get the velocity U(y) and second derivative U''(y) as interpolatable_field
  # functions
  y = eta(fssol)
  u_interp = LinearInterpolation(y, u(fssol))
  udd_interp = LinearInterpolation(y, d2udy2(fssol))

  # Re-scale the given Re and α from δ* to δ = sqrt(νx/Uinf)
  Re_delta = Rstar/dstar(fssol)
  αdelta = αstar/dstar(fssol)

  pinit = sqrt(αdelta^2 + im*αdelta*Re_delta*(1-cinit))
  if (real(pinit) < 0)
    pinit = -pinit
  end

  par = Dict("ymaxos" => ymaxos,
            "tol" => toler,
            "delv" => [0.001;0.001],
            "alpha" => αdelta,
            "Re" => Re_delta,
            "c" => cinit, "p" => pinit,
            "u_interp" => u_interp,
            "udd_interp" => udd_interp)

  vinit = [real(cinit);imag(cinit)]

  sol = grad_optimization(computeFos,vinit, par)

  return c = sol.v[1]+im*sol.v[2]

end


"""
    computeFos(v,par)

Given the value ``v``, a vector containing the real and imaginary parts
of the complex wave speed, calculate the value ``F(v)``
for the Orr-Sommerfeld equations.
"""
function computeFos(v,par)

    α = par["alpha"]
    Re = par["Re"]
    ymaxos = par["ymaxos"]

    cr, ci = v
    c = cr + im*ci
    par["c"] = c

    p = sqrt(α^2 + im*α*Re*(1-c))
    if (real(p) < 0)
        p = -p
    end

    Finit = zeros(ComplexF64,6,1);
    Finit[1] = 1
    Finit[2] = -(α+p)
    Finit[3] = α^2 + α*p + p^2
    Finit[4] = α*p
    Finit[5] = -α*p*(α+p)
    Finit[6] = α^2*p^2

    yspan = (ymaxos,0.0)

    prob = ODEProblem(orrsommrhs!,Finit,yspan,par)
    sol = OrdinaryDiffEq.solve(prob,Tsit5(),reltol=1e-10,abstol=1e-10)

    ulast = sol.u[end]
    Dr = real(ulast[1])/maximum(abs.(ulast))
    Di = imag(ulast[1])/maximum(abs.(ulast))
    Fvec = [Dr;Di]

    return Fvec

end

"""
    orrsommrhs!(df,f,par,y)

Return the right hand side `df` (in place) for the Orr-Sommerfeld
equations in compound matrix form. The parameters `par` are
a Dict of parameters.
"""
function orrsommrhs!(df,f,par,y)
    α = par["alpha"]
    Re = par["Re"]
    c = par["c"]
    p = par["p"]
    u_interp = par["u_interp"]
    udd_interp = par["udd_interp"]

    uy = u_interp(y)
    uddy = udd_interp(y)

    a1 = 0.0
    a2 = 2*α^2 + im*α*Re*(uy - c)
    a3 = 0.0
    a4 = -α^4 - im*α*Re*(α^2*(uy-c)+uddy)

    #fact = kx + p;
    fact = 0

    df[1] = fact*f[1] + f[2]
    df[2] = fact*f[2] + f[3] + f[4]
    df[3] = a3*f[1] + a2*f[2] + (a1+fact)*f[3] + f[5]
    df[4] = fact*f[4] + f[5]
    df[5] = -a4*f[1]+a2*f[4] + (a1+fact)*f[5] + f[6]
    df[6] = -a4*f[2] - a3*f[4] + (a1+fact)*f[6]

    return df

end


#### Support routines

"""
    jac(Ffcn,v0,F0,par; delv0 = 0.001)

Calculate the Jacobian of the function ``F(v)`` at
``v_0`` using finite difference methods. If
`delv` is included as a key in the `par` Dict,
then use this to set the differencing intervals;
otherwise use the default `delv0` in each direction.
"""
function jac(computeF,v,F,par; delv0 = 0.001)

    delv = similar(v)
    fill!(delv,1.0)

    delv = get(par,"delv",delv0*delv)

    n = length(v)
    vnew = zero(v)
    J = zeros(eltype(F),n,n)
    for k in 1:length(F)

      vnew .= v
      vnew[k] += delv[k]
      Fnew = computeF(vnew,par)

      J[:,k] = (Fnew - F)/delv[k]
    end
    return J

end

struct OptimizationSolution{VT,ET}
    v :: VT
    err :: ET
    check :: Integer
end

"""
    grad_optimization(Ffcn,vi,par) -> OptimizationSolution

Return the solution of a gradient descent, in which the
function attempts to find the root ``v`` that satisfies
``F(v) = 0``, where ``F`` is the error vector
and ``v`` is the left-side boundary guess vector.
It uses a line search along the gradient direction to improve the
method convergence. The root ``v`` is given as the `v` field
of the output structure.
"""
function grad_optimization(computeF,v,par)
    stpmx = 100
    tolf = 5e-6
    tolmin = 1e-6
    tolv = 1e-7
    maxits = 200

    n = length(v)

    # compute initial Fvec
    Fvec = computeF(v,par)
    f = 0.5*real(transpose(conj(Fvec))*Fvec)

    # Check if initial guess is sufficient
    err = norm(Fvec,Inf)
    earray = [];
    if (err < 0.01*tolf)
        push!(earray,err)
        check = 0
        vf = v
        return
    end

    # Compute maximum step size
    stpmax = stpmx*max(norm(v),n)

    vold = zero(v)
    vf = zero(v)

    # Iterate
    for it in 1:maxits
        J = jac(computeF,v,Fvec,par)

        # Compute grad(f)
        g = Fvec'*J

        # Save old values
        vold .= v
        fold = f

        # Compute full Newton step
        p = -J\Fvec

        # Call line search algorithm.  This calculates new V and f, as well
        # as new Fvec.
        out = lnsrch(computeF,vold,fold,g,p,stpmax,par)
        check = out.check
        f = out.f
        Fvec .= out.Fvec
        v .= out.v

        # Test for convergence on function values
        err = norm(Fvec,Inf)
        push!(earray,err)

        if (err < tolf)
            check = 0
            vf .= v
            return OptimizationSolution(vf,earray,check)
        end

        # Check for spurious convergence at grad(f) = 0
        if (check == 1)
            den = max(f,0.5*n)
            test = 0
            for k in 1:n
                test = max(test,abs(g[k])*max(abs(v[k]),1)/den)
            end
            if (test < tolmin)
              check = 1
            else
              check = 0
            end
            vf .= v
            return OptimizationSolution(vf,earray,check)
        end

        # Test for convergence on dV
        test = 0
        for k in 1:n
            test = max(test,abs(v[k]-vold[k])/max(abs(v[k]),1))
        end

    end

    return OptimizationSolution(vf,earray,-1)

end


struct LineSearchSolution{VT,FT}
    v :: VT
    f :: Float64
    Fvec :: FT
    check :: Integer
end


function lnsrch(computeF,vold,fold,g,p,stpmax,par)

    alf = 1e-4
    tolv = 1e-7

    check = 0
    n = length(vold)

    # Rescale step size if too big
    magp = norm(p)
    if (magp > stpmax)
      p *= stpmax/magp
    end

    slope = g*p

    if (slope >= 0)
      error("roundoff problem")
    end

    # Compute min lambda
    test = 0
    for k in 1:n
      test = max(test,abs(p[k])/max(abs(vold[k]),1))
    end
    alamin = tolv/test

    # Start with full step
    alam = 1

    v = zero(vold)

    rhs = zeros(2)

    while true # Continue until 'return' reached

      v .= vold + alam*p

      Fvec = computeF(v,par)
      f = 0.5*real(transpose(conj(Fvec))*Fvec)

      if (f < fold + alf*alam*slope)
        return LineSearchSolution(v,f,Fvec,check)
      else
        if (alam == 1)
          # For first time, use quadratic approximation
          tmplam = -slope/2/(f-fold-slope)
        else
          # Use cubic approximation
          rhs[1] = f - fold - alam*slope
          rhs[2] = f2 - fold - alam2*slope

          a = (rhs[1]/alam^2-rhs[2]/alam2^2)/(alam-alam2)
          b = (-alam2*rhs[1]/alam^2+alam*rhs[2]/alam2^2)/(alam-alam2)
          if (a == 0)
            tmplam = -slope/2/b
          else
            disc = b^2-3*a*slope
            if (disc < 0)
              tmplam = 0.5*alam
            elseif (b <= 0)
              tmplam = (-b+sqrt(disc))/(3*a)
            else
              tmplam = -slope/(b+sqrt(disc))
            end
          end
          if (tmplam > 0.5*alam)
            tmplam = 0.5*alam
          end
        end
      end
      alam2 = alam
      f2 = f
      alam = max(tmplam,0.1*alam)

    end

    if (alam < alamin)
      v .= vold
      check = 1
    end

    return LineSearchSolution(v,f,Fvec,check)
end
