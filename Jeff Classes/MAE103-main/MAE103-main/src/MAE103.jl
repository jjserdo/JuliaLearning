module MAE103

  using Pkg, InteractiveUtils, IJulia

  using Reexport

  using Requires
  using Conda

  @reexport using ThermofluidQuantities
  import ThermofluidQuantities: Unitful

  @reexport using Gasdynamics1D


  @reexport using Statistics

  using Roots

  repo_directory = joinpath(@__DIR__,"..")

  proj_file = Base.active_project() #Pkg.project().path
  #proj_dir = dirname(proj_file)
  #notebook_dir = joinpath(proj_dir,"notebook")

  include("quantities.jl")
  include("fluidstatics.jl")
  include("pipeflow.jl")


  #const localunits = Unitful.basefactors

  function __init__()

    #merge!(Unitful.basefactors, localunits)
    #Unitful.register(MAE103)

    @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin


      if isdefined(Main, :IJulia) && Main.IJulia.inited
        # The Pkg.build does not work if non-development package, so need to
        # ensure JIT install of matplotlib using Conda
        _hasmatplotlib() || Conda.add("matplotlib")
      else
        # For develop (e.g., CI), force re-build of PyCall with internal Python dist,
        # to make sure matplotlib is installed:
        if _iswritable(proj_file)
          ENV["PYTHON"] = ""
          Pkg.build("PyCall")
        else
           _hasmatplotlib() || error("Project file is not writable. Cannot build PyCall")
        end
      end


      # Get LaTeXStrings from PyPlot
      #using PyPlot: LaTeXStrings
      using PyCall
      @reexport using LaTeXStrings

      Plots.pyplot()
      rcParams = Plots.PyPlot.PyDict(Plots.PyPlot.matplotlib."rcParams")

      # Ensure that LaTeX stuff is handled
      rcParams["mathtext.fontset"] = "cm"

      Plots.default(markerstrokealpha = 0, legend = false,
        dpi = 100, size = (400, 300), grid = false, widen=false)

      include("plot_recipes.jl")
      include("arrows.jl")


    end

#=
    @require ViscousFlow="103da179-b3e4-57c1-99a4-586354eb2c5a" begin
      import ViscousFlow: Edges, NavierStokes, interpolatable_field, GridData,
                        VectorGridData, ScalarGridData, limits, cellsize

      include("viscousflow/fileio.jl")
      include("viscousflow/trajectories.jl")
    end
=#

  end

  _iswritable(file) = (uperm(file) >> 1) & 1 != 0
  _hasmatplotlib() = haskey(Conda._installed_packages_dict(),"matplotlib")


  function open_notebooks()
    Base.eval(Main, Meta.parse("import IJulia"))
    path = joinpath(repo_directory,"notebook")
    IJulia.notebook(;dir=path)
  end







end
