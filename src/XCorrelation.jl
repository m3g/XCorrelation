
module XCorrelation

  using Statistics
  using ProgressMeter
  using PDBTools # from https://github.com/m3g/PDBTools

  include("./CorrelationData.jl")
  include("./index4D.jl")
  include("./contact_data.jl")
  include("./contact_correlations.jl")
  include("./getcorr.jl")
  include("./getconstraints_Luciano.jl")
  include("./meancorr.jl")
  include("./cbin.jl")

  include("./censoni_information.jl")

  include("./biserial_correlation.jl")

end


