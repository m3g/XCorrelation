
module XCorrelation

  using Statistics
  using ProgressMeter
  using JPleXL # from https://github.com/mcubeg/JPleXL

  include("./CorrelationData.jl")
  include("./index4D.jl")
  include("./contact_correlations.jl")
  include("./getcorr.jl")
  include("./getconstraints_Luciano.jl")
  include("./meancorr.jl")

  include("./censoni_information.jl")

end


