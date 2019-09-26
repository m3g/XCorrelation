
module XCorrelation

  using Statistics
  using ProgressMeter
  using JPleXL # from https://github.com/mcubeg/JPleXL

  include("./src/CorrelationData.jl")
  include("./src/index4D.jl")
  include("./src/contact_correlations.jl")
  include("./src/getcorr.jl")
  include("./src/getconstraints_Luciano.jl")
  include("./src/meancorr.jl")

end


