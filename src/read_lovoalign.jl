#
# Read lovoalign file
#
function read_lovoalign(directory)
  lovoalign_file_name = directory*"/lovoalign.log"
  scores = Vector{Float64}(undef,1000)
  lovoalign_file = open(lovoalign_file_name,"r")
  n50 = 0
  maxscore = 0.
  for line in eachline(lovoalign_file)
    if line[1:1] == "#"
        continue
    end
    data = split(line)
    model_name = data[1]
    imodel = parse(Int64,model_name[3:10])+1
    score = parse(Float64,data[3])
    if score > 0.5
      n50 = n50 + 1
    end
    maxscore = max(score,maxscore)
    scores[imodel] = score
  end
  close(lovoalign_file)
  return scores, n50, maxscore
end

