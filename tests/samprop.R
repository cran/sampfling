library(sampfling)

N <- 1000

# loop until valid probs
repeat
  {
    probs <- runif(4)
    if (max(probs*3) < 1) break
  }

freq <- array(0,c(4,4,4))
for (i in 1:4)
  for (j in 1:4)
  for (k in 1:4)
  if (i<j&&j<k) freq[i,j,k] <- probs[i]*probs[j]*probs[k]
# Expected freqs
freq*N/sum(freq) 

freq <- array(0,c(4,4,4))
for (i in 1:N)
{
  j <- sort(samprop(4,3,probs))
  j1 <- j[1]; j2 <- j[2]; j3 <- j[3]
  freq[j1,j2,j3] <- freq[j1,j2,j3]+1
}
# Sampled freqs
freq
