library(sampfling)

## samples of size 1
iter <- 1000
n <- 3
probs<-c(1,2,2)

freq <- numeric(n)
system("date")
for (i in 1:iter)
{
  j <- sampfle(n,1,probs)
  freq[j] <- freq[j]+1
}
freq

freq <- numeric(n)
for (i in 1:iter) {
  j <- samprop(3,1,probs)
  freq[j] <- freq[j]+1
}
freq
