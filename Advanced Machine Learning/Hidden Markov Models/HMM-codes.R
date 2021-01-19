library(HMM)
library(entropy)

# Question 1

states <- symbols <- 1:10
symbols <- c("A","B","C","D","E","F","G","H","I","J")
emitprob <- matrix(c(.2,.2,.2,0,0,0,0,0,.2,.2,
                   .2,.2,.2,.2,0,0,0,0,0,.2,
                   .2,.2,.2,.2,.2,0,0,0,0,0,
                   0,.2,.2,.2,.2,.2,0,0,0,0,
                   0,0,.2,.2,.2,.2,.2,0,0,0,
                   0,0,0,.2,.2,.2,.2,.2,0,0,
                   0,0,0,0,.2,.2,.2,.2,.2,0,
                   0,0,0,0,0,.2,.2,.2,.2,.2,
                   .2,0,0,0,0,0,.2,.2,.2,.2,
                   .2,.2,0,0,0,0,0,.2,.2,.2), ncol = 10)


transprob <- matrix(c(.5,0,0,0,0,0,0,0,0,.5,
                      .5,.5,0,0,0,0,0,0,0,0,
                     0,.5,.5,0,0,0,0,0,0,0,
                     0,0,.5,.5,0,0,0,0,0,0,
                     0,0,0,.5,.5,0,0,0,0,0,
                     0,0,0,0,.5,.5,0,0,0,0,
                     0,0,0,0,0,.5,.5,0,0,0,
                     0,0,0,0,0,0,.5,.5,0,0,
                     0,0,0,0,0,0,0,.5,.5,0,
                     0,0,0,0,0,0,0,0,.5,.5), ncol = 10)

hmm <- initHMM(States = states, Symbols = symbols, emissionProbs = emitprob, transProbs = transprob)


# Question 2

simulated <- simHMM(hmm, 100)

obs <- simulated$observation

# Question 3

# the forward step probabilities (alpha)
frwrd <- forward(hmm, obs)

# the backward step probabilities (beta)
bkwrd <- backward(hmm, obs)

# Probabilities are converted back from log-scale
backward <- exp(bkwrd)
forward <- exp(frwrd)

# obtaining filtering distribution for each t 
filtering <- matrix(0, 10, 100)

for (i in 1:100) {
  filtering[,i] <- forward[,i]/sum(forward[,i])  
}

# obtaining smoothing distribution for each t

smoothing <- matrix(0, 10, 100)

for(i in 1:100) {
  smoothing[,i] <- (forward[,i] * backward[,i])/sum(forward[,i] * backward[,i])
}

# Using the viterbi algorithm we find the most probable path

mostProb <- viterbi(hmm, obs)


# Question 4

accuracy <- function(method_prob,sim) {
  conf <- table(sim, method_prob)
  misclass <- 1 - sum(diag(conf))/sum(conf)
  return(list(ConfusionMatrix=conf, MisclassificationRate=misclass))
}


most_probable <- function(mat,sim) {
  method_prob <- apply(mat, 2, FUN = which.max)
  accuracy(method_prob, sim)
}


most_probable(smoothing,simulated$states)
most_probable(filtering,simulated$states)

# Misclassification of the most probable path
1 - sum(mostProb == simulated$states)/length(simulated$states)

################# Question 5

nSim <- 200

simulated2 <- simHMM(hmm, nSim)
obs2 <- simulated2$observation

frwrd2 <- forward(hmm, obs2)
bkwrd2 <- backward(hmm, obs2)

backward2 <- exp(bkwrd2)
forward2 <- exp(frwrd2)

filtering2 <- matrix(0, 10, nSim)

for (i in 1:nSim) {
  filtering2[,i] <- forward2[,i]/sum(forward2[,i])  
}

smoothing2 <- matrix(0, 10, nSim)

for(i in 1:nSim) {
  smoothing2[,i] <- (forward2[,i] * backward2[,i])/sum(forward2[,i] * backward2[,i])
}

# Using the viterbi algorithm we find the most probable path

mostProb2 <- viterbi(hmm, obs2)

#### 

most_probable(smoothing2,simulated2$states)
most_probable(filtering2,simulated2$states)

1 - sum(mostProb2 == simulated2$states)/length(simulated2$states)

# The smoothing distribution is generally more accurate since it uses all the available data to make or update
# the predictions on the hidden states. This is while for making prediction on the hidden state at time t using
# filtering, we use only observations up to time t so the model learns less as it takes a smaller amount of data
# to make prediction. 

# The most probable path, obtained through the Viterbi algorithm, tries to find the hidden states that were 
# most likely to occur. For this, the algorithm maximizes the likelihood of the robot being in a state at time
# t (z_t) given the next state (z_t+1). In smoothing however, the states are updated after observing all the 
# data (into the future), decreasing the error and thus resulting in a higher accuracy.

##### Question 6


# Additionally, we evaluate this with *entropy.empirical()* function. 

entropy1 <- c()

for (i in 1:ncol(filtering)) {
  entropy1[i] <- entropy.empirical(filtering[,i])
}

entropy2 <- c()
for (i in 1:ncol(filtering2)) {
  entropy2[i] <- entropy.empirical(filtering2[,i])  
}


plot(1:200, entropy2, type="l", col="red", lwd=2, 
     main = "Comparison of entropy of filtering distributions for both models",
     xlab = "time step", ylab = "entropy")
lines(1:100, entropy1, type="l", col="blue", lwd=2)


cat("Average entropy for model 1 filtering distributions: ", mean(entropy1),"\n")
cat("Average entropy for model 2 filtering distributions: ", mean(entropy2),"\n")


# Question 7

coefs2 <- filtering[,100]

t(transprob) %*% coefs2

