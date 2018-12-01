# Parameters
# write.csv(parameters, "../parameters.csv", row.names = FALSE)

# Standard

## Experiment 1
N      <- c(1000, 10000)
D      <- c(5, 10, 30, 200, 600)
P      <- c(.05, .1, .2, .3, .4, .5)
Lambda <- c(0.8, 0.9, 1, 1.1, 1.2)

Sign <- c("pos", "neg")

sampSize <- 100
graphNum <- length(N)*length(D)*length(P)*length(Lambda)*length(Sign)*sampSize


### Parameter matrix
parameters  <- NULL
graph_stats <- NULL

for (n in N) {
  for (d in D) {
    for (p in P) {
      for (lambda in Lambda){
        for (sign in Sign){

          # Parameters
          paramRow <- paramGenerator(n, d, p, lambda, sign)
          graphRow <- graph_stats(paramRow)

          # Append
          parameters  <- rbind(parameters, paramRow)
          graph_stats <- rbind(graph_stats, graphRow)

        }
      }
    }
  }
}


## sampSize = 100
parameters_100  <- parameters[rep(seq_len(nrow(parameters)), each = sampSize),]
graph_stats_100 <- graph_stats[rep(seq_len(nrow(graph_stats)), each = sampSize),]


### Graph indices
graphNums <- seq(1, graphNum)

row.names(parameters_100) <- graphNums
parameters_100 <- cbind(graphNums, parameters_100)

### Complete graph_stats
graph_stats     <- cbind(parameters, graph_stats)
graph_stats_100 <- cbind(parameters_100, graph_stats_100)


## memReq

memReq <- read.csv("data-raw/memory_requirements.csv")
memReq <- memReq$MEM



#############################################

# Overlap

## Experiment 1
N_overlap <- c(1000, 10000)
S         <- c(0, 0.01, 0.1, 0.33, 0.5, 1)
p_overlap <- 0.5
D_1k      <- c(5, 10, 30, 200)
D_10k     <- c(5, 10, 30, 200, 600)
cLength   <- 10

sampSize_overlap <- 100
graphNum_overlap <- (length(D_1k) + length(D_10k))*length(S)*cLength*sampSize_overlap


### Parameter matrix
parameters_overlap  <- NULL

parameters_overlap_base <- rbind(expand.grid(n = N_overlap[1], d = D_1k),
                                 expand.grid(n = N_overlap[2], d = D_10k))


for (i in 1:nrow(parameters_overlap_base)) {

  baseRow <- parameters_overlap_base[i,]
  n <- baseRow$n
  d <- baseRow$d

  C <- seq(0, 2*d, length.out = cLength)

  for (c in C) {

    paramRow <- data.frame(n, p = p_overlap, d, s = S, c)
    parameters_overlap  <- rbind(parameters_overlap, paramRow)

  }

}

parameters_overlap <- dplyr::arrange(parameters_overlap, n, d, s, c)  # to match Vaishaki's parameter list

### Graph stats matrix
graph_stats_overlap <- NULL

for (i in 1:nrow(parameters_overlap)) {

  paramRow <- parameters_overlap[i,]
  graphRow <- graph_stats_overlap(paramRow)

  # Append
  graph_stats_overlap <- rbind(graph_stats_overlap, graphRow)

}


## sampSize = 100
parameters_overlap_100  <- parameters_overlap[rep(seq_len(nrow(parameters_overlap)), each = sampSize_overlap),]
graph_stats_overlap_100 <- graph_stats_overlap[rep(seq_len(nrow(graph_stats_overlap)), each = sampSize_overlap),]


### Graph indices
graphNums_overlap <- seq(1, graphNum_overlap)

row.names(parameters_overlap_100) <- graphNums_overlap
parameters_overlap_100 <- cbind(graphNum = graphNums_overlap, parameters_overlap_100)


### Complete parameters
param_vai <- read.csv("data-raw/parameters_vai.csv")
ks <- param_vai$ks

parameters_overlap_100 <- cbind(parameters_overlap_100, ks)


### Complete graph_stats
graph_stats_overlap     <- cbind(parameters_overlap, graph_stats_overlap)
graph_stats_overlap_100 <- cbind(parameters_overlap_100, graph_stats_overlap_100)


# memReq
S_sans_1 <- S[!S %in% 1.00]

memReq_overlap <- rbind(expand.grid(n = N_overlap[1], d = D_1k,  s = S_sans_1),
                        expand.grid(n = N_overlap[2], d = D_10k, s = S_sans_1))
memReq_overlap <- dplyr::arrange(memReq_overlap, n, d, s)


mem_overlap_1k  <- c(1, 1, 1, 1, 1,  # 1k, 5
                     1, 1, 1, 1, 1,  # 1k, 10
                     2, 2, 2, 3, 3,  # 1k, 30
                     3, 5, 5, 6, 6)  # 1k, 200

mem_overlap_10k <- c( 2,  2,  3,  3,  3,  # 10k, 5
                      3,  3,  3,  3,  5,  # 10k, 10
                      6,  6,  6,  6,  6,  # 10k, 30
                     12, 16, 16, 26, 32,  # 10k, 200
                     30, 32, 38, 60, 70)  # 10k, 600


memReq_overlap <- cbind(memReq_overlap, mem = c(mem_overlap_1k, mem_overlap_10k))


#############################################

# 3 groups

## Experiment 2


N <- 1000
D <- c(5, 10, 30, 200)
P <- list(c(0.334, 0.333, 0.333),
          c(0.8, 0.1, 0.1),
          c(0.45, 0.45, 0.1),
          c(0.6, 0.3, 0.1))

#RB | create grid
#RB | add lambdas
#RB | call graph generator for 3 group
#RB | in matrix generator, include the set of calculations; pick the largest root.

