
# Port of original DeltaCon to R
# Baoxu(Dash) Shi
# Data Sciense Group
# University of Notre Dame
#
# Citationx
#   D. Koutra, J. T. Vogelstein, and C. Faloutsos:
#   DeltaCon: A Principled Massive-Graph Similarity Function.
#   SIAM 2013: 162-170.
#   D. Koutra, T. Ke, U. Kang, D. H. Chau, H. K. Pao, C. Faloutsos:
#   Unifying Guilt-by-Association Approaches: Theorems and Fast Algorithms.
#   ECML/PKDD (2) 2011: 245-260

library(Matrix)
library(SparseM)
library(pracma)

.MAX_POWER = 20
.p = 0.51

output_time <- function(debug, tim, s) {
  if(debug) {
    print(paste(s, "user time:", tim[1], "system time:", tim[2], "elapsed time:", tim[3]))
  }
}

inverse_lbp <- function(graph, nnodes, priors=NULL, times = 10, debug = FALSE) {

  # Sparse identity matrix
  I = NULL
  tim <- system.time(
    {
      I <- .sparseDiagonal(nnodes, x = 1, shape = "t")
    })
  output_time(debug, tim, "Create sparse identity matrix")

  # Sparse degree-diagonal matrix, D[i,i] = sum(graph[i,])
  x = NULL
  tim <- system.time(
    {
      x <- rowSums(graph, sparseResult = TRUE)
    })
  output_time(debug, tim, "Calculate degree of each node")

  D = NULL
  tim <- system.time(
    {
      D <- sparseMatrix(c(1:nnodes), c(1:nnodes), x = x, dims = c(nnodes, nnodes))
    })
  output_time(debug, tim, "Create degree diagonal matrix")

  # Compute about-half homophily factor to guarantee covergence
  c1 = sum(D) + 2
  c2 = sum(D^2) - 1
  h_h = sqrt((-c1 + sqrt(c1^2 + 4 * c2)) / (8 * c2))

  # Compute constant ah and ch
  ah = 4 * h_h^2 / (1 - 4 * h_h^2)
  ch = 2 * h_h / (1 - 4 * h_h^2)

  # Invert matrices M1 and M2
  M = NULL
  tim <-system.time({
    M = ch * graph  - ah * D
  })
  output_time(debug, tim, "Initialize Invert matrix")



  # Calculate inverse of M
  if (is.null(priors)) {
    inv_ = I
    mat_ = M
    pow = 1
    tim <- system.time({
      while(max(mat_) > 1e-09 && pow < .MAX_POWER) {
        inv_ = inv_ + mat_
        mat_ = mat_ %*% M
        pow = pow + 1
      }
    })
    output_time(debug, tim, "Invert of matrix")
    return(inv_)
  } else {
    final_mat <- NULL
    tim <- system.time({
      for(i in c(1:times)) {
        inv_ = matrix(priors[, i], nnodes, 1)
        mat_ = matrix(priors[, i], nnodes, 1)
        pow = 1
        while(max(mat_) > 1e-09 && pow < .MAX_POWER) {
          mat_ = M %*% mat_
          inv_ = inv_ + mat_
          pow = pow + 1
        }
        if (i == 1) {
          final_mat <- matrix(inv_)
        } else {
          final_mat <- cbind(final_mat, matrix(inv_))
        }
      }
    })
    output_time(debug, tim, "Approximate invert of matrix")
    return(final_mat)
  }

}

init_priors_percent <- function(percent, nnodes) {
  times <- ceiling(1 / percent)
  init_nodes <- floor(percent * nnodes)

  rand_vector <- runif(nnodes)

  rand_mat <- repmat(matrix(rand_vector, nnodes, 1), 1, times)

  for(i in c(1:times)) {
    rand_mat[rand_mat[,i] >= (i-1)*percent & rand_mat[,i] < i *percent, i] <- 1
    rand_mat[rand_mat[,i] != 1, i] <- 0
  }

  return(rand_mat)
}

delta_con <- function(graph1, graph2, nnodes,
                      method = "naive", percent = 0.1, debug = FALSE, symmetrical = TRUE) {
  if(ncol(graph1)!=2 || ncol(graph2)!=2) {
    print("Input file should be data.frame with two cols(src dst).")
    return(0);
  }
  colnames(graph1) <- c("src", "dst")
  colnames(graph2) <- c("src", "dst")

  # Construct sparse adjacent matrix from edge list
  node_vector <- NULL
  tim <- system.time({
    g1 <- sparseMatrix(graph1$src, graph1$dst, x=1, dims=c(nnodes, nnodes))
    g2 <- sparseMatrix(graph2$src, graph2$dst, x=1, dims=c(nnodes, nnodes))
  })
  output_time(debug, tim, "Construct sparse adjacent matrix")

  # Change directed graph to undirected
  if(symmetrical) {
    tim <- system.time({
      g1 <- g1 + t(g1)
      g2 <- g2 + t(g2)
    })
    output_time(debug, tim, "Change to directed graph")
  }

  if (method == "fast") {
    # Number of groups to be initialized
    ngroups <- ceiling(1 / percent)

    repetitions <- 10
    t_all <- rep(0, repetitions)
    sim <- rep(0, repetitions)
    for(i in c(1:repetitions)) {
      priors <- NULL
      tim <- system.time({
        priors <- init_priors_percent(percent, nnodes)
      })
      output_time(debug, tim, "Calculate priors")

      inv1 <- inverse_lbp(g1, nnodes, priors, ngroups, debug = debug) * (.p - 0.5)
      inv2 <- inverse_lbp(g2, nnodes, priors, ngroups, debug = debug) * (.p - 0.5)
      sim[i] <- 1 / 1 / (1 + sqrt(sum( (sqrt(inv1) - sqrt(inv2))^2 )))
    }
    delta_con <- mean(sim)
    return(delta_con)
  } else {
    # Naive FaBP
    inv1 <- inverse_lbp(g1, nnodes, debug = debug) * (.p - 0.5)
    inv2 <- inverse_lbp(g2, nnodes, debug = debug) * (.p - 0.5)

    # Compute DeltaCon similarity score
    delta_con <- 1 / (1 + sqrt(sum( (sqrt(inv1) - sqrt(inv2))^2 )))
    return(delta_con)
  }
}

#delta_con_example()
delta_con_example <- function() {
  g1 <- as.data.frame(cbind(c(1,1,2,2,3,3,4,8,5,"a"),
                            c(2,3,4,5,6,7,8,9,10,"b")))
  g2 <- as.data.frame(cbind(c(1,1,2,2,3,3,4,8,5,9,10,5,6),
                            c(2,3,4,5,6,7,8,9,10,11,11,12,12)))

  print(delta_con(g1,g1,max(g1,g2), debug = TRUE))
  print(delta_con(g1,g2,max(g1,g2), debug = TRUE))
  print(delta_con(g2,g1,max(g1,g2), debug = TRUE))
  print(delta_con(g2,g2,max(g1,g2), debug = TRUE))

  #print(delta_con(g1,g1,max(g1,g2), node_list = c(1,2,3,4,5,6,7,8,9,10), method="fast", debug = TRUE))
  #print(delta_con(g1,g2,max(g1,g2), node_list = c(1,2,3,4,5,6,7,8,9,10), method="fast", debug = TRUE))
  print(delta_con(g2,g1,max(g1,g2), method="fast", debug = TRUE))
  print(delta_con(g2,g2,max(g1,g2), method="fast", debug = TRUE))

  print(delta_con(g1,g1,max(g1,g2), method="fast", symmetrical = FALSE))
  print(delta_con(g1,g2,max(g1,g2), method="fast", symmetrical = FALSE))
  print(delta_con(g2,g1,max(g1,g2), method="fast", symmetrical = FALSE))
  print(delta_con(g2,g2,max(g1,g2), method="fast", symmetrical = FALSE))

  print(delta_con(g1,g1,max(g1,g2), symmetrical = FALSE))
  print(delta_con(g1,g2,max(g1,g2), symmetrical = FALSE))
  print(delta_con(g2,g1,max(g1,g2), symmetrical = FALSE))
  print(delta_con(g2,g2,max(g1,g2), symmetrical = FALSE))

  g3 <- as.data.frame(cbind(c(2,3),
                            c(1,1)))
  g4 <- as.data.frame(cbind(c(2,3),
                            c(4,4)))

  g5 <- as.data.frame(cbind(c(2,3,3,1,5,6,6),
                            c(1,1,2,4,4,4,5)))
  g6 <- as.data.frame(cbind(c(2,3,3,5,6,6),
                            c(1,1,2,4,4,5)))
}


delta_con_example()

deltacon_igraph= function(g1_net,g2_net){

  #g1_net= hfpef_net

  g1= igraph::as_data_frame(g1_net, what= "edges")
  g1= g1[,1:2]

  #g2_net= hfref_net

  g2= igraph::as_data_frame(g2_net, what= "edges")
  g2= g2[,1:2]

  n1= unique(c(rbind(g2, g1)[,1]))
  n2= unique(rbind(g2, g1)[,2])

  nn = sort(unique(c(n1,n2)))

  length(nn)
  dic1= rank(nn)
  names(dic1)=nn

  g1[,1] = dic1[g1[,1]]
  g1[,2] = dic1[g1[,2]]

  g2[,1] = dic1[g2[,1]]
  g2[,2] = dic1[g2[,2]]

  simil= delta_con(graph1 = g1, graph2 = g2,nnodes = length(nn), debug = F,
                   symmetrical = F)

  return(simil)
}



