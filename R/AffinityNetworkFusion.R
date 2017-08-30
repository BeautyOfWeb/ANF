#' Calculate k-nearest-neighbor graph from affinity matrix
#' normalize it as transition matrix
#'
#' @param W affinity matrix
#' @param K the number of k nearest neighbors
#'
#' @return an transition matrix of the same shape as W
#' @export
#'
#' @examples
#' D = matrix(runif(400),20)
#' W = affinity_matrix(D, 5)
#' S = kNN_graph(W, 5)
kNN_graph = function(W, K) {
    n = nrow(W)
    # k-nearest neighbors does not include itself
    # set similarity = 0 for n-K-1 farest neighbors
    if (n - K - 1 < 1) {
        print("Error: K is too big: n-K-1 < 1. Assuming n > 2, 
          set K = floor(n/2)")
        K = floor(n / 2)
    }
    idx = t(apply(W, 1, order)[1:(n - K - 1),])
    for (i in 1:n) {
        W[i, idx[i,]] = 0
    }
    # row normalization --> transition matrix
    S = W / rowSums(W)
    return(S)
}

#' Fuse affinity networks through one-step or two-step random walk
#'
#' @param Wall a list of affinity matrices of the same shape
#' @param K the number of k nearest neighbors for function kNN_graph
#' @param weight the weight for each view
#' @param type 1: one-step random walk; 2: two-step random walk
#' @param alpha used as weights for terms in type 2 (two-step random walk)
#' @param verbose Boolean; if true, print some information
#'
#' @return a fused transition matrix (representing fused network)
#' @export
#' @examples
#' D1 = matrix(runif(400),20)
#' W1 = affinity_matrix(D1, 5)
#' D2 = matrix(runif(400),20)
#' W2 = affinity_matrix(D1, 5)
#' W = ANF(list(W1, W2), K=10)
ANF = function(Wall,
               K = 20,
               weight = NULL,
               type = 2,
               alpha = c(1, 1, 0, 0, 0, 0, 0, 0),
               verbose = FALSE) {
    # sanity check Wall
    if (!is.list(Wall)) {
        stop("Wall must be a list of matrices")
    }
    # n: number of views
    n = length(Wall)
    if (n == 0) {
        stop("Wall is empty")
    }
    if (n == 1) {
        if (verbose) {
            print('length(Wall) = 1, return the only element.')
        }
        return(Wall[[1]])
    }
    # <UnhandledAssertion>
    # 1. Elements of Wall are matrices of the same dimension
    # 2. The rows and columns are aligned across all matrices
    for (i in 1:n) {
        # row normalization: similarity matrix --> transition matrix
        Wall[[i]] = Wall[[i]] / rowSums(Wall[[i]])
    }
    
    # N: number of points (objects)
    N = nrow(Wall[[1]])
    # normalize alpha
    alpha = alpha / sum(alpha)
    
    # <?> should I ensure non-negative weights?
    if (is.null(weight)) {
        # assign uniform weight
        weight = rep(1 / n, n)
    }
    if (length(weight) > n) {
        print("WARNING: length(weight) > length(Wall).
          Discard the last Ws")
        weight = weight[1:n]
    }
    if (length(weight) < n) {
        print("WARNING: length(weight) < length(Wall).
          Assume the last Ws have 0 weight")
        weight = c(weight, rep(0, n - length(weight)))
    }
    if (sum(weight != 0) == 1) {
        print("only one view has a non-zero weight")
        return(Wall[[which(weight != 0)]])
    } else {
        # normalize weight
        if (sum(weight) != 0) {
            weight = weight / sum(weight)
        } else {
            stop("sum(weight)==0")
        }
    }
    
    # construct k-nearest-neighbor graph
    # <?> Should it change during iterations?
    Sall = list()
    for (j in 1:n) {
        Sall[[j]] = kNN_graph(Wall[[j]], K)
    }
    
    if (type == 1) {
        W = matrix(0, N, N)
        for (j in 1:j) {
            # simplified version (does not include Wall[[j]])
            W = W + weight[j] * Sall[[j]]
        }
    } else if (type == 2) {
        for (j in 1:n) {
            W = matrix(0, N, N)
            S = matrix(0, N, N)
            sumWeight = sum(weight) - weight[j]
            if (sumWeight == 0) {
                print("WARNING: sum(weight) - weight[j] == 0.
              The result is probably unreliable.")
                sumWeight = sum(weight)
            }
            for (k in 1:n) {
                if (k != j) {
                    # complementary view
                    W = W + (weight[k] / sumWeight) * Wall[[k]]
                    S = S + (weight[k] / sumWeight) * Sall[[k]]
                }
            }
            
            Wall[[j]] = alpha[1] * Sall[[j]] %*% S + alpha[2] * S %*% Sall[[j]] +
                alpha[3] * Sall[[j]] %*% W + alpha[4] * W %*% Sall[[j]] +
                alpha[5] * Wall[[j]] %*% S  + alpha[6] * S %*% Wall[[j]] +
                alpha[7] * Wall[[j]] %*% W + alpha[8] * W %*% Wall[[j]]
        }
    } else {
        stop("type other than 1 or 2 is undefined")
    }
    W = matrix(0, N, N)
    for (j in 1:n) {
        W = W + weight[j] * Wall[[j]]
    }
    rownames(W) = rownames(Wall[[1]])
    colnames(W) = colnames(Wall[[1]])
    return(W)
}