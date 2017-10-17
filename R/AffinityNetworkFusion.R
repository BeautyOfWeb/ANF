#' Calculate k-nearest-neighbor graph from affinity matrix
#' and normalize it as transition matrix
#'
#' @param W affinity matrix (its elements are non-negative real numbers)
#' @param K the number of k nearest neighbors
#'
#' @return a transition matrix of the same shape as W
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
        warning("K is too big: n-K-1 < 1. Assuming n > 2, 
          set K = floor(n/2)")
        K = floor(n / 2)
    }
    idx = t(apply(W, 1, order)[1:(n - K - 1),])
    for (i in seq_len(n)) {
        W[i, idx[i,]] = 0
    }
    # row normalization --> transition matrix
    S = W / rowSums(W)
    return(S)
}



#' Fuse affinity networks (i.e., matrices) through one-step or two-step random walk
#'
#' @param Wall a list of affinity matrices of the same shape.
#' @param K the number of k nearest neighbors for function kNN_graph
#' @param weight a list of non-negative real numbers (which will be 
#' normalized internally so that it sums to 1) that one-to-one correspond to 
#' the affinity matrices included in `Wall`. If not set, internally uniform weights 
#' are assigned to all affinity matrices in `Wall`.
#' @param type choose one of the two options: perform "one-step" random walk, 
#' or "two-step" random walk on the list of affinity matrices in `Wall`` to 
#' generate a fused affinity matrix. Default: "two-step" random walk
#' @param alpha a list of eight non-negative real numbers (which will be normalized
#' internally to make it sums to 1). Only used when "two-step" (default value 
#' of `type`) random walk is used. `alpha` is the weights for eight terms in the 
#' "two-step" random walk formula (check research paper for more explanations about 
#' the terms). Default value: (1, 1, 0, 0, 0, 0, 0, 0), i.e., only use the 
#' first two terms (since they are most effective in practice).
#' @param verbose logical(1); if true, print some information
#' 
#' @import Biobase
#' @import RColorBrewer
#' 
#' @return a fused transition matrix (representing a fused network)
#' @export
#' @examples
#' D1 = matrix(runif(400), nrow=20)
#' W1 = affinity_matrix(D1, 5)
#' D2 = matrix(runif(400), nrow=20)
#' W2 = affinity_matrix(D1, 5)
#' W = ANF(list(W1, W2), K=10)
ANF = function(Wall,
               K = 20,
               weight = NULL,
               type = c("two-step", "one-step"),
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
            message('length(Wall) = 1, return the only element.')
        }
        return(Wall[[1]])
    }
    
    # <UnhandledAssertion>
    # 1. Elements of Wall are matrices of the same dimension
    # 2. The rows and columns are aligned across all matrices
    
    # row normalization: similarity matrix --> transition matrix
    Wall <- lapply(Wall, function(elt) elt/rowSums(elt))
    
    # N: number of objects (e.g., patients)
    N = nrow(Wall[[1]])
    
    # normalize alpha to make it sum to 1
    alpha = alpha / sum(alpha)
    
    if (is.null(weight)) {
        # assign uniform weight
        weight = rep(1 / n, n)
    }
    if (length(weight) > n) {
        warning("length(weight) > length(Wall).
          Discard the last Ws")
        weight = weight[1:n]
    }
    if (length(weight) < n) {
        warning("length(weight) < length(Wall).
          Assume the last Ws have zero weights")
        weight = c(weight, rep(0, n - length(weight)))
    }
    if (sum(weight != 0) == 1) {
        message("Only one view has a non-zero weight")
        return(Wall[[which(weight != 0)]])
    } else {
        # normalize weight to make it sum to 1
        if (sum(weight) != 0) {
            weight = weight / sum(weight)
        } else {
            stop("sum(weight)==0")
        }
    }
    

    # Construct kNN graphs
    Sall <- lapply(Wall, function(elt) kNN_graph(elt, K))
    
    type <- match.arg(type)
    switch(type,
           "one-step" = type <- 1,
           "two-step" = type <- 2)
    
    if (type == 1) {
        # "one-step" random walk
        W = matrix(0, N, N)
        for (j in seq_len(n)) {
            # return a weighted average of kNN affinity matrices
            W = W + weight[j] * Sall[[j]]
        }
    } else if (type == 2) {
        # "two-step" random walk
        for (j in seq_len(n)) {
            W = matrix(0, N, N)
            S = matrix(0, N, N)
            sumWeight = sum(weight) - weight[j]
            if (sumWeight == 0) {
                # This could happen only when input parameter `weight` contains
                # both positive and negative numbers
                warning("sum(weight) - weight[j] == 0.
              The result is probably unreliable.")
                sumWeight = sum(weight)
            }
            for (k in seq_len(n)) {
                if (k != j) {
                    # complementary view
                    W = W + (weight[k] / sumWeight) * Wall[[k]]
                    S = S + (weight[k] / sumWeight) * Sall[[k]]
                }
            }
            
            # The eight terms for "two-step" random walk
            Wall[[j]] = alpha[1] * Sall[[j]] %*% S + alpha[2] * S %*% Sall[[j]] +
                alpha[3] * Sall[[j]] %*% W + alpha[4] * W %*% Sall[[j]] +
                alpha[5] * Wall[[j]] %*% S  + alpha[6] * S %*% Wall[[j]] +
                alpha[7] * Wall[[j]] %*% W + alpha[8] * W %*% Wall[[j]]
        }
    } else {
        stop("type other than 1 or 2 is undefined")
    }
    W = matrix(0, N, N)
    for (j in seq_len(n)) {
        W = W + weight[j] * Wall[[j]]
    }
    rownames(W) = rownames(Wall[[1]])
    colnames(W) = colnames(Wall[[1]])
    return(W)
}