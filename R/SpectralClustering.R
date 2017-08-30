#' Generate a symmetric affinity matrix based on distance matrix 
#' using 'local' Gaussian kernel
#'
#' @param D distance matrix
#' @param k the number of k-nearest neighbors
#' @param alpha coefficient for local diameters
#' @param beta coefficient for pair-wise distance 
#'
#' @return an affinity matrix 
#' @export
#'
#' @examples
#' D = matrix(runif(400),20)
#' A = affinity_matrix(D, 5)
affinity_matrix <- function(D, k, alpha=1/6, beta=1/6)
{
    n = nrow(D)
    # make it symmetric in case it is not
    D = (D + t(D)) / 2
    diag(D) = 0
    D[D<0] = 0
    finiteMean = function(x){
        return(mean(x[is.finite(x)]))
    }
    # <Tricky> k nearest neighbors should exclude itself: 1:k+1
    # <?> here we divide by k+1, maybe k is better
    d = apply(apply(D, 2, sort)[1:k+1,], 2, finiteMean)
    sigma = alpha * (outer(d, d, '+')) + beta * D + .Machine$double.eps
    return(stats::dnorm(D, 0, sigma))
}
#' Finding optimal discrete solutions for spectral clustering
#' 
#' @references Stella, X. Yu, and Jianbo Shi. "Multiclass spectral clustering." ICCV. IEEE, 2003.
#' @param Y first k eigenvectors of graph Laplacian corresponding to k smallest eigenvalues
#' @param verbose Boolean; if true, print some information
#'
#' @return class assignment matrix (0-1)
#' @export
#'
#' @examples
#' D = matrix(runif(400),20)
#' A = affinity_matrix(D, 5)
#' d = rowSums(A)
#' L = diag(d) - A
#' NL = diag(1/d) %*% L
#' e = eigen(NL)
#' Y = Re(e$vectors[,-1:-17])
#' X = pod(Y)
pod = function(Y, verbose=FALSE) {
    if(!is.matrix(Y)) {
        stop("Y should be a matrix")
    } else {
        # normalize in case the input is not
        Y = Y / sqrt(rowSums(Y^2))
    }
    N = nrow(Y)
    K = ncol(Y)
    if ( N<K) {
        stop("Y should have more rows than columns")
    }
    R = matrix(0, K, K)
    # chose any row of Y, here the middle row
    R[,1] = Y[round(N/2),]
    c = rep(0, N)
    for (i in 2:K) {
        c = c + abs(Y %*% R[,i-1])
        R[,i] = Y[which.min(c),]
    }
    e = 0
    error = 1
    while (error > .Machine$double.eps) {
        idx = apply(Y %*% R, 1, which.max)
        X = matrix(0, N, K)
        X[cbind(1:N, idx)] = 1
        s = svd(t(X)%*%Y)
        R = s$v %*% t(s$u)
        error = abs(sum(s$d)-e)
        e = sum(s$d)
    } 
    if(verbose) {
        print(paste("e =", e))
    }
    return(X)
}
#' spectral_clustering
#'
#' @param A affinity matrix
#' @param k the number of clusters
#' @param type choice of normalized graph Laplacian
#' @param verbose Boolean; if true, print user-friendly information
#'
#' @return a numeric vector as class labels
#' @export
#'
#' @examples
#' D = matrix(runif(400),20)
#' A = affinity_matrix(D, 5)
#' labels = spectral_clustering(A, 2)
spectral_clustering <- function(A, k, type=2, verbose=FALSE)
{
    n = nrow(A)
    d = rowSums(A)
    d[d == 0] = .Machine$double.eps
    L = diag(d) - A
    if(type==1) {
        NL = L
    } else if(type==2) {
        NL = diag(1/d) %*% L
    } else if(type == 3) {
        NL = diag(1/sqrt(d)) %*% L %*% diag(1/sqrt(d))
    } else {
        stop("type should be 1, 2, or 3")
    }
    e = eigen(NL)
    if(verbose) {
        print(paste("the k+3 smallest eigen values:", e$values[-1:-(n-k-3)]))
    }
    Y = Re(e$vectors[,-1:-(n-k)])
    if (type == 3) {
        Y = Y / sqrt(rowSums(Y ^ 2))
    }
    labels = apply(pod(Y),1,which.max)
    return(labels)
}

#' Evaluate cluster result 
#'
#' @param true_class A named vector of true class labels
#' @param w affinity matrix
#' @param d distance matrix; 
#' if w is NULL, calcuate w using d
#' @param k an integer; if w is null, w = affinity_matrix(d, k)
#' @param num_clu an integer; number of clusters; 
#' if NULL, set num_clu to be the number of classes using true_class
#' @param surv a data.frame with at least two columns: 
#' time (days_to_death or days_to_last_follow_up), and censored (Boolean)
#' @param type_L an integer of 1, 2, or 3; choice of normalized graph Laplacian
#' @param verbose Boolean; if true, print some information
#'
#' @return a vector of Normalized Mutual Informatin (NMI), Adjusted Rand Index (ARI),
#' and -log10(p) of log rank test of survival distributions of patient clusters
#' @import igraph MASS Biobase
#' @export
#' @examples 
#' true.class = rep(c(1,2),each=100)
#' feature.mat1 = MASS::mvrnorm(100, rep(0, 20), diag(runif(20,0.2,2)))
#' feature.mat2 = MASS::mvrnorm(100, rep(0.5, 20), diag(runif(20,0.2,2)))
#' feature1 = rbind(feature.mat1, feature.mat2)
#' d = dist(feature1)
#' d = as.matrix(d)
#' A = affinity_matrix(d, 10)
#' res = eval_clu(true_class=true.class, w=A)
eval_clu = function(true_class, w= NULL, d = NULL, 
                    k = 10, num_clu = NULL, surv=NULL, 
                    type_L=2, verbose=TRUE) {
    if(is.null(w) & is.null(d)){
        stop("Either distance or affinity matrix must be provided")
    }
    if(is.null(w) & !is.null(d)) {
        w = affinity_matrix(d, min(k, nrow(d)-1))
    }
    # <InputRequirement> 
    if (is.null(names(true_class))) {
        names(true_class) = rownames(w)
    }
    if(!all(rownames(w) %in% names(true_class))) {
        stop("some rownames(w) is unknown to true_class")
    }
    true_class = as.factor(true_class[rownames(w)])
    if(is.null(num_clu)) {
        num_clu = length(unique(true_class))
    }
    labels = spectral_clustering(w, num_clu, type = type_L)
    nmi = igraph::compare(true_class, labels, method='nmi')
    adjusted.rand = igraph::compare(true_class, labels, method='adjusted.rand')
    
    if(!is.null(surv) & 
       length(intersect(rownames(w), rownames(surv))) > nrow(w)/2) {
        surv = surv[rownames(w),]
        f = survival::Surv(surv$time, !surv$censored)
        fit = survival::survdiff(f~labels)
        pval = stats::pchisq(fit$chisq, df=length(fit$n)-1, lower.tail = FALSE)
    } else {
        pval = NA
    }
    if(verbose) {
        print(paste('nmi:', nmi))
        print(paste('adjusted.rand:', adjusted.rand))
        print(table(true_class, labels))
        print(paste('survdiff p value:', pval))
    }
    
    clu.res = c(nmi, adjusted.rand, -log10(pval))
    names(clu.res) = c("NMI", "ARI", "-log10(p)")
    return(list(w=w, clu.res=clu.res, labels=labels))
}