library(ANF)
context("Affinity Network Fusion")

library(ExperimentHub)
eh <- ExperimentHub()
myfiles <- query(eh, "HarmonizedTCGAData")
Wall <- myfiles[[1]]
project_ids <- myfiles[[2]]
surv.plot <- myfiles[[3]]

test_that("ANF", {
    fused.mat <- ANF(Wall = Wall$uterus$raw.all)
    res <- eval_clu(true_class = project_ids[rownames(fused.mat)], w = fused.mat, 
                    verbose = FALSE)
    expect_equal(as.numeric(res$clu.res), c(0.4852676, 0.6842040, NA) + 1e-8)
    
    fused.mat <- ANF(Wall = Wall$kidney$log.sel)
    res <- eval_clu(true_class = project_ids[rownames(fused.mat)], w = fused.mat, 
                    verbose = FALSE)
    expect_equal(as.numeric(res$clu.res), c(0.7877615, 0.8653626, NA) + 1e-8)
    
    fused.mat <- ANF(Wall = Wall$lung$normalized)
    res <- eval_clu(true_class = project_ids[rownames(fused.mat)], w = fused.mat, 
                    surv = surv.plot, verbose = FALSE)
    expect_equal(as.numeric(res$clu.res), c(0.7182888, 0.8121813, 0.07952455) - 2e-8)
    
    fused.mat <- ANF(Wall = Wall$colorectal$log.sel)
    res <- eval_clu(true_class = project_ids[rownames(fused.mat)], w = fused.mat, 
                    surv = surv.plot, verbose = FALSE)
    expect_equal(as.numeric(res$clu.res), c(0.05393634, -0.03659459, 0.37974029))
})
