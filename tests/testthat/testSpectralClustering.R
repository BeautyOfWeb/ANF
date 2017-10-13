library(ANF)
context("Spectral Clustering")

library(ExperimentHub)
eh <- ExperimentHub()
myfiles <- query(eh, "HarmonizedTCGAData")
Wall <- myfiles[[1]]
project_ids <- myfiles[[2]]
surv.plot <- myfiles[[3]]


test_that("spectral_clustering", {
    affinity.mat <- Wall[["adrenal_gland"]][["log.sel"]][["fpkm"]]
    true.disease.types <- as.factor(project_ids[rownames(affinity.mat)])
    
    labels <- spectral_clustering(affinity.mat, k = 2)
    expect_equal(as.numeric(table(labels, true.disease.types)), c(0, 76, 176, 1))
    
    labels <- spectral_clustering(affinity.mat, k = 3)
    expect_equal(as.numeric(table(labels, true.disease.types)), 
                 c(0, 0, 76, 155, 21, 1))
    
    labels <- spectral_clustering(affinity.mat, k = 2, type = 'sym')
    expect_equal(as.numeric(table(labels, true.disease.types)), c(0, 76, 176, 1))
    
    labels <- spectral_clustering(affinity.mat, k = 2, type = 'unnormalized')
    expect_equal(as.numeric(table(labels, true.disease.types)), c(0, 76, 176, 1))
    
})



test_that("eval_clu", {
    affinity.mat <- Wall[["adrenal_gland"]][["log.sel"]][["fpkm"]]
    true.disease.types <- as.factor(project_ids[rownames(affinity.mat)])
    
    res <- eval_clu(project_ids, w = affinity.mat, surv = surv.plot, verbose = FALSE)
    expect_equal(as.numeric(table(res$labels, true.disease.types)), c(0, 76, 176, 1))
    expect_equal(res$w, affinity.mat)
    expect_equal(as.numeric(res$clu.res), c(0.9628813, 0.9838114, 6.8336591))
    
    affinity.mat <- Wall[["uterus"]][["raw.all"]][["fpkm"]]
    true.disease.types <- as.factor(project_ids[rownames(affinity.mat)])
    res <- eval_clu(project_ids, w = affinity.mat, verbose = FALSE)
    expect_equal(as.numeric(table(res$labels, true.disease.types)), c(153, 268, 3, 51))
    expect_equal(as.numeric(res$clu.res), c(0.05633316, -0.05506271, NA))
})

