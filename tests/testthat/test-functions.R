library(DOAGDC)

context("Test Download and Concatenate functions")

test_that("Parse metadata automatically", {

    # testDownload function
    download_gdc("gene", "CHOL", "legacy", work_dir = "~/Desktop")

    path_2_check <- "~/Desktop/DOAGDC/CHOL/gene_data/"

    # test1
    testthat::expect_true(dir.exists(path_2_check))

    # test2
    testthat::expect_equal(length(dir(path_2_check)), 92)

    # testConcatenate function
    concatenate_files("gene",
                        name = "HIF3A",
                        data_base = "legacy",
                        tumor = "CHOL",
                        work_dir = "~/Desktop")

    # test3
    data <- CHOL_LEGACY_gene_tumor_data$gene_tumor_normalized
    testthat::expect_equal(dim(data), c(20531, 36))

    # test4
    data <- CHOL_LEGACY_gene_tumor_data$gene_tumor_normalized_selected_HIF3A
    testthat::expect_equal(dim(data), c(36, 1))
})

