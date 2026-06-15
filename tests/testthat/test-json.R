test_that("AF3 job has expected AlphaFold 3 fields", {
    job <- af3_job(
        name = "example",
        dna_sequence = "GGGTTAGGGTTAGGGTTAGGG",
        ions = list(list(code = "K", copies = 2)),
        seeds = c(1, 2)
    )
    expect_equal(job$dialect, "alphafold3")
    expect_equal(job$version, 4L)
    expect_length(job$sequences, 2L)
    expect_equal(job$sequences[[2L]]$ligand$ccdCodes[[1L]], "K")
})

test_that("table writer allows configurable sequence and description columns", {
    x <- data.frame(
        candidate = "g4_1",
        score = -1.5,
        seq = "CCCTAACCCTAACCCTAACCCTAACCC",
        stringsAsFactors = FALSE
    )
    out_dir <- file.path(tempdir(), "af3input-test")
    paths <- write_af3_jsons_from_table(
        x,
        out_dir = out_dir,
        sequence_col = "seq",
        description_cols = c("candidate", "score"),
        name_cols = "candidate",
        score_col = "score",
        reverse_complement_negative = TRUE,
        ions = list(list(code = "K", copies = 1))
    )
    expect_length(paths, 1L)
    parsed <- jsonlite::read_json(paths[[1L]])
    expect_equal(parsed$sequences[[1L]]$dna$sequence,
                 "GGGTTAGGGTTAGGGTTAGGGTTAGGG")
    expect_match(parsed$sequences[[1L]]$dna$description,
                 "candidate=g4_1")
})
