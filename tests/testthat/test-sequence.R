test_that("reverse complement works for DNA", {
    expect_equal(
        reverse_complement_dna("CCCTAACCCTAACCCTAACCCTAACCC"),
        "GGGTTAGGGTTAGGGTTAGGGTTAGGG"
    )
})

test_that("invalid DNA is rejected by default", {
    expect_error(af3_dna_entity("A", "ACMGT"), "A/C/G/T")
})
