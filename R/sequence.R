#' Reverse-complement a DNA sequence
#'
#' @param sequence Character scalar containing A, C, G, T, or IUPAC
#'   ambiguous DNA bases.
#'
#' @return A character scalar with the reverse complement.
#'
#' @examples
#' reverse_complement_dna("CCCTAACCCTAACCCTAACCCTAACCC")
#' @export
reverse_complement_dna <- function(sequence) {
    stopifnot(is.character(sequence), length(sequence) == 1L)
    sequence <- validate_dna_sequence(sequence, allow_ambiguous = TRUE)
    as.character(
        Biostrings::reverseComplement(Biostrings::DNAString(sequence))
    )
}

validate_dna_sequence <- function(sequence, allow_ambiguous = FALSE) {
    if (!is.character(sequence) ||
            length(sequence) != 1L ||
            is.na(sequence)) {
        stop(
            "'sequence' must be a non-missing character scalar.",
            call. = FALSE
        )
    }
    sequence <- toupper(sequence)
    pattern <- if (allow_ambiguous) "^[ACGTRYSWKMBDHVN]+$" else "^[ACGT]+$"
    if (!grepl(pattern, sequence)) {
        allowed <- if (allow_ambiguous) "IUPAC DNA symbols." else "A/C/G/T."
        stop("'sequence' contains bases outside ", allowed, call. = FALSE)
    }
    sequence
}

chain_ids <- function(n, start = 1L) {
    if (n < 1L) {
        character()
    } else {
        vapply(seq.int(start, start + n - 1L), spreadsheet_id, character(1L))
    }
}

spreadsheet_id <- function(index) {
    letters <- LETTERS
    out <- character()
    while (index > 0L) {
        index <- index - 1L
        out <- c(letters[(index %% 26L) + 1L], out)
        index <- index %/% 26L
    }
    paste(out, collapse = "")
}
