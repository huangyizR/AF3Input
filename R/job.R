#' Create an AlphaFold 3 job object
#'
#' @param name Job name.
#' @param dna_sequence Main DNA sequence.
#' @param dna_copies Number of identical DNA copies.
#' @param seeds Integer model seeds.
#' @param ions Optional list of ion specs. Each ion spec is a list with `code`
#'   and `copies`.
#' @param add_complement Whether to add a complementary DNA entity.
#' @param complement_type Either `"reverse_complement"` or `"complement"`.
#' @param modifications Optional DNA modifications for the main DNA entity.
#' @param description Optional main DNA description.
#' @param version AlphaFold 3 input JSON version.
#'
#' @return A list ready to serialize as AlphaFold 3 JSON.
#'
#' @examples
#' af3_job(
#'     name = "g4_example",
#'     dna_sequence = "GGGTTAGGGTTAGGGTTAGGG",
#'     ions = list(list(code = "K", copies = 2)),
#'     seeds = c(1, 2)
#' )
#' @export
af3_job <- function(name, dna_sequence, dna_copies = 1L, seeds = 1L,
        ions = NULL, add_complement = FALSE,
        complement_type = c("reverse_complement", "complement"),
        modifications = NULL, description = NULL, version = 4L) {
    complement_type <- match.arg(complement_type)
    dna_sequence <- validate_dna_sequence(dna_sequence)
    dna_copies <- validate_positive_integer(dna_copies, "'dna_copies'")
    seeds <- validate_seeds(seeds)

    state <- initial_dna_state(
        dna_sequence = dna_sequence,
        dna_copies = dna_copies,
        modifications = modifications,
        description = description
    )

    if (add_complement) {
        state <- append_complement_entity(
            state,
            dna_sequence,
            complement_type
        )
    }
    state <- append_ion_entities(state, ions)

    list(
        name = as.character(name),
        modelSeeds = as.list(seeds),
        sequences = state$sequences,
        dialect = "alphafold3",
        version = as.integer(version)
    )
}

validate_positive_integer <- function(value, name) {
    value <- as.integer(value)
    if (length(value) != 1L || is.na(value) || value < 1L) {
        stop(name, " must be a positive integer.", call. = FALSE)
    }
    value
}

validate_seeds <- function(seeds) {
    seeds <- as.integer(seeds)
    if (length(seeds) < 1L || anyNA(seeds)) {
        stop("'seeds' must contain at least one integer seed.", call. = FALSE)
    }
    seeds
}

initial_dna_state <- function(dna_sequence, dna_copies, modifications,
        description) {
    dna_ids <- chain_ids(dna_copies)
    sequences <- list(af3_dna_entity(
        id = single_or_list_id(dna_ids),
        sequence = dna_sequence,
        modifications = modifications,
        description = description
    ))
    list(sequences = sequences, next_id = dna_copies + 1L)
}

append_complement_entity <- function(state, dna_sequence, complement_type) {
    comp <- complement_dna(
        dna_sequence,
        reverse = complement_type == "reverse_complement"
    )
    state$sequences[[length(state$sequences) + 1L]] <- af3_dna_entity(
        id = chain_ids(1L, state$next_id),
        sequence = comp,
        description = paste(complement_type, "of main DNA")
    )
    state$next_id <- state$next_id + 1L
    state
}

append_ion_entities <- function(state, ions) {
    for (ion in ions %||% list()) {
        copies <- validate_positive_integer(
            ion$copies %||% 1L,
            "Ion copies"
        )
        ion_ids <- chain_ids(copies, state$next_id)
        state$sequences[[length(state$sequences) + 1L]] <- af3_ion_entity(
            code = ion$code,
            id = single_or_list_id(ion_ids),
            description = paste(copies, toupper(ion$code), "ion copy/copies")
        )
        state$next_id <- state$next_id + copies
    }
    state
}

single_or_list_id <- function(ids) {
    if (length(ids) == 1L) ids else as.list(ids)
}

complement_dna <- function(sequence, reverse = TRUE) {
    sequence <- validate_dna_sequence(sequence)
    if (reverse) {
        return(reverse_complement_dna(sequence))
    }
    chartr("ACGT", "TGCA", sequence)
}

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}
