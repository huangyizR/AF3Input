#' Create an AlphaFold 3 DNA entity
#'
#' @param id Chain identifier, or a character vector of identifiers for copies.
#' @param sequence DNA sequence.
#' @param modifications Optional data.frame with columns
#'   `modificationType` and `basePosition`, or compatible aliases `type` and
#'   `position`.
#' @param description Optional entity description.
#' @param allow_ambiguous Whether to allow IUPAC ambiguous DNA bases.
#'   AlphaFold 3 DNA input itself expects A/C/G/T, so the default is `FALSE`.
#'
#' @return A nested list representing one AlphaFold 3 DNA sequence entity.
#'
#' @examples
#' af3_dna_entity(
#'     id = "A",
#'     sequence = "GGGTTAGGGTTAGGGTTAGGG",
#'     description = "candidate G4-forming strand"
#' )
#' @export
af3_dna_entity <- function(id, sequence, modifications = NULL,
        description = NULL, allow_ambiguous = FALSE) {
    sequence <- validate_dna_sequence(
        sequence,
        allow_ambiguous = allow_ambiguous
    )
    dna <- list(id = id, sequence = sequence)
    mods <- normalize_dna_modifications(modifications)
    if (length(mods) > 0L) {
        validate_modification_positions(mods, sequence)
        dna$modifications <- mods
    }
    if (!is.null(description)) {
        dna$description <- as.character(description)
    }
    list(dna = dna)
}

#' Create an AlphaFold 3 ion entity
#'
#' @param code CCD code, for example `"K"`, `"NA"`, or `"MG"`.
#' @param id Ligand identifier, or vector of identifiers for copies.
#' @param description Optional entity description.
#'
#' @return A nested list representing one AlphaFold 3 ligand entity.
#'
#' @examples
#' af3_ion_entity(code = "K", id = "B", description = "potassium ion")
#' @export
af3_ion_entity <- function(code, id, description = NULL) {
    if (!is.character(code) || length(code) != 1L || is.na(code)) {
        stop("'code' must be a non-missing character scalar.", call. = FALSE)
    }
    ligand <- list(id = id, ccdCodes = list(toupper(code)))
    if (!is.null(description)) {
        ligand$description <- as.character(description)
    }
    list(ligand = ligand)
}

normalize_dna_modifications <- function(modifications) {
    if (is.null(modifications)) {
        return(list())
    }
    if (!is.data.frame(modifications)) {
        stop("'modifications' must be NULL or a data.frame.", call. = FALSE)
    }
    names(modifications) <- sub(
        "^type$",
        "modificationType",
        names(modifications)
    )
    names(modifications) <- sub(
        "^position$",
        "basePosition",
        names(modifications)
    )
    required <- c("modificationType", "basePosition")
    missing <- setdiff(required, names(modifications))
    if (length(missing) > 0L) {
        stop(
            "Missing modification columns: ",
            paste(missing, collapse = ", "),
            call. = FALSE
        )
    }
    lapply(seq_len(nrow(modifications)), function(i) {
        list(
            modificationType = as.character(
                modifications$modificationType[[i]]
            ),
            basePosition = as.integer(modifications$basePosition[[i]])
        )
    })
}

validate_modification_positions <- function(modifications, sequence) {
    n <- nchar(sequence)
    for (modification in modifications) {
        pos <- modification$basePosition
        if (is.na(pos) || pos < 1L || pos > n) {
            stop(
                "DNA modification basePosition is outside the sequence.",
                call. = FALSE
            )
        }
    }
}
