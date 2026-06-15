#' Read a candidate table
#'
#' @param path Input table path.
#' @param sep Field separator. G4Hunter-style files often use tabs even when
#'   named `.csv`.
#' @param ... Additional arguments passed to [utils::read.table()].
#'
#' @return A data.frame.
#'
#' @examples
#' path <- system.file("extdata", "g4_candidates.tsv", package = "AF3Input")
#' candidates <- read_candidate_table(path)
#' names(candidates)
#' @export
read_candidate_table <- function(path, sep = "\t", ...) {
    utils::read.table(
        file = path,
        header = TRUE,
        sep = sep,
        quote = "\"",
        comment.char = "",
        stringsAsFactors = FALSE,
        check.names = FALSE,
        ...
    )
}

#' Write one AlphaFold 3 JSON file
#'
#' @param job AlphaFold 3 job list.
#' @param path Output JSON path.
#' @param pretty Whether to pretty-print JSON.
#'
#' @return Invisibly returns `path`.
#'
#' @examples
#' job <- af3_job("example", "GGGTTAGGGTTAGGGTTAGGG")
#' path <- tempfile(fileext = ".json")
#' write_af3_json(job, path)
#' file.exists(path)
#' @export
write_af3_json <- function(job, path, pretty = TRUE) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    jsonlite::write_json(
        job,
        path = path,
        pretty = pretty,
        auto_unbox = TRUE,
        null = "null"
    )
    invisible(path)
}

#' Write AlphaFold 3 JSON files from a candidate table
#'
#' @param table Data.frame containing candidate sequences.
#' @param out_dir Output directory.
#' @param sequence_col Name of the column containing DNA sequence.
#' @param description_cols Column names to include in DNA entity description.
#' @param name_cols Column names used to build job/file names.
#' @param score_col Optional score column used for reverse-complementing
#'   negative-score rows.
#' @param reverse_complement_negative If `TRUE`, reverse-complement sequences
#'   whose `score_col` value is negative.
#' @param filter Optional function taking `table` and returning a logical
#'   vector.
#' @param limit Optional maximum number of rows to write after filtering.
#' @param dna_copies Number of identical DNA copies.
#' @param ions Optional list of ion specs with `code` and `copies`.
#' @param seeds Integer model seeds.
#' @param add_complement Whether to add a complementary DNA entity.
#' @param modifications Optional data.frame of modifications applied to each
#'   generated main DNA entity.
#' @param skip_invalid Whether to skip rows with invalid DNA sequence
#'   characters.
#'
#' @return Character vector of written JSON paths.
#'
#' @examples
#' candidates <- data.frame(
#'     id = "candidate_1",
#'     score = -1.5,
#'     sequence = "CCCTAACCCTAACCCTAACCCTAACCC"
#' )
#' paths <- write_af3_jsons_from_table(
#'     candidates,
#'     out_dir = tempdir(),
#'     sequence_col = "sequence",
#'     description_cols = c("id", "score"),
#'     name_cols = "id",
#'     score_col = "score",
#'     reverse_complement_negative = TRUE,
#'     ions = list(list(code = "K", copies = 1))
#' )
#' file.exists(paths)
#' @export
write_af3_jsons_from_table <- function(table, out_dir, sequence_col,
        description_cols = character(), name_cols = NULL, score_col = NULL,
        reverse_complement_negative = FALSE, filter = NULL, limit = NULL,
        dna_copies = 1L, ions = NULL, seeds = 1L,
        add_complement = FALSE, modifications = NULL,
        skip_invalid = TRUE) {
    options <- list(
        sequence_col = sequence_col,
        description_cols = description_cols,
        name_cols = name_cols,
        score_col = score_col,
        reverse_complement_negative = reverse_complement_negative,
        filter = filter,
        limit = limit,
        dna_copies = dna_copies,
        ions = ions,
        seeds = seeds,
        add_complement = add_complement,
        modifications = modifications,
        skip_invalid = skip_invalid
    )
    table <- prepare_candidate_table(table, options)

    paths <- vector("list", nrow(table))
    for (i in seq_len(nrow(table))) {
        paths[[i]] <- write_af3_table_row(
            row = table[i, , drop = FALSE],
            row_index = i,
            out_dir = out_dir,
            options = options
        )
    }
    unlist(paths, use.names = FALSE)
}

prepare_candidate_table <- function(table, options) {
    stopifnot(is.data.frame(table))
    require_columns(table, options$sequence_col)
    require_columns(table, options$description_cols)
    require_columns(table, options$name_cols)
    require_columns(table, options$score_col)

    table <- filter_candidate_table(table, options$filter)
    if (!is.null(options$limit)) {
        table <- utils::head(table, options$limit)
    }
    table
}

filter_candidate_table <- function(table, filter) {
    if (is.null(filter)) {
        return(table)
    }
    keep <- filter(table)
    if (!is.logical(keep) || length(keep) != nrow(table)) {
        stop(
            "'filter' must return one logical value per table row.",
            call. = FALSE
        )
    }
    table[keep, , drop = FALSE]
}

write_af3_table_row <- function(row, row_index, out_dir, options) {
    sequence <- row_sequence(row, options)
    if (!is_valid_dna(sequence$value)) {
        if (options$skip_invalid) {
            return(NULL)
        }
        stop("Invalid DNA sequence in row ", row_index, ".", call. = FALSE)
    }

    description <- row_description(
        row,
        options$description_cols,
        sequence$reverse_applied
    )
    job_name <- row_name(row, row_index, options$name_cols)
    job <- af3_job(
        name = job_name,
        dna_sequence = sequence$value,
        dna_copies = options$dna_copies,
        seeds = options$seeds,
        ions = options$ions,
        add_complement = options$add_complement,
        modifications = options$modifications,
        description = description
    )
    path <- file.path(out_dir, paste0(sanitize_filename(job_name), ".json"))
    write_af3_json(job, path)
    path
}

row_sequence <- function(row, options) {
    sequence <- toupper(as.character(row[[options$sequence_col]]))
    reverse_applied <- FALSE
    if (options$reverse_complement_negative && !is.null(options$score_col)) {
        score <- parse_numeric_score(row[[options$score_col]])
        if (!is.na(score) && score < 0) {
            sequence <- reverse_complement_dna(sequence)
            reverse_applied <- TRUE
        }
    }
    list(value = sequence, reverse_applied = reverse_applied)
}

is_valid_dna <- function(sequence) {
    tryCatch({
        validate_dna_sequence(sequence)
        TRUE
    }, error = function(e) FALSE)
}

require_columns <- function(table, columns) {
    if (is.null(columns) || length(columns) == 0L) {
        return(invisible(TRUE))
    }
    missing <- setdiff(columns, names(table))
    if (length(missing) > 0L) {
        stop(
            "Missing column(s): ",
            paste(missing, collapse = ", "),
            call. = FALSE
        )
    }
    invisible(TRUE)
}

row_description <- function(row, columns, reverse_applied) {
    parts <- character()
    for (column in columns) {
        parts <- c(parts, paste0(column, "=", as.character(row[[column]])))
    }
    parts <- c(parts, paste0("reverse_complemented=", reverse_applied))
    paste(parts, collapse = "; ")
}

row_name <- function(row, index, name_cols) {
    if (is.null(name_cols) || length(name_cols) == 0L) {
        return(paste0("af3_job_", index))
    }
    paste(
        vapply(
            name_cols,
            function(column) as.character(row[[column]]),
            character(1L)
        ),
        collapse = "_"
    )
}

sanitize_filename <- function(x) {
    x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
    x <- gsub("^_+|_+$", "", x)
    if (nzchar(x)) x else "af3_job"
}

parse_numeric_score <- function(x) {
    if (is.numeric(x)) {
        return(as.numeric(x))
    }
    x <- as.character(x)
    if (!grepl("^-?[0-9]+(\\.[0-9]+)?([eE][-+]?[0-9]+)?$", x)) {
        return(NA_real_)
    }
    as.numeric(x)
}
