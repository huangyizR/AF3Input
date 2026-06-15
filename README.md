# AF3Input

AF3Input creates AlphaFold 3 input JSON files from tabular DNA candidate calls.
It is designed first for G-quadruplex candidate tables produced by R workflows,
while keeping the sequence column and metadata columns configurable.

## Example

```r
library(AF3Input)

candidates <- read_candidate_table(
    "chr1_test_25_1.5.csv",
    sep = "\t"
)

write_af3_jsons_from_table(
    candidates,
    out_dir = "af3_json",
    sequence_col = "SEQUENCE",
    description_cols = c("ID", "POSITION", "LENGTH", "SCORE", "ABS_SCORE"),
    score_col = "SCORE",
    reverse_complement_negative = TRUE,
    ions = list(list(code = "K", copies = 3)),
    seeds = c(1, 2, 3),
    limit = 10
)
```

Ions are written as AlphaFold 3 ligand entities using CCD codes.
