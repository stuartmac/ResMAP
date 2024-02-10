01_mutlib-initial-processing
================
Stuart MacGowan
2024-02-10

``` r
# Read the datasets
paths <- c("../data/PfKRS-mutlib-final-calculations-wildtype_231117.tsv",
           "../data/PfKRS-mutlib-final-calculations-clado_231117.tsv",
           "../data/PfKRS-mutlib-final-calculations-706_231117.tsv")

# Read the datasets
mutlib <- paths |>
  tidytable::map_dfr(read_delim, delim = "\t", escape_double = FALSE,
          trim_ws = TRUE, col_types = cols("Name" = col_character(), "Sequence" = col_character(), "Count" = col_double()), .id = "dataset")

# Parse Name column to get codon and position
# format is KRS[position][codon], examples: KRS1TTT, KRS20ATG
# codon is last 3 characters, position is the number between KRS and codon
mutlib <- mutlib |>
  mutate(codon = str_sub(Name, -3, -1),
         position = as.numeric(str_sub(Name, 4, -4)))
```

    ## Warning: There was 1 warning in `mutate()`.
    ## ℹ In argument: `position = as.numeric(str_sub(Name, 4, -4))`.
    ## Caused by warning:
    ## ! NAs introduced by coercion

``` r
# replace NA in position with 0 (KRS_1[codon])
mutlib$position[is.na(mutlib$position)] <- 0

# Translate codons to amino acids
mutlib <- mutlib |>
  left_join(codon_table, by = "codon") |>
  relocate(amino_acid, .after = codon)


# Validate
# Codon should be 3 characters long and only contain A, T, G, or C
stopifnot(str_length(mutlib$codon) == 3)
stopifnot(all(str_detect(mutlib$codon, "[ATGC]")))

# Position should be a number between 0 and 21
stopifnot(all(mutlib$position >= 0 & mutlib$position <= 21))

# Expect 64 codons at each position, once each per dataset
stopifnot(ncol(table(mutlib$position, mutlib$codon)) == 64)
stopifnot(all(table(mutlib$position, mutlib$codon) == length(paths)))

# Mark the wildtype codons
# Here's the sequence from 0 to 21:
# AAA   GTA TTT AGA AAT GAA GGT ATA GAT AAT ACA CAT AAT CCT GAA TTT ACT TCG TGT GAA TTT TAT
wildtype_codon_sequence <- c("AAA", "GTA", "TTT", "AGA", "AAT", "GAA", "GGT", "ATA", "GAT", "AAT", "ACA", "CAT", "AAT", "CCT", "GAA", "TTT", "ACT", "TCG", "TGT", "GAA", "TTT", "TAT")

mutlib <- mutlib |>
  mutate(wildtype_codon = wildtype_codon_sequence[position + 1],
         is_wildtype = codon == wildtype_codon)

# Add a column for the wildtype amino acid
wildtype_residue_sequence <- mutlib |>
  filter(dataset == 1, is_wildtype) |>
  arrange(position) |>
  pull(amino_acid)

mutlib <- mutlib |>
  mutate(wildtype_residue = wildtype_residue_sequence[position + 1],
         is_wildtype_residue = amino_acid == wildtype_residue)

# Validate wildtype
stopifnot(all(table(mutlib$is_wildtype, mutlib$is_wildtype_residue)[2, ] == c(0, 66)))

# Label conditions
mutlib$condition <- factor(mutlib$dataset, levels = 1:3,
                           labels = c("Control", "Cladosporin", "706"))

# Add a single column that labels wildtype, synonymous, non-synonymous or stop
mutlib <- mutlib |>
  mutate(mutation_type = case_when(is_wildtype ~ "Wildtype",
                                   is_wildtype_residue ~ "Synonymous",
                                   amino_acid == "Stop" ~ "Stop",
                                   TRUE ~ "Mutant"))

# Add residue labels
offset = 327
mutlib <- mutlib |>
  mutate(residue = paste0(wildtype_residue, position + offset),
         mutation = paste0(residue, amino_acid))

# Write to file with hash for reproducibility
hash <- digest::digest(mutlib)
write_delim(mutlib, paste0("../data/mutlib-", hash, ".tsv"), delim = "\t")

mutlib
```

    ## # A tidytable: 4,224 × 15
    ##    dataset Name    Sequence       Count codon amino_acid position wildtype_codon
    ##      <int> <chr>   <chr>          <dbl> <chr> <chr>         <dbl> <chr>         
    ##  1       1 KRS1TTT TATGAAATTGGTA…  8945 TTT   Phe               1 GTA           
    ##  2       1 KRS1TTC TATGAAATTGGTA…  2152 TTC   Phe               1 GTA           
    ##  3       1 KRS1TTA TATGAAATTGGTA…  1682 TTA   Leu               1 GTA           
    ##  4       1 KRS1TTG TATGAAATTGGTA…  2826 TTG   Leu               1 GTA           
    ##  5       1 KRS1CTT TATGAAATTGGTA…  2245 CTT   Leu               1 GTA           
    ##  6       1 KRS1CTC TATGAAATTGGTA…   863 CTC   Leu               1 GTA           
    ##  7       1 KRS1CTA TATGAAATTGGTA…   850 CTA   Leu               1 GTA           
    ##  8       1 KRS1CTG TATGAAATTGGTA…  1099 CTG   Leu               1 GTA           
    ##  9       1 KRS1ATT TATGAAATTGGTA…  2254 ATT   Ile               1 GTA           
    ## 10       1 KRS1ATC TATGAAATTGGTA…  1378 ATC   Ile               1 GTA           
    ## # ℹ 4,214 more rows
    ## # ℹ 7 more variables: is_wildtype <lgl>, wildtype_residue <chr>,
    ## #   is_wildtype_residue <lgl>, condition <fct>, mutation_type <chr>,
    ## #   residue <chr>, mutation <chr>

``` r
# Report where we save the file, what the hash was and other file stats
cat("Saved to ../data/mutlib-", hash, ".tsv\n", sep = "")
```

    ## Saved to ../data/mutlib-1187736dfa7e94e0c60d8ef17aa56979.tsv

``` r
cat("Number of rows: ", nrow(mutlib), "\n")
```

    ## Number of rows:  4224

``` r
cat("Number of columns: ", ncol(mutlib), "\n")
```

    ## Number of columns:  15
