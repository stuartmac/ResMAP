02_mutlib-validation
================
Stuart MacGowan
2024-02-10

# Introduction

The mutlib data is a Plasmodium KRS deep mutational scanning experiment
designed to identify potential resistance mutants against wildtype
inhibitors. The data contains 3 experiments, 1 control and 2 treatment
conditions.

# Results

## Comparing pre- and post-selection libraries

``` r
# Control Median
mutlib |>
  filter(condition == "Control") |>
  group_by(condition) |> mutate(CPM = (Count / sum(Count) * 10^6)) |>
  pull(CPM) |>
  median() ->
  control_median

mutlib |>
  # Rename treatment 1_Control: Pre-selection, 2_Cladosporin: Low dose, 1_Cladosporin: High dose
  mutate(treatment = factor(condition, levels = c("Control", "706", "Cladosporin"),
                            labels = c("Pre-selection", "DDD01510706", "Cladosporin"))) |>
  filter(position != 0 & position != 21) |>
  filter(mutation_type != "Wildtype") |>
  # Counts per million
  group_by(condition) |> mutate(CPM = (Count / sum(Count) * 10^6)) |>
  mutate(mutation_type = factor(mutation_type, levels = c("Stop", "Mutant", "Synonymous"))) |>
  
  # reduce mutation type to Synonymous and Non-synonymous
  mutate(mutation_type = ifelse(mutation_type %in% c("Stop", "Mutant"), "Non-synonymous", "Synonymous")) |>
  ggplot(aes(x = CPM, fill = mutation_type, colour = mutation_type)) +
  geom_density(colour = "black", alpha = 0.5, bw = 0.4) +
  geom_rug(alpha = 0.5) +
  geom_vline(xintercept = control_median, linetype = 2,
             colour = "black") +
  # facet_grid(condition ~ batch, scales = "free_y", labeller = ) +
  facet_grid(treatment ~ ., scales = "free_y", labeller = ) +
  # scale_x_log10() +
  scale_x_log10(breaks = c(1, 10, 100, 1000, 10000, 100000),
                labels = c("", "10", "", "1,000", "", "100,000"),
                minor_breaks = NULL, ) +
  scale_fill_manual(values = mutation_type_colours) +
  scale_colour_manual(values = mutation_type_colours) +
  xlab("Counts per million") +
  ylab("Desnity") +
  guides(fill = FALSE, colour = FALSE) +
  theme_bw() +
  theme(aspect.ratio = 0.33, axis.text.x = element_text(angle = 30)) ->
  density_plot
```

    ## Warning: The `<scale>` argument of `guides()` cannot be `FALSE`. Use "none" instead as
    ## of ggplot2 3.3.4.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
density_plot
```

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning: Removed 412 rows containing non-finite values (`stat_density()`).

![](02_mutlib-validation_files/figure-gfm/pre%20vs.%20post%20selection%20density%20plots-1.png)<!-- -->

``` r
# Save SVG
ggsave("figures/02-density-plot.svg", plot = density_plot,
       width = 7.5, height = 7.5, units = "cm")
```

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning: Removed 412 rows containing non-finite values (`stat_density()`).

## Basic statistics

``` r
mutlib |>
  # de-duplicate sequences because wildtype is duplicated
  distinct(dataset, Sequence, .keep_all = TRUE) |>
  group_by(condition) |>
  summarize(total_count = sum(Count))
```

    ## # A tibble: 3 × 2
    ##   condition   total_count
    ##   <chr>             <dbl>
    ## 1 706             3120265
    ## 2 Cladosporin     2413583
    ## 3 Control         2209362

``` r
mutlib |>
  filter(!is_wildtype) |>  # exclude wildtype sequence
  group_by(condition, position) |>
  summarize(min = min(Count),
            max = max(Count),
            mean = mean(Count),
            median = median(Count),
            sd = sd(Count),
            n = n()) |>
  arrange(condition, position)
```

    ## # A tibble: 66 × 8
    ## # Groups:   condition [3]
    ##    condition position   min    max   mean median      sd     n
    ##    <chr>        <dbl> <dbl>  <dbl>  <dbl>  <dbl>   <dbl> <int>
    ##  1 706              0     0    354   15.6      0    56.1    63
    ##  2 706              1     3 186445 8716.     200 29760.     63
    ##  3 706              2     0  18925  950.      25  3456.     63
    ##  4 706              3     0    209   10.7      1    32.9    63
    ##  5 706              4     0    276   15.3      1    48.6    63
    ##  6 706              5     0    278   21.0      7    50.1    63
    ##  7 706              6     0   2149   52.4      8   270.     63
    ##  8 706              7     0   1427   65        6   215.     63
    ##  9 706              8     0    225   14.7      3    37.8    63
    ## 10 706              9     0   7588  247.       8  1048.     63
    ## # ℹ 56 more rows

``` r
mutlib |>
  # Shift position by 0.2 to distinguish the wildtype residues
  mutate(position = position + -0.4 * (-0.5 * is_wildtype_residue)) |>
  ggplot(aes(x = position, y = Count, colour = is_wildtype_residue)) +
  facet_grid(condition ~ .) +
  geom_point() +
  scale_y_log10() +
  # scale_x_continuous(breaks = 0:21, labels = unique(mutlib$residue))
  guides(colour = FALSE) +
  theme_bw()
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](02_mutlib-validation_files/figure-gfm/Manhatten%20style%20plot-1.png)<!-- -->

## Potential codon bias

``` r
codon_rank_control <- mutlib |>
  filter(position != 0 & position != 21) |>
  filter(mutation_type != "Wildtype") |>
  filter(condition == "Control") |>
  group_by(codon) |>
  summarise(mean_count = mean(Count)) |>
  arrange(desc(mean_count))

mutlib |>
  filter(position != 0 & position != 21) |>
  filter(mutation_type != "Wildtype") |>
  filter(condition == "Control") |>
  mutate(codon = factor(codon, levels = codon_rank_control$codon)) |>
  mutate(gc_content = stringr::str_count(codon, "G|C") / 3) |>
  ggplot(aes(x = Count, y = codon, fill = gc_content)) +
  geom_density_ridges(scale = 2, rel_min_height = 0.001) +
  scale_y_discrete(expand = c(0, 0), guide = guide_axis(n.dodge=4)) +     # will generally have to set the `expand` option
  scale_x_log10(expand = c(0, 0)) +
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  guides(fill = FALSE) +
  theme_ridges()
```

    ## Picking joint bandwidth of 0.147

![](02_mutlib-validation_files/figure-gfm/codon%20bias%20ridge%20plot-1.png)<!-- -->

``` r
ggsave('figures/02-codon-ridge-plot.svg', width = 12, height = 20, units = 'cm')
```

    ## Picking joint bandwidth of 0.147
