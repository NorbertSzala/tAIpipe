# Static plotting configuration guide

## Audit of the supplied scripts

The supplied plotting scripts are internally consistent, but the consistency is achieved through repeated hard-coded values rather than reusable settings.

### Repeated patterns

- `theme_stata()` is used in most older scripts; newer scripts use `theme_minimal()` or `theme_bw()`.
- Most plots are exported as `8 x 6 in, 300 dpi`; large faceted plots use `12 x 18 in, 600 dpi`.
- Binary categories are usually represented by dark and light grey.
- Multi-category plots frequently use `colorRampPalette(brewer.pal(12, "Paired"))(n)`.
- Violin plots are commonly overlaid with notched boxplots.
- Several final figures are combined with `cowplot`.

## What belongs in the config

| Put in `plotting.yaml`                               | Keep in plotting/analysis scripts                      |
| ---------------------------------------------------- | ------------------------------------------------------ |
| Font family and font sizes                           | Which variables are mapped to x, y, colour or facets   |
| Standard output sizes, formats and DPI               | Statistical transformations and tests                  |
| Theme, grid, margins and legend defaults             | Correlation method, smoothing method and model formula |
| Point, line, boxplot and annotation defaults         | Histogram bin width when it has analytical meaning     |
| Reusable palettes and exact category-colour mappings | Axis limits derived from data or scientific thresholds |
| Standard layout profiles (`paired`, `2x2`, stacked)  | Category order unless it is globally meaningful        |
| Reusable shape and linetype sequences                | Titles, subtitles and biological interpretation        |
| Standard plot presets and export directories         | Outlier handling and filtering                         |


## Colour strategy

### Palette selection by variable type

| Data                            | Recommended default        | Notes                                                               |
| ------------------------------- | -------------------------- | ------------------------------------------------------------------- |
| Binary state                    | Grey + blue                | Keep absence/reference neutral and highlight presence               |
| 3-8 unordered groups            | Okabe-Ito                  | Add shape for points and linetype for lines                         |
| 9-10 unordered groups           | Petroff 10                 | Prefer facets or grouped panels; never interpolate                  |
| Ordered categories              | Sequential palette         | Lightness should change monotonically                               |
| Continuous values               | Viridis/cividis            | Perceptually ordered and suitable for grayscale                     |
| Values around a meaningful zero | Diverging blue-neutral-red | Fix the midpoint and symmetric limits when comparisons require them |
| Missing values                  | Neutral grey               | Missingness must not resemble an extreme value                      |



## Recommended changes to existing figures

### Binary feature versus tAI

Current approach: violin + notched boxplot.

Preferred default: boxplot with jittered raw observations, or a half-violin/raincloud plot when distribution shape is important. Remove notches unless median-confidence interpretation is needed and group sizes are adequate.

### Histograms comparing groups

Overlapping or dodged histograms become difficult to compare. Prefer:

- faceted histograms with the same bin boundaries;
- density curves for smooth comparison;
- ECDF curves for distribution shifts;
- ridgeline plots for several groups.

### LCR/protein length versus tAI

For large gene-level datasets, replace dense scatterplots with `geom_hex()` or 2D binning. Overlay a transparent trend and report the correlation/model separately. LOESS is useful for exploration but should not be a default inferential model.

### Ordered count or length bins

Do not assign a distinct hue to every ordered bin. Use one sequential scale, or display median/mean with confidence intervals as a line across ordered bins. This makes the order explicit.

### Correlation heatmap

The lower-triangle layout is appropriate. Keep a fixed `[-1, 1]` scale and a neutral zero midpoint. Consider:

- abbreviating variable labels;
- separating effect size from significance using text weight, a dot or an outline;
- clustering variables only for exploratory figures, not when a predefined biological order is required;
- exporting as PDF when many labels are present.

### tRNA codon profiles

The current vertically stacked 12 x 18 inch faceted barplots are hard to scan. Better alternatives:

- heatmap: genomes as rows, codons as columns, normalized abundance as fill;
- small multiples by amino acid, with codons within each panel;
- dot plot: position = codon/genome, size = count, colour = normalized abundance;
- a focused profile plot for selected genomes plus a separate summary heatmap for all genomes.

### Top/bottom tAI comparisons

Use dumbbell or forest plots for paired top-versus-bottom summaries. Plot effect size and confidence interval rather than only separate group distributions. This directly shows direction and magnitude.

### GO enrichment

Use a dot/lollipop plot with:

- x = enrichment ratio or gene ratio;
- y = reordered GO term;
- point size = count;
- colour = adjusted p-value on a sequential scale;
- facets = ontology or comparison.

Limit labels to the most relevant terms per panel and wrap long descriptions.

### Genome-level tAI by phylum/lifestyle

Use boxplot + jitter or raincloud plots with nested facets. When phylum has many categories, facet by phylum and colour by the smaller lifestyle variable rather than colouring every phylum.

## Layout policy

Use `patchwork` and the presets in the YAML file:

- `paired`: two related plots sharing one legend;
- `stacked`: vertically aligned plots with shared x semantics;
- `grid_2x2`: four comparable panels;
- `grid_3col`: compact small multiples.

Do not centralize data-dependent facet counts completely. Use `recommended_facet_size()` to derive a minimum canvas from the number of panels, then override it when labels are unusually long.

## Minimal usage

```r
library(ggplot2)
source("R/plotting_utils.R")
cfg <- read_plot_config("config/plotting.yaml")

levels_group <- levels(df$group)

p <- ggplot(df, aes(x, y, colour = group, shape = group)) +
  geom_point(
    size = cfg$geoms$point$size,
    alpha = cfg$geoms$point$alpha
  ) +
  scale_colour_project(
    cfg,
    palette = "petroff_10",
    levels = levels_group
  ) +
  scale_shape_project(cfg, levels = levels_group) +
  project_theme(cfg)

save_plot(p, "my_plot", cfg, size = "double_column")
```

For stable binary labels:

```r
p +
  scale_fill_project(
    cfg,
    mapping = "signal_presence",
    levels = c("Signal absence", "Signal presence")
  )
```

For a continuous scale:

```r
p + scale_continuous_project(cfg, aesthetic = "fill", palette = "sequential")
```

For a correlation scale:

```r
p + scale_continuous_project(
  cfg,
  aesthetic = "fill",
  palette = "diverging",
  limits = c(-1, 1),
  midpoint = 0
)
```