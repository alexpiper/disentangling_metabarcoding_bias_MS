#Run with SAVE_FIGURES = TRUE to save figures in figures/.

SAVE_FIGURES = TRUE
#Libraries and paths
library(here)
library(tidyverse)
library(ggthemes)
library(cowplot) 
theme_set(theme_grey()) 
library(ggbeeswarm)
# Functions for bias estimation and calibration
library(metacal)
packageVersion("metacal")

# This package
devtools::load_all(here())

#Options for ggplot:
  
  ase_theme <- theme_tufte() + 
  theme(
    text = element_text(size=9, family = ""),
    strip.text = element_text(size=9, family = ""),
    # axis.text = element_text(size=8, family = ""),
    legend.position = "none"
  )
tax_theme <- theme(
  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
  axis.title.x = element_blank())

update_geom_defaults("point", list(size = 1))
update_geom_defaults("text", list(size = 2.5))
update_geom_defaults("hline", list(color = "grey"))
update_geom_defaults("abline", list(color = "grey"))
# Colors for color scales
colors.OkabeIto <-  tribble(
  ~Name, ~Hex,
  "white", "#FFFFFF",
  "silver", "#C0C0C0",
  "gray", "#808080",
  "black", "#000000",
  "red", "#FF0000",
  "maroon", "#800000",
  "yellow", "#FFFF00",
  "olive", "#808000",
  "lime", "#00FF00",
  "green", "#008000",
  "aqua", "#00FFFF",
  "teal", "#008080",
  "blue", "#0000FF",
  "navy", "#000080",
  "fuchsia", "#FF00FF",
  "purple",  "#800080")

colors.taxon <- c(
  "Acizzia_alternata" = "#FFFFFF",
  "Acizzia_solanicola" = "#C0C0C0",
  "Aphidius_colemani" = "#808080",
  "Bactrocera_tryoni" = "#000000",
  "Bradysia_ocellaris" = "#FF0000",
  "Carpophilus_davidsoni" = "#800000",
  "Carpophilus_dimidiatus" = "#FFFF00",
  "Diuraphis_noxia" = "#808000",
  "Drosophila_hydeii" = "#00FF00"  ,
  "Drosophila_melanogaster" = "#008000",
  "Drosophila_simulans" = "#00FFFF",
  "Lysiphlebus_testaceipes" = "#008080",
  "Metopolophium_dirhodum" = "#0000FF",
  "Psyllid_sp" = "#000080",
  "Rhopalosiphum_padi" = "#FF00FF",
  "Scaptodrosophila_lativittata" = "#800080"
) 
colors.num_species <- c(
  "14" = "#E69F00",
  "15" = "#56B4E9",
  "16" = "#CC79A7"
)


###START HERE

#The vast majority of reads were classified above threshold,

main0 %>%
  group_by(Table) %>%
  summarize(Sum = sum(Count)) %>%
  mutate(Proportion = Sum / sum(Sum))

#In preparing this dataframe, we have grouped all reads classified to non-mock taxa as “Other”. The vast majority of above and below threshold reads were classified to the mock taxa,

main0 %>%
  group_by(Table, Taxon != "Other") %>%
  summarize(Sum = sum(Count)) %>%
  mutate(Proportion = Sum / sum(Sum))

#We will follow Brooks2015 in using the above threshold reads only, and will restrict to just the mock taxa. We will also create a new dataframe that will be our main data frame for our analysis going forward.

main <- main0 %>%
  filter(Table == "above", Taxon != "Other") %>%
  select( ï..Sample, Taxon, Count)


####Prep for main analysis
####Our main analysis will use Observed and Actual variables, with the following simplifications. Our downstream analysis requires that the samples only have reads from taxa actually in the samples. We will therefore remove the out-of-sample reads and treat the actual compositions as an even mixture of the expected taxa. We will also ignore the extra precision information in the observed read counts, and treat Observed as well as Actual as compositional vectors. We will further filter to samples with more than one species. We will generally store relative abundances as proportions for plotting purposes (we use the close_elts function within samples to convert compositional vectors to vectors of proportions).

main0 <- main %>%
  filter(Expected, Num_species > 1) %>%
  mutate_by(ï..Sample, 
            Observed = close_elts(Count),
            Actual = close_elts(Expected),
            Error = Observed / Actual) %>%
  select(-Count, -Expected) %>%
  select(ï..Sample, Taxon, Observed, Actual, everything())

main0 %>%
  mutate_at(vars(Actual, Observed), logit) %>%
  ggplot(aes(Actual, Observed, color = Taxon)) + 
  geom_quasirandom() +
  facet_grid(Mixture_type ~ .) +
  geom_rangeframe(sides = "bl", color= "black") +
  scale_color_manual(values = colors.taxon) + 
  labs(title = "Observed vs. Actual (log-odds)")

main0 %>%
  mutate_at(vars(Actual, Observed), logit) %>%
  ggplot(aes(Actual, Observed, color = Taxon)) + 
  geom_quasirandom() +
  facet_grid(Mixture_type ~ Taxon) +
  geom_rangeframe(sides = "bl", color= "black") +
  labs(title = "Observed vs. Actual (log-odds)") +
  scale_color_manual(values = colors.taxon) + 
  theme(legend.position = "none")

gvs <- c("Mixture_type", "ï..Sample", "Num_species")
ratios <- main0 %>%
  select(-Plate, -Barcode) %>%
  compute_ratios(group_vars = gvs) %>%
  filter(!is.na(Actual)) %>%
  mutate(Pair = paste(Taxon.x, Taxon.y, sep = ":"))

cmb <- combn(taxa, 2, simplify = FALSE) %>% 
  map_chr( ~paste(.[1], .[2], sep = ":"))


p.pw <- ratios %>% 
  filter(Pair %in% cmb) %>%
  ggplot(aes(Pair, Error, color = as.factor(Num_species))) +
  geom_hline(yintercept = 1) +
  facet_grid(Mixture_type ~ .) +
  geom_rangeframe(sides = "l", color= "black") + 
  scale_y_log10() +
  scale_x_discrete() + 
  scale_color_manual(values = colors.num_species) +
  tax_theme +
  theme(legend.position = "bottom")


set.seed(20190517)
temp <- main0 %>%
  select(Mixture_type, ï..Sample, Taxon, Observed, Actual) %>%
  mutate(Error = Observed / Actual) %>%
  group_by(Mixture_type) %>%
  nest %>%
  mutate(
    Error_matrix = map(data, build_matrix, 
                       ï..Sample, Taxon, Error, fill = NaN),
    Estimate = map(Error_matrix, center, enframe = TRUE),
    Bootreps = map(Error_matrix, bootrep_center, R = 4e3)
  )
est <- temp %>%
  unnest(Estimate) %>%
  rename(Bias = Center)
bootreps <- temp %>%
  unnest(Bootreps)
stats <- bootreps %>%
  group_by(Mixture_type, Taxon) %>%
  summarize(gm_mean = gm_mean(Center), gm_se = gm_sd(Center)) %>%
  ungroup
bias <- left_join(est, stats, by = c("Mixture_type", "Taxon"))
rm(temp)




The Observed and Actual compositions are normalized to proportions for plotting in the next section.

Add copy-number bias for each taxon + experiment pair:
  
  main0 <- main0 %>%
  left_join(cnbias %>% select(Mixture_type, Taxon, CN_bias),
            by = c("Mixture_type", "Taxon"))
print(main0, n = 10)
## # A tibble: 591 x 11
##    Sample Taxon Observed Actual Plate Barcode Mixture_type Num_species
##    <chr>  <chr>    <dbl>  <dbl> <dbl>   <dbl> <chr>              <int>
##  1 s1-1   Lact…   0.602   0.333     1       1 Cells                  3
##  2 s1-1   Prev…   0.332   0.333     1       1 Cells                  3
##  3 s1-1   Stre…   0.0663  0.333     1       1 Cells                  3
##  4 s1-10  Atop…   0.0391  0.333     1      10 Cells                  3
##  5 s1-10  Lact…   0.392   0.333     1      10 Cells                  3
##  6 s1-10  Snea…   0.569   0.333     1      10 Cells                  3
##  7 s1-11  Gard…   0.0154  0.333     1      11 Cells                  3
##  8 s1-11  Lact…   0.337   0.333     1      11 Cells                  3
##  9 s1-11  Lact…   0.647   0.333     1      11 Cells                  3
## 10 s1-12  Atop…   0.0616  0.5       1      12 Cells                  2
## # … with 581 more rows, and 3 more variables: Species_list <chr>,
## #   Error <dbl>, CN_bias <dbl>
View the errors
Let’s take a look at the observed vs. actual proportions, taking a log-odds transform to avoid compressing the variation near p=0 and p=1.

main0 %>%
  mutate_at(vars(Actual, Observed), logit) %>%
  ggplot(aes(Actual, Observed, color = Taxon)) + 
  geom_quasirandom() +
  facet_grid(Mixture_type ~ .) +
  geom_rangeframe(sides = "bl", color= "black") +
  scale_color_manual(values = colors.taxon, labels = tax_labeller) + 
  labs(title = "Observed vs. Actual (log-odds)")


There is significant error. Moreover, it is not consistent even for a particular species and actual abundance (except for the 7-species mixtures), as our model predicts:
  
  main0 %>%
  mutate_at(vars(Actual, Observed), logit) %>%
  ggplot(aes(Actual, Observed, color = Taxon)) + 
  geom_quasirandom() +
  facet_grid(Mixture_type ~ Taxon) +
  geom_rangeframe(sides = "bl", color= "black") +
  labs(title = "Observed vs. Actual (log-odds)") +
  scale_color_manual(values = colors.taxon, labels = tax_labeller) + 
  theme(legend.position = "none")


What about after copy-number correction?
  
  main0 <- main0 %>%
  mutate_by(Sample, CN_corrected = close_elts(Observed / CN_bias))
main0 %>%
  mutate_at(vars(Actual, CN_corrected), logit) %>%
  ggplot(aes(Actual, CN_corrected, color = Taxon)) + 
  geom_quasirandom() +
  facet_grid(Mixture_type ~ .) +
  geom_rangeframe(sides = "bl", color= "black") +
  scale_color_manual(values = colors.taxon, labels = tax_labeller) + 
  labs(title = "Observed with CN correction vs. Actual (log-odds)")


The error in the DNA mixtures looks to be reduced. We return to the efficacy of CN correction later on.

Now let’s look at the taxon ratios, which we predict should show a consistent error. To do so, we first need a data frame with all the pairwise ratios for each sample,

gvs <- c("Mixture_type", "Sample", "Num_species")
ratios <- main0 %>%
  select(-Plate, -Barcode) %>%
  compute_ratios(group_vars = gvs) %>%
  filter(!is.na(Actual)) %>%
  mutate(Pair = paste(Taxon.x, Taxon.y, sep = ":"))
Pick a non-redundant set of pairs for plotting,

cmb <- combn(taxa, 2, simplify = FALSE) %>% 
  map_chr( ~paste(.[1], .[2], sep = ":"))
and a reduced table for plotting the copy-number prediction (which only depends on the taxa pair and not the sample)

cn_ratios <- ratios %>%
  select(Mixture_type, Pair, CN_bias) %>%
  distinct %>%
  filter(Pair %in% cmb)
Plot the error ratios in along with the error predicted by 16S CN (dark red cross),

p.pw <- ratios %>% 
  filter(Pair %in% cmb) %>%
  ggplot(aes(Pair, Error, color = as.factor(Num_species))) +
  geom_hline(yintercept = 1) +
  facet_grid(Mixture_type ~ .) +
  geom_rangeframe(sides = "l", color= "black") + 
  scale_y_log10() +
  scale_x_discrete(labels = tax_labeller) + 
  scale_color_manual(values = colors.num_species) +
  tax_theme +
  theme(legend.position = "bottom")
p.pw +
  geom_point(data = cn_ratios, aes(Pair, CN_bias), 
             color = "darkred", shape = 3, size = 3) +
  geom_quasirandom()


The error in ratios is consistent across samples, as our model predicts.

Bias estimation
Estimate bias w/ standard error,

set.seed(20190517)
temp <- main0 %>%
  select(Mixture_type, Sample, Taxon, Observed, Actual) %>%
  mutate(Error = Observed / Actual) %>%
  group_by(Mixture_type) %>%
  nest %>%
  mutate(
    Error_matrix = map(data, build_matrix, 
                       Sample, Taxon, Error, fill = NaN),
    Estimate = map(Error_matrix, center, enframe = TRUE),
    Bootreps = map(Error_matrix, bootrep_center, R = 4e3)
  )
est <- temp %>%
  unnest(Estimate) %>%
  rename(Bias = Center)
bootreps <- temp %>%
  unnest(Bootreps)
stats <- bootreps %>%
  group_by(Mixture_type, Taxon) %>%
  summarize(gm_mean = gm_mean(Center), gm_se = gm_sd(Center)) %>%
  ungroup
bias <- left_join(est, stats, by = c("Mixture_type", "Taxon"))
rm(temp)
Check the results.

bias %>%
  print(n = Inf)
## # A tibble: 21 x 5
##    Mixture_type Taxon                     Bias gm_mean gm_se
##    <chr>        <chr>                    <dbl>   <dbl> <dbl>
##  1 Cells        Atopobium_vaginae        0.285   0.285  1.04
##  2 Cells        Gardnerella_vaginalis    0.160   0.160  1.05
##  3 Cells        Lactobacillus_crispatus  2.29    2.28   1.03
##  4 Cells        Lactobacillus_iners      4.68    4.68   1.02
##  5 Cells        Prevotella_bivia         1.79    1.79   1.04
##  6 Cells        Sneathia_amnii           4.59    4.59   1.04
##  7 Cells        Streptococcus_agalactiae 0.250   0.250  1.03
##  8 DNA          Atopobium_vaginae        1.06    1.06   1.03
##  9 DNA          Gardnerella_vaginalis    0.403   0.403  1.05
## 10 DNA          Lactobacillus_crispatus  0.536   0.536  1.04
## 11 DNA          Lactobacillus_iners      2.31    2.31   1.02
## 12 DNA          Prevotella_bivia         0.392   0.393  1.03
## 13 DNA          Sneathia_amnii           2.39    2.39   1.04
## 14 DNA          Streptococcus_agalactiae 2.01    2.01   1.02
## 15 PCR_product  Atopobium_vaginae        1.03    1.03   1.04
## 16 PCR_product  Gardnerella_vaginalis    0.817   0.818  1.05
## 17 PCR_product  Lactobacillus_crispatus  0.907   0.907  1.03
## 18 PCR_product  Lactobacillus_iners      1.19    1.19   1.05
## 19 PCR_product  Prevotella_bivia         0.927   0.926  1.05
## 20 PCR_product  Sneathia_amnii           1.30    1.30   1.03
## 21 PCR_product  Streptococcus_agalactiae 0.911   0.911  1.04
bias <- select(bias, -gm_mean)
Add the bias estimate, predicted proportions, and residual compositional error

main0 <- main0 %>%
  left_join(bias, by = c("Mixture_type", "Taxon")) %>%
  mutate_by(Sample, 
            Predicted = close_elts(Actual * Bias),
            Residual = Observed / Predicted
  )
And also get the predicted ratios,

ratios.pred <- bias %>%
  select(Mixture_type, Taxon, Bias) %>%
  compute_ratios(group_vars = c("Mixture_type")) %>%
  rename(Predicted = Bias) %>%
  mutate(Pair = paste(Taxon.x, Taxon.y, sep = ":"))
Ratio errors are explained by the bias (Supplemental Figure SBrooksRatios)

# Put predictions underneath the geom layer showing the observed errors
p.pw +
  geom_point(data = ratios.pred %>% filter(Pair %in% cmb),
             aes(Pair, Predicted), inherit.aes = FALSE,
             shape = 3, size = 4, color = "black") +
  geom_quasirandom() +
  labs(color = "# species\nin mixture", y = "Observed / Actual") +
  base_theme +
  tax_theme +
  theme(legend.position = "right",
        plot.margin = margin(l = 0.2, unit = "in")
  )


Get the prediction of our measurements given the estimated bias. For comparison, also get the prediction under no bias and under copy-number bias,

pred <- main0 %>%
  mutate(`No bias` = 1) %>%
  rename(`Copy-number bias` = CN_bias, `Estimated bias` = Bias) %>%
  gather("Bias_type", "Bias", 
         `No bias`, `Copy-number bias`, `Estimated bias`) %>%
  mutate_by(c(Sample, Bias_type), Predicted = close_elts(Actual * Bias)) %>%
  mutate(Bias_type = factor(Bias_type, 
                            c("No bias", "Copy-number bias", "Estimated bias")))
Check that the predictions accounting for bias reduce the error by various proportion-based (non-compositional) measures.

error <- pred %>%
  group_by(Mixture_type, Bias_type) %>%
  summarize(
    MSE.prop = mean((Observed - Predicted)^2),
    MSE.logit = mean((logit(Observed) - logit(Predicted))^2),
    RMSE.logit = sqrt(mean( (logit(Observed) - logit(Predicted))^2 )),
  )
error
## # A tibble: 9 x 5
## # Groups:   Mixture_type [3]
##   Mixture_type Bias_type        MSE.prop MSE.logit RMSE.logit
##   <chr>        <fct>               <dbl>     <dbl>      <dbl>
## 1 Cells        No bias           0.0853     3.49        1.87 
## 2 Cells        Copy-number bias  0.0859     3.25        1.80 
## 3 Cells        Estimated bias    0.00106    0.0506      0.225
## 4 DNA          No bias           0.0419     1.11        1.06 
## 5 DNA          Copy-number bias  0.0198     0.546       0.739
## 6 DNA          Estimated bias    0.00136    0.0423      0.206
## 7 PCR_product  No bias           0.00585    0.120       0.346
## 8 PCR_product  Copy-number bias  0.00585    0.120       0.346
## 9 PCR_product  Estimated bias    0.00323    0.0696      0.264
Summary for the manuscript: The proportions predicted from the fitted model reduced the mean squared error by 98.8%.

Also check the reduction in compositional error using the Aitchison distance,

pred %>%
  group_by(Mixture_type, Bias_type, Sample) %>%
  summarize(Adist = anorm(Observed / Predicted)) %>%
  summarize(Adist = mean(Adist))
## # A tibble: 9 x 3
## # Groups:   Mixture_type [3]
##   Mixture_type Bias_type        Adist
##   <chr>        <fct>            <dbl>
## 1 Cells        No bias          1.77 
## 2 Cells        Copy-number bias 1.70 
## 3 Cells        Estimated bias   0.178
## 4 DNA          No bias          0.995
## 5 DNA          Copy-number bias 0.702
## 6 DNA          Estimated bias   0.179
## 7 PCR_product  No bias          0.297
## 8 PCR_product  Copy-number bias 0.297
## 9 PCR_product  Estimated bias   0.235
Plot the log-odds with mean-squared errors of the log-odds proportion,

p <- ggplot(pred, aes(logit(Predicted), logit(Observed), color = Taxon)) +
  geom_abline(intercept = 0, slope = 1, color = "grey") +
  geom_jitter(width = 0.1, height = 0) +
  geom_rangeframe(color = "black") + 
  facet_grid(Mixture_type ~ Bias_type) +
  base_theme +
  labs(x = "log-odds(Predicted proportion)", 
       y = "log-odds(Observed proportion)") +
  coord_fixed() +
  scale_color_manual(values = colors.taxon, labels = tax_labeller) + 
  theme(
    panel.spacing.x = unit(1, "lines"),
    legend.position = "bottom",
  ) + 
  geom_text(data = error, 
            aes(x = 2.5, y = -4, 
                label = paste("MSE:", round(MSE.logit, 2))),
            color = "black")
p


Main text 4-panel figure
These are the taxa used in the second row of the main figure.

plot_taxa <- c("Gardnerella_vaginalis", "Lactobacillus_crispatus",
               "Lactobacillus_iners", "Sneathia_amnii", "Streptococcus_agalactiae")
Observed proportion (odds) vs. expected (Upper left)
Observed proportion vs. model prediction (upper right)
Fold error in proportion (lower left)
Fold error in ratios (lower right)
Top row: Observed vs. expected proportions before and after model fitting

# Tibble for plotting
pred0 <- pred %>%
  filter(Mixture_type == "Cells", Bias_type == "Estimated bias")
# Pick a fixed xy range to use in the top row
range.top <- c(0, 1)
# for y=x ref line
ref_line = tibble(x = range.top[1], y = range.top[1], 
                  xend = range.top[2], yend = range.top[2])
# for (original) formatting the x-axis labels in left panel
xtb <- tibble(labels = c("0", "1/7", "1/4", "1/3", "1/2", "1")) %>%
  rowwise() %>%
  mutate(breaks = eval(parse(text = labels)))
# For setting the rangeframe axes to run from 0 to 1 instead of giving the data
# range
rftb <- tibble(Actual = c(0, 1), Observed = c(0, 1), Predicted = c(0, 1))
# Panel A
p.A <- ggplot(pred0, aes(Actual, Observed, color = Taxon)) +
  geom_segment(data=ref_line, aes(x=x, xend=xend, y=y, yend=yend), 
               colour="grey") +
  # geom_abline(slope = 1, intercept = 0, color = "grey") + 
  geom_quasirandom() +
  geom_rangeframe(data = rftb, color = "black") +
  scale_color_manual(values = colors.taxon, labels = tax_labeller) + 
  scale_x_continuous(limits = range.top,
                     breaks = xtb$breaks, labels = xtb$labels) +
  # scale_x_continuous(limits = range.top, 
  #     breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(limits = range.top, 
                     breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  labs(x = "Actual", y = "Observed",
       title = "Observed proportion vs. actual") +
  base_theme
# Panel B
p.B <- ggplot(pred0, aes(Predicted, Observed, 
                         color = fct_reorder(Taxon, -Observed))) +
  geom_segment(data=ref_line, aes(x=x, xend=xend, y=y, yend=yend), 
               colour="grey") +
  geom_point() +
  geom_rangeframe(data = rftb, color = "black") +
  scale_color_manual(values = colors.taxon, labels = tax_labeller) + 
  scale_x_continuous(limits = range.top, 
                     breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(limits = range.top, 
                     breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  labs(x = "Model prediction", y = "Observed", color = "Taxon",
       title = "Observed proportion vs. model prediction") +
  base_theme
Bottom row: Error in taxon proportions and error in taxon ratios. To give both plots the same error (y-axis) scale, we odds transform the proportions.

# Ratios, for just five of the seven taxa to simplify the plot:
plot_pairs <- combn(plot_taxa, 2, simplify = FALSE) %>% 
  map_chr( ~paste(.[1], .[2], sep = ":"))
ratios0 <- ratios %>% 
  filter(Mixture_type == "Cells", Pair %in% plot_pairs)
ratios0.pred <- ratios.pred %>%
  filter(Mixture_type == "Cells", Pair %in% plot_pairs)
# Add shape keys for figure legend
ratios0 <- ratios0 %>%
  mutate(Shape = "Observed")
ratios0.pred <- ratios0.pred %>%
  mutate(Shape = "Predicted")
# Odds transformed proportions
pred0.odds <- pred0 %>%
  mutate_at(vars(Observed, Actual, Predicted), odds)
# Choose y-axis range
pred0.odds %>%
  {range(.$Observed / .$Actual)}
## [1]  0.02874036 32.32727273
range.bot <- c(0.02,33)
# Panel C: Error in proportions (odds transformed)
p.C <- ggplot(pred0.odds, aes(x = Taxon, y = Observed / Actual,
                              color = as.factor(Num_species))) +
  geom_hline(yintercept = 1) +
  geom_quasirandom() +
  geom_rangeframe(color = "black", sides = "l") +
  labs(title = "Fold error in taxon proportions (odds)",
       y = "Observed / Actual", 
       color = "# species\nin mixture") +
  scale_y_log10(limits = range.bot, breaks = c(0.02, 0.1, 1, 10, 30), 
                labels = log_formatter) +
  scale_x_discrete(labels = tax_labeller) + 
  scale_color_manual(values = colors.num_species) +
  base_theme +
  tax_theme
p.D <- ggplot(ratios0, aes(Pair, Error, color = as.factor(Num_species), 
                           shape = Shape)) +
  # Reference line
  geom_hline(yintercept = 1) +
  # Prediction from bias estimate
  geom_point(data = ratios0.pred,
             aes(Pair, Predicted),
             size = 4, color = "black") +
  # Observed
  # TODO: figure out why the points are smaller in just this panel
  geom_quasirandom(size = 1.4) +
  # Other stuff
  geom_rangeframe(sides = "l", color= "black") + 
  # scale_y_log10(breaks = c(0.02, 0.1, 1, 10, 20),
  #     labels = log_formatter) +
  scale_y_log10(limits = range.bot, breaks = c(0.02, 0.1, 1, 10, 30), 
                labels = log_formatter) +
  scale_x_discrete(labels = tax_labeller) + 
  scale_shape_manual(name = "", values = c(16, 3)) +
  scale_color_manual(values = colors.num_species) +
  labs(title = "Fold error in taxon ratios with model prediction",
       y = "Observed / Actual", 
       color = "# species\nin mixture") +
  base_theme +
  tax_theme
Make the plot!
  
  l.A <- get_legend(p.B + theme(legend.position = "right"))
l.D <- get_legend(p.D + theme(legend.position = "right"))

row1 <- plot_grid(p.A, p.B,
                  labels = c("A", "B"), align = "hv") %>%
  plot_grid(l.A, rel_widths = c(1, 0.14))
row2 <- plot_grid(p.C, p.D, 
                  labels = c("C", "D"), align = "hv") %>%
  plot_grid(l.D, rel_widths = c(1, 0.14))
plot_grid(row1, row2, ncol=1, rel_heights = c(1, 1.1))


Fraction of squared Aitchison error due to bias
err <- main0 %>%
  group_by(Mixture_type, Sample) %>%
  summarize(
    Error2 = anorm(Observed / Actual)^2,
    Bias2 = anorm(Bias)^2,
    Resid2 = anorm(Observed / (Actual * Bias))^2
  ) %>%
  summarize_at(vars(-Sample), sum)
err <- err %>% 
  mutate(Frac_Bias = Bias2 / Error2, Frac_Resid = Resid2 / Error2)
err
## # A tibble: 3 x 6
##   Mixture_type Error2  Bias2 Resid2 Frac_Bias Frac_Resid
##   <chr>         <dbl>  <dbl>  <dbl>     <dbl>      <dbl>
## 1 Cells        267.   263.     3.70     0.986     0.0139
## 2 DNA           87.0   83.4    3.65     0.958     0.0419
## 3 PCR_product    9.18   3.47   5.71     0.378     0.622
Bias of each step in the workflow
Compute Bias for each step:
  
  bias_steps <- bias %>% 
  select(Mixture_type, Taxon, Bias) %>%
  spread(Mixture_type, Bias) %>%
  mutate(Extraction = Cells / DNA,
         PCR = DNA / PCR_product,
         Sequencing = PCR_product
  )
bias_steps
## # A tibble: 7 x 7
##   Taxon                 Cells   DNA PCR_product Extraction   PCR Sequencing
##   <chr>                 <dbl> <dbl>       <dbl>      <dbl> <dbl>      <dbl>
## 1 Atopobium_vaginae     0.285 1.06        1.03       0.268 1.03       1.03 
## 2 Gardnerella_vaginalis 0.160 0.403       0.817      0.396 0.493      0.817
## 3 Lactobacillus_crispa… 2.29  0.536       0.907      4.27  0.591      0.907
## 4 Lactobacillus_iners   4.68  2.31        1.19       2.03  1.94       1.19 
## 5 Prevotella_bivia      1.79  0.392       0.927      4.55  0.423      0.927
## 6 Sneathia_amnii        4.59  2.39        1.30       1.92  1.85       1.30 
## 7 Streptococcus_agalac… 0.250 2.01        0.911      0.124 2.20       0.911
Add copy-number correction as a hypothetical final step in the bioinformatics pipeline. CN correction involves dividing by the CN for each taxon, and thus has a bias given by the reciprocal of the CNs.

cn_corr <- species_info %>% 
  transmute(Taxon, CN_correction = center_elts(1 / Copy_number))
bias_steps <- left_join(bias_steps, cn_corr, by = c("Taxon")) %>%
  mutate(Total = Extraction * PCR * Sequencing * CN_correction)
Gather back into tidy form

bias_steps <- bias_steps %>%
  gather("Type", "Bias", -Taxon)
Get standard errors from the bootreps,

bootreps0 <- bootreps %>%
  spread(Mixture_type, Center) %>%
  mutate(
    Extraction = Cells / DNA,
    PCR = DNA / PCR_product,
    Sequencing = PCR_product
  )
ses <- bootreps0 %>%
  gather("Type", "Bias", -.id, -Taxon) %>%
  group_by(Taxon, Type) %>%
  summarize(gm_se = gm_sd(Bias)) %>%
  ungroup
Join Bias w/ the standard errors:
  
  bias_steps <- left_join(bias_steps, ses, by = c("Type", "Taxon"))
Because CN correction is deterministic, the geometric standard error is 1.

bias_steps <- bias_steps %>%
  mutate(gm_se = ifelse(Type == "CN_correction", 1, gm_se))
Create the 2-panel figure for the main text, showing the bias estimates in the left panel and the predicted composition through the workflow in the right panel.

# We will order the taxa by the in the legend by which are observed as most
# abundant in an actual even mixture; this will be the taxa with the largest
# "Total" bias.
lvls.taxon <- bias_steps %>%
  filter(Type == "Total") %>%
  arrange(Bias) %>%
  .$Taxon
# Order for the steps (left panel)
lvls.step <- c("Extraction", "PCR", "Sequencing", "CN_correction")
labs.step <- c("Extraction", "PCR", "Seq. + Inf.", "CN correction")
names(labs.step) <- c("Extraction", "PCR", "Sequencing", "CN_correction")
# Labels for the positions (right panel)
# labs.position <- c("Actual", "Post extraction", "Post PCR", "Reads", 
#             "CN corrected")
labs.position <- c("Cells", "DNA", "PCR product", "Reads", "Reads / CN")
names(labs.position) <- paste0("T", 0:4)
# Tibble of the relative abundances through the workflow (right panel)
ra_steps <- bias_steps %>%
  select(Taxon, Type, Bias) %>%
  spread(Type, Bias) %>%
  transmute(Taxon,
            T0 = 1,
            T1 = T0 * Extraction, 
            T2 = T1 * PCR, 
            T3 = T2 * Sequencing,
            T4 = T3 * CN_correction,
  ) %>%
  mutate(Taxon = fct_reorder(Taxon, T4)) %>%
  gather("Position", "Abundance", -Taxon)
# shared params
# ylimits <- ra_steps$Abundance %>% range
ylimits <- c(0.1, 6)
xexpand <- c(0.1, 0.1)
## Left panel
# Tibble for reference line
ref.left <- tibble(Step = factor(lvls.step, lvls.step), Bias = 1)
# make plot
p.left <- bias_steps %>% 
  filter(Type %in% lvls.step) %>%
  mutate(
    Step = factor(Type, lvls.step),
    Taxon = factor(Taxon, levels(ra_steps$Taxon))
  ) %>%
  ggplot(aes(Step, Bias, color = Taxon)) +
  geom_path(data = ref.left, aes(group = Bias), color = "grey") +
  # geom_quasirandom(width = 0.2) +
  geom_pointrange(aes(ymin = Bias / gm_se^2, ymax = Bias * gm_se^2),
                  fatten = 0.5,
                  position = position_dodge(width = 0.3)) +
  scale_y_log10(limits = ylimits, labels = log_formatter) +
  scale_x_discrete(expand = xexpand,labels = labs.step) +
  geom_rangeframe(color = "black") +
  scale_color_manual(values = colors.taxon, labels = tax_labeller) + 
  labs(x = "Step", y = "Efficiency / gm. mean", 
       title = "Bias of each step") +
  guides(color = guide_legend(reverse = TRUE)) +
  base_theme +
  theme(axis.text.x = element_text(size=8, family = ""))
## Right panel
# Tibble for reference line
ref.right <- tibble(Position = names(labs.position), Abundance = 1)
# make plot
p.right <- ra_steps %>%
  ggplot(aes(Position, Abundance, color = Taxon)) +
  geom_path(data = ref.right, aes(group = Abundance), color = "grey") +
  geom_path(aes(group = Taxon)) +
  geom_point() +
  scale_y_log10(limits = ylimits, breaks = c(0.1, 0.3, 1, 3, 10), 
                labels = log_formatter) +
  scale_x_discrete(labels = labs.position, expand = xexpand) +
  geom_rangeframe(color = "black") +
  scale_color_manual(values = colors.taxon, labels = tax_labeller) + 
  # labs(x = "Position in the workflow", y = "Abundance / gm. mean", 
  labs(x = "Biological substance", y = "Abundance / gm. mean", 
       title = "Composition through the workflow") +
  guides(color = guide_legend(reverse = TRUE)) +
  base_theme +
  theme(axis.text.x = element_text(size=8, family = ""))
l <- get_legend(p.right + theme(legend.position = c(0.45, 0.55), 
                                plot.margin = margin(t=5, b=5)))
plot_grid(
  # p.left + theme(plot.margin = margin(l=5, t=5, b=5)),
  # p.right + theme(plot.margin = margin(t=5, b=5)),
  p.left,
  p.right,
  l,
  labels = c("A", "B", NULL),
  nrow = 1, rel_widths = c(1, 1, 0.3))


Bias versus 16S CN
How much PCR bias and total bias is explained by 16S CN variation?
  
  Note: The “Total” variable in bias_steps includes CN correction; here we are interested in the total bias observed in the cell mixtures without CN correction being applied.

For this analysis I’m going to ignore the uncertainty in the estimate of the total bias and the PCR bias, since the standard errors are small relative to the magnitude of the bias,

bias_steps %>%
  filter(Type == "Cells")
## # A tibble: 7 x 4
##   Taxon                    Type   Bias gm_se
##   <chr>                    <chr> <dbl> <dbl>
## 1 Atopobium_vaginae        Cells 0.285  1.04
## 2 Gardnerella_vaginalis    Cells 0.160  1.05
## 3 Lactobacillus_crispatus  Cells 2.29   1.03
## 4 Lactobacillus_iners      Cells 4.68   1.02
## 5 Prevotella_bivia         Cells 1.79   1.04
## 6 Sneathia_amnii           Cells 4.59   1.04
## 7 Streptococcus_agalactiae Cells 0.250  1.03
bias_steps %>%
  filter(Type == "PCR")
## # A tibble: 7 x 4
##   Taxon                    Type   Bias gm_se
##   <chr>                    <chr> <dbl> <dbl>
## 1 Atopobium_vaginae        PCR   1.03   1.05
## 2 Gardnerella_vaginalis    PCR   0.493  1.08
## 3 Lactobacillus_crispatus  PCR   0.591  1.05
## 4 Lactobacillus_iners      PCR   1.94   1.06
## 5 Prevotella_bivia         PCR   0.423  1.06
## 6 Sneathia_amnii           PCR   1.85   1.05
## 7 Streptococcus_agalactiae PCR   2.20   1.04
First let’s get a table with all the various bias estimates as well as CN and genome size, with elements centered (divided by the geometric mean over taxa).

tb <- bias_steps %>%
  select(-gm_se) %>%
  spread(Type, Bias) %>%
  left_join(species_info, by = "Taxon") %>%
  select(Taxon, Cells, PCR, Copy_number, Genome_size) %>%
  mutate_if(is.numeric, center_elts)
tb
## # A tibble: 7 x 5
##   Taxon                    Cells   PCR Copy_number Genome_size
##   <chr>                    <dbl> <dbl>       <dbl>       <dbl>
## 1 Atopobium_vaginae        0.285 1.03        0.568       0.837
## 2 Gardnerella_vaginalis    0.160 0.493       0.568       0.954
## 3 Lactobacillus_crispatus  2.29  0.591       1.14        1.19 
## 4 Lactobacillus_iners      4.68  1.94        1.42        0.742
## 5 Prevotella_bivia         1.79  0.423       1.14        1.46 
## 6 Sneathia_amnii           4.59  1.85        0.852       0.773
## 7 Streptococcus_agalactiae 0.250 2.20        1.99        1.26
Next, we compute the predicted bias based on copy-number variation. The expected PCR bias is given by the CN per bp for the PCR step, and by the CN per genome for the total bias measured in the cell mixtures.

tb <- tb %>%
  mutate(
    Cells.pred = Copy_number, 
    PCR.pred = Copy_number / Genome_size,
    Cells.resid = Cells / Cells.pred,
    PCR.resid = PCR / PCR.pred
  )
How much of the bias is explained by CN variation? We will quantify the bias and the residual bias by the squared error as measured by the Aitchison norm, and compute a coefficient of determination quantifying how much of the bias is explained by the predicted bias from CN.

an2 <- tb %>%
  summarize_if(is.numeric, ~anorm(.)^2) %>% 
  unlist %>%
  as.list
The coefficient of determination for CN variation as an explanation of PCR bias is

1 - an2$PCR.resid / an2$PCR
## [1] 0.5998252
and of total bias is

1 - an2$Cells.resid / an2$Cells
## [1] 0.0991634
These numbers are equivalant to computing the standard R^2 measure on the log-transformed bias vectors. E.g.,

tb0 <- tb %>%
  mutate_if(is.numeric, log)
mu <- mean(tb0$PCR)
1 - sum((tb0$PCR - tb0$PCR.pred)^2) / sum((tb0$PCR - mu)^2)
## [1] 0.5998252
Supplemental figure plotting bias against copy number:
  
  p.pcr <- tb %>%
  mutate(Taxon.abbrev = tax_abbrev(Taxon)) %>%
  ggplot(aes(Copy_number / Genome_size, PCR, label = Taxon.abbrev)) +
  geom_abline() +
  geom_text() +
  scale_x_log10(expand = expand_scale(mult = 0.1)) +
  scale_y_log10() +
  geom_rangeframe() +
  coord_fixed() +
  labs(title = "PCR bias vs. copy number", y = "Efficiency / gm. mean", 
       x = "16S copies per bp / gm. mean") +
  base_theme
p.tot <- tb %>%
  mutate(Taxon.abbrev = tax_abbrev(Taxon)) %>%
  ggplot(aes(Copy_number, Cells, label = Taxon.abbrev)) +
  geom_abline() +
  geom_text() +
  scale_x_log10(expand = expand_scale(mult = 0.25)) +
  scale_y_log10() +
  geom_rangeframe() +
  coord_fixed() +
  labs(title = "Total bias vs. copy number", y = "Efficiency / gm. mean", 
       x = "16S copies per genome / gm. mean") +
  base_theme
plot_grid(p.pcr, p.tot, nrow = 1, labels = c("A", "B"), rel_widths = c(1, 1))


Permutation tests
I use a permutation test to see whether the fraction of the bias explained by copy number could be due to chance, for PCR bias and the total (Cells) bias. In particular, I compute the p-value as the fraction of permutations of the predicted bias that gives a squared error (residual bias, measured by the Aitchison norm) less than or equal to the observed error,

K <- nrow(tb)
perms <- gtools::permutations(K, K)
# the residual bias we observe
pcr_obs <- anorm(tb$PCR / tb$PCR.pred)^2
cells_obs <- anorm(tb$Cells / tb$Cells.pred)^2
# the residual bias for each permutation
pcr_perms  <- apply(perms, 1, 
                    function (x) { anorm(tb$PCR / tb$PCR.pred[x])^2 })
cells_perms  <- apply(perms, 1, 
                      function (x) { anorm(tb$Cells / tb$Cells.pred[x])^2 })
# the p values
pcr_pval <- mean(pcr_perms <= pcr_obs) %>% round(3)
cells_pval <- mean(cells_perms <= cells_obs) %>% round(3)
p.pcr <- qplot(pcr_perms) + 
  geom_vline(xintercept = pcr_obs) +
  labs(title = paste0("PCR bias vs. copy number; p = ", pcr_pval), 
       x = "Squared Aitchison error")
p.cells <- qplot(cells_perms) + 
  geom_vline(xintercept = cells_obs) +
  labs(title = paste0("Total bias vs. copy number; p = ", cells_pval), 
       x = "Squared Aitchison error")
plot_grid(p.pcr, p.cells)


Table with estimated bias and summary statistics
Compute pairwise-error summary statistics
Mixture experiments: Bias and noise
Get the ratios with the estimated bias and the residuals,

ratios1 <- main0 %>%
  select(Sample, Taxon, Mixture_type, Observed, Actual, Error, Bias,
         Residual) %>%
  mutate(Taxon = tax_abbrev(Taxon)) %>%
  compute_ratios(group_vars = c("Mixture_type", "Sample")) %>%
  filter(!is.na(Actual)) %>%
  mutate(Pair = paste(Taxon.x, Taxon.y, sep = ":"))
Compute pairwise error summary statistics:
  
  avg_errors <- ratios1 %>%
  # Take the geometric absolute values of the errors
  mutate_at(vars(Error, Bias, Residual), gm_abs) %>%
  # Average the errors for each pair of different taxa
  filter(Taxon.x != Taxon.y) %>%
  group_by(Mixture_type, Taxon.x, Taxon.y) %>%
  summarize_at(vars(Error, Bias, Residual), gm_mean) %>%
  # Average over pairs of taxa for each mixture type
  group_by(Mixture_type) %>%
  summarize_at(vars(Error, Bias, Residual), gm_mean) %>%
  rename_at(vars(Error, Bias, Residual), paste0, ".avg")
max_errors <- ratios1 %>%
  # Take the geometric absolute values of the errors
  mutate_at(vars(Error, Bias, Residual), gm_abs) %>%
  # Find the max fold error in a ratio for each mixture type
  group_by(Mixture_type) %>%
  summarize_at(vars(Error, Bias, Residual), gm_range) %>%
  rename_at(vars(Error, Bias, Residual), paste0, ".max")
pw_summary <- bind_cols(avg_errors, max_errors) %>%
  select(Type = Mixture_type, Bias.max, Bias.avg, Residual.avg)
pw_summary
## # A tibble: 3 x 4
##   Type        Bias.max Bias.avg Residual.avg
##   <chr>          <dbl>    <dbl>        <dbl>
## 1 Cells          29.3      5.57         1.21
## 2 DNA             6.10     2.65         1.19
## 3 PCR_product     1.59     1.22         1.27
Workflow steps: Bias only
Similarly, compute pairwise bias summary statistics for each protocol step

bias_steps.ratios <- bias_steps %>%
  filter(Type %in% c("Extraction", "PCR", "Sequencing")) %>%
  select(Type, Taxon, Bias) %>%
  compute_ratios(group_vars = c("Type"))
avg_bias_steps <- bias_steps.ratios %>%
  mutate_at(vars(Bias), gm_abs) %>%
  filter(Taxon.x != Taxon.y) %>%
  group_by(Type, Taxon.x, Taxon.y) %>%
  summarize_at(vars(Bias), gm_mean) %>%
  group_by(Type) %>%
  summarize_at(vars(Bias), gm_mean) %>%
  rename_at(vars(Bias), paste0, ".avg")
max_bias_steps <- bias_steps.ratios %>%
  mutate_at(vars(Bias), gm_abs) %>%
  group_by(Type) %>%
  summarize_at(vars(Bias), gm_range) %>%
  rename_at(vars(Bias), paste0, ".max")
pw_summary_steps <- bind_cols(avg_bias_steps, max_bias_steps) %>%
  select(Type, Bias.max, Bias.avg)
pw_summary_steps
## # A tibble: 3 x 3
##   Type       Bias.max Bias.avg
##   <chr>         <dbl>    <dbl>
## 1 Extraction    36.6      5.53
## 2 PCR            5.21     2.32
## 3 Sequencing     1.59     1.22
Combine with bias estimates in a single table
# Function to create string for formatting numbers with glue::glue()
fmt_string <- function(var, digits) {
  paste0("{format(", var, ", digits = ", digits, ")}")
}
bias_tab <- bias_steps %>%
  select(-gm_se) %>%
  mutate(Bias = glue::glue(fmt_string("Bias", 1))) %>%
  mutate_all(as.character) %>%
  spread(Type, Bias) %>%
  select(Taxon, Cells, DNA, PCR_product, Extraction, PCR, Sequencing)
pw_summary_tab <- bind_rows(pw_summary, pw_summary_steps) %>%
  mutate(
    Bias.max = glue::glue(fmt_string("Bias.max", 2)),
    Bias.avg = glue::glue(fmt_string("Bias.avg", 2)),
    Residual.avg = glue::glue(fmt_string("Residual.avg", 2))
  ) %>%
  mutate_all(as.character) %>%
  gather("Taxon", "Value", -Type) %>%
  spread(Type, Value) %>%
  mutate_all(~ifelse(str_detect(., "NA"), "---", .))
# Order taxa from highest to lowest efficiency in the cell mixtures
taxa.ranked <- bias %>%
  filter(Mixture_type == "Cells") %>%
  arrange(-Bias) %>%
  {.$Taxon}
lvls.taxon <- c(taxa.ranked, "Bias.max", "Bias.avg", "Residual.avg")
lbls.taxon <- c(taxa.ranked, "Max pairwise bias", 
                "Avg. pairwise bias", "Avg. pairwise noise")
bias_tab0 <- bind_rows(bias_tab, pw_summary_tab) %>%
  mutate(
    Taxon = factor(Taxon, lvls.taxon, lbls.taxon),
  ) %>%
  arrange(Taxon)
bias_tab0
## # A tibble: 10 x 7
##    Taxon               Cells DNA    PCR_product Extraction PCR   Sequencing
##    <fct>               <chr> <chr>  <chr>       <chr>      <chr> <chr>     
##  1 Lactobacillus_iners 4.7   2.3    1.2         2.0        1.9   1.2       
##  2 Sneathia_amnii      4.6   2.4    1.3         1.9        1.8   1.3       
##  3 Lactobacillus_cris… 2.3   0.5    0.9         4.3        0.6   0.9       
##  4 Prevotella_bivia    1.8   0.4    0.9         4.6        0.4   0.9       
##  5 Atopobium_vaginae   0.3   1.1    1.0         0.3        1.0   1.0       
##  6 Streptococcus_agal… 0.2   2.0    0.9         0.1        2.2   0.9       
##  7 Gardnerella_vagina… 0.2   0.4    0.8         0.4        0.5   0.8       
##  8 Max pairwise bias   29.3  " 6.1" " 1.6"      36.6       " 5.… " 1.6"    
##  9 Avg. pairwise bias  5.6   2.7    1.2         5.5        2.3   1.2       
## 10 Avg. pairwise noise 1.2   1.2    1.3         ---        ---   ---
Latex version for the manuscript:
  
  tex <- bias_tab0 %>%
  rename(`PCR prod.` = PCR_product, `Seq. + Inf.` = Sequencing) %>%
  mutate(
    Taxon = kableExtra::cell_spec(
      str_replace(Taxon, "_", " "), 
      "latex", 
      italic = Taxon %in% taxa)
  ) %>%
  knitr::kable(format="latex", booktabs = TRUE, linesep = "",
               escape = FALSE, align = c("l", rep("r", 6))) %>%
  kableExtra::add_header_above(c(" ", "Mixtures" = 3, 
                                 "Steps" = 3)) %>%
  kableExtra::row_spec(length(taxa), extra_latex_after = "\\midrule")
# tex
# clipr::write_clip(tex)
Calibration
To illustrate how one can calibrate samples with compositions different from those used to measure bias, we estimate bias from just the 7-species samples and use it to correct all samples.

Estimate bias as before,

bias.7sp <- main0 %>%
  filter(Num_species == 7) %>%
  select(Mixture_type, Sample, Taxon, Observed, Actual) %>%
  mutate(Error = Observed / Actual) %>%
  group_by(Mixture_type) %>%
  nest %>%
  mutate(
    Error_matrix = map(data, build_matrix, 
                       Sample, Taxon, Error, fill = NaN),
    Estimate = map(Error_matrix, center, enframe = TRUE),
  ) %>%
  unnest(Estimate) %>%
  rename(Bias_7sp = Center)
Compare the two bias estimates:
  
  left_join(bias, bias.7sp, by = c("Mixture_type", "Taxon")) %>%
  filter(Mixture_type == "Cells")
## # A tibble: 7 x 5
##   Mixture_type Taxon                     Bias gm_se Bias_7sp
##   <chr>        <chr>                    <dbl> <dbl>    <dbl>
## 1 Cells        Atopobium_vaginae        0.285  1.04    0.377
## 2 Cells        Gardnerella_vaginalis    0.160  1.05    0.223
## 3 Cells        Lactobacillus_crispatus  2.29   1.03    2.10 
## 4 Cells        Lactobacillus_iners      4.68   1.02    4.02 
## 5 Cells        Prevotella_bivia         1.79   1.04    1.38 
## 6 Cells        Sneathia_amnii           4.59   1.04    4.23 
## 7 Cells        Streptococcus_agalactiae 0.250  1.03    0.240
Get calibrated compositions from both:
  
  cal <- main0 %>%
  left_join(bias.7sp, by = c("Mixture_type", "Taxon")) %>%
  mutate_by(Sample,
            Calibrated_all = close_elts(Observed / Bias),
            Calibrated_7sp = close_elts(Observed / Bias_7sp)
  ) %>%
  gather("Estimate_type", "Estimate", 
         Observed, Calibrated_all, Calibrated_7sp) %>%
  mutate(Estimate_type = factor(Estimate_type, 
                                c("Observed", "Calibrated_all", "Calibrated_7sp")))
Check the reduction in error in proportions,

ggplot(cal, aes(Taxon, odds(Estimate) / odds(Actual), 
                color = as.factor(Num_species))) +
  geom_quasirandom() +
  scale_y_log10() + 
  facet_grid(Mixture_type ~ Estimate_type) +
  scale_color_manual(values = colors.num_species) +
  tax_theme


As expected, the improvement isn’t as great but is still substantial for the cell and DNA mixtures. In the cell mixtures, there seems to be a systematic underestimate of G. vaginalis and A. vaginae and a corresponding over-estimate of others such as P bivia. This could be due to chance over-representation of Gv and Av in the 7sp mixtures, e.g. due to sample construction error.

Check that calibration reduces the error in the proportions on the samples not used to estimate bias:
  
  cal_error.props <- cal %>%
  filter(Num_species < 7) %>%
  group_by(Mixture_type, Estimate_type) %>%
  summarize(
    MSE.prop = mean((Estimate - Actual)^2),
    MSE.logit = mean((logit(Estimate) - logit(Actual))^2),
    RMSE.logit = sqrt(mean( (logit(Estimate) - logit(Actual))^2 )),
  )
cal_error.props
## # A tibble: 9 x 5
## # Groups:   Mixture_type [3]
##   Mixture_type Estimate_type  MSE.prop MSE.logit RMSE.logit
##   <chr>        <fct>             <dbl>     <dbl>      <dbl>
## 1 Cells        Observed        0.0906     3.59        1.90 
## 2 Cells        Calibrated_all  0.00237    0.0450      0.212
## 3 Cells        Calibrated_7sp  0.00666    0.139       0.373
## 4 DNA          Observed        0.0445     1.14        1.07 
## 5 DNA          Calibrated_all  0.00212    0.0447      0.211
## 6 DNA          Calibrated_7sp  0.00229    0.0480      0.219
## 7 PCR_product  Observed        0.00619    0.121       0.348
## 8 PCR_product  Calibrated_all  0.00345    0.0693      0.263
## 9 PCR_product  Calibrated_7sp  0.00494    0.0986      0.314
Also check the reduction in the average Bray-Curtis dissimilarity and Aitchison distance to the actual compositions after calibration:
  
  cal_error.dist <- cal %>%
  filter(Num_species < 7) %>%
  group_by(Mixture_type, Estimate_type, Sample) %>%
  summarize(
    Dist.BC = xydist(Estimate, Actual, method = "bray"),
    Dist.Ai = xydist(Estimate, Actual, method = "aitchison", trim = TRUE)
  ) %>%
  summarize_at(vars(Dist.BC, Dist.Ai), mean)
cal_error.dist 
## # A tibble: 9 x 4
## # Groups:   Mixture_type [3]
##   Mixture_type Estimate_type  Dist.BC Dist.Ai
##   <chr>        <fct>            <dbl>   <dbl>
## 1 Cells        Observed        0.348    1.73 
## 2 Cells        Calibrated_all  0.0467   0.166
## 3 Cells        Calibrated_7sp  0.0830   0.305
## 4 DNA          Observed        0.243    0.966
## 5 DNA          Calibrated_all  0.0479   0.175
## 6 DNA          Calibrated_7sp  0.0499   0.181
## 7 PCR_product  Observed        0.0792   0.285
## 8 PCR_product  Calibrated_all  0.0620   0.225
## 9 PCR_product  Calibrated_7sp  0.0752   0.272
Summary for the manuscript: Calibration (using just the 7-species samples) reduced the mean squared error of the proportions in the calibrated (non-7-species) samples by 92.6% and the average Bray-Curtis dissimilarity between the actual and observed compositions from 0.35 to 0.08.

Model comparisons
Starting point for fitting the linear models. For fitting the linear models on the proportions, keep single-species samples and rows where Actual == 0. Will still set Observed == 0 for taxa not supposed to be in the samples (i.e., remove cross-sample contamination).

main1 <- main %>%
  mutate_by(Sample, 
            Observed = close_elts(Count * Expected),
            Actual = close_elts(Expected)) %>%
  select(-Count, -Expected) %>%
  select(Sample, Taxon, Observed, Actual, everything())
Simple linear regression on the proportions
Simple linear model of Actual ~ Observed, without and with intercept term. Each taxon is fit independently.

fits.slm <- main1 %>%
  group_by(Mixture_type, Taxon) %>%
  nest %>%
  mutate(
    fit0 = map(data, ~lm(Observed ~ 0 + Actual, data = .)),
    fit1 = map(data, ~lm(Observed ~ 1 + Actual, data = .)),
  )
The fits in fit0 and fit1 are without and with an intercept term.

preds.slm <- fits.slm %>%
  gather("Model", "Fit", fit0, fit1) %>%
  unnest(data, map(Fit, broom::augment))
ggplot(preds.slm, aes(.fitted, Observed, color = Taxon)) +
  geom_quasirandom() +
  facet_grid(Mixture_type ~ Model) +
  scale_color_manual(values = colors.taxon, labels = tax_labeller) + 
  labs(c = "Model prediction")


The fits are very similar, so we’ll use the simpler model (with no intercept term) from here on.

preds.slm <- preds.slm %>%
  filter(Model == "fit0")
Brooks2015 model
The model is decribed in the Methods section “Experimental design andmixture effect models” of Brooks et al 2015. For a given taxon in a given mixture experiment, the model is specified by the formula

fmla <- as.formula(paste(
  "Observed ~ 0 + (", 
  paste(taxa, collapse = "+"),
  ")^3"))
fmla
## Observed ~ 0 + (Atopobium_vaginae + Gardnerella_vaginalis + Lactobacillus_crispatus + 
##     Lactobacillus_iners + Prevotella_bivia + Sneathia_amnii + 
##     Streptococcus_agalactiae)^3
where Observed is the observed proportion of the given taxon and the taxon names on the right-hand side (the predictors) are the actual proportions of the taxa. To fit this model, we need a data frame with each row corresponding to a given taxon’s observed abundance in a sample, and the actual abundances of all 7 taxa.

actual <- main1 %>%
  select(Sample, Taxon, Actual) %>%
  spread(Taxon, Actual)
tb <- main1 %>%
  left_join(actual, by = "Sample")
tb %>% glimpse
## Observations: 1,680
## Variables: 16
## $ Sample                   [3m[90m<chr>[39m[23m "s1-1", "s1-1", "s1-1", "s1-1", "s1-1",…
## $ Taxon                    [3m[90m<chr>[39m[23m "Atopobium_vaginae", "Gardnerella_vagin…
## $ Observed                 [3m[90m<dbl>[39m[23m 0.00000000, 0.00000000, 0.60167254, 0.0…
## $ Actual                   [3m[90m<dbl>[39m[23m 0.0000000, 0.0000000, 0.3333333, 0.0000…
## $ Plate                    [3m[90m<dbl>[39m[23m 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, …
## $ Barcode                  [3m[90m<dbl>[39m[23m 1, 1, 1, 1, 1, 1, 1, 10, 10, 10, 10, 10…
## $ Mixture_type             [3m[90m<chr>[39m[23m "Cells", "Cells", "Cells", "Cells", "Ce…
## $ Num_species              [3m[90m<int>[39m[23m 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, …
## $ Species_list             [3m[90m<chr>[39m[23m "Lactobacillus_crispatus;Prevotella_biv…
## $ Atopobium_vaginae        [3m[90m<dbl>[39m[23m 0.0000000, 0.0000000, 0.0000000, 0.0000…
## $ Gardnerella_vaginalis    [3m[90m<dbl>[39m[23m 0.0000000, 0.0000000, 0.0000000, 0.0000…
## $ Lactobacillus_crispatus  [3m[90m<dbl>[39m[23m 0.3333333, 0.3333333, 0.3333333, 0.3333…
## $ Lactobacillus_iners      [3m[90m<dbl>[39m[23m 0.0000000, 0.0000000, 0.0000000, 0.0000…
## $ Prevotella_bivia         [3m[90m<dbl>[39m[23m 0.3333333, 0.3333333, 0.3333333, 0.3333…
## $ Sneathia_amnii           [3m[90m<dbl>[39m[23m 0.0000000, 0.0000000, 0.0000000, 0.0000…
## $ Streptococcus_agalactiae [3m[90m<dbl>[39m[23m 0.3333333, 0.3333333, 0.3333333, 0.3333…
We again fit this model with lm for each experiment and taxon, and extract the predictions for plotting.

fits.brooks <- tb %>%
  group_by(Mixture_type, Taxon) %>%
  nest %>%
  mutate( fit = map(data, ~lm(fmla, data = .)) )
preds.brooks <- fits.brooks %>%
  unnest(data, map(fit, broom::augment))
Plot with all three models
Finally, we’ll plot the fits from the simple linear model (no intercept), the Brooks model, and our model. We will restrict to just the samples with at least two species and where the Actual proportion is > 0 in order to filter most of the cases where the linear models prediction proportions that are outside of the [0, 1] interval.

preds <- left_join(
  preds.slm %>% select(Mixture_type, Taxon, Sample, Observed, Actual,
                       Num_species, Species_list, SLM = .fitted),
  preds.brooks %>% select(Taxon, Sample, Brooks = .fitted),
  by = c("Sample", "Taxon")
) %>%
  left_join(
    main0 %>% select(Taxon, Sample, Ours = Predicted),
    by = c("Sample", "Taxon")
  ) %>%
  gather("Model", "Predicted", SLM, Brooks, Ours) %>%
  filter(Actual > 0, Num_species > 1) %>%
  mutate(
    Model = factor(Model, c("SLM", "Brooks", "Ours")),
    Model = fct_recode(Model, 
                       `Simple linear model` = "SLM",
                       `Brooks et al model` = "Brooks",
                       `Our model` = "Ours",
    )
  )
ggplot(preds, aes(logit(Predicted), logit(Observed), color = Taxon)) +
  geom_abline(intercept = 0, slope = 1, color = "grey") +
  geom_jitter(width = 0.1, height = 0) +
  geom_rangeframe(color = "black") + 
  facet_grid(Mixture_type ~ Model) +
  base_theme +
  labs(x = "log-odds(Predicted proportion)", 
       y = "log-odds(Observed proportion)") +
  coord_fixed() +
  scale_color_manual(values = colors.taxon, labels = tax_labeller) + 
  theme(
    panel.spacing.x = unit(1, "lines"),
    legend.position = "bottom",
  )
## Warning in log(x): NaNs produced

## Warning in log(x): NaNs produced

## Warning in log(x): NaNs produced
## Warning: Removed 2 rows containing missing values (geom_point).


The few remaining NaNs arise during the logit transform of the few remaining cases where the Brooks model predicts a proportion less than zero.

preds %>%
  filter(Predicted < 0) %>%
  select(-Species_list)
## # A tibble: 2 x 8
##   Mixture_type Taxon    Sample Observed Actual Num_species Model  Predicted
##   <chr>        <chr>    <chr>     <dbl>  <dbl>       <int> <fct>      <dbl>
## 1 Cells        Gardner… s1-23    0.0141  0.143           7 Brook…  -0.00238
## 2 Cells        Gardner… s2-14    0.0222  0.143           7 Brook…  -0.00238
Session info
sessionInfo()
## R version 3.6.1 (2019-07-05)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Arch Linux
## 
## Matrix products: default
## BLAS:   /usr/lib/libblas.so.3.8.0
## LAPACK: /usr/lib/liblapack.so.3.8.0
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] BiasManuscript_0.0.0.9000 metacal_0.1.0            
##  [3] ggbeeswarm_0.6.0          cowplot_0.9.99           
##  [5] ggthemes_4.2.0            forcats_0.4.0            
##  [7] stringr_1.4.0             dplyr_0.8.3              
##  [9] purrr_0.3.2               readr_1.3.1              
## [11] tidyr_0.8.3               tibble_2.1.3             
## [13] ggplot2_3.2.0             tidyverse_1.2.1          
## [15] here_0.1                  rmarkdown_1.13           
## [17] nvimcom_0.9-82            usethis_1.5.0            
## 
## loaded via a namespace (and not attached):
##  [1] nlme_3.1-140      fs_1.3.1          lubridate_1.7.4  
##  [4] devtools_2.0.2    webshot_0.5.1     httr_1.4.0       
##  [7] rprojroot_1.3-2   useful_1.2.6      tools_3.6.1      
## [10] backports_1.1.4   utf8_1.1.4        R6_2.4.0         
## [13] vipor_0.4.5       lazyeval_0.2.2    colorspace_1.4-1 
## [16] withr_2.1.2       tidyselect_0.2.5  prettyunits_1.0.2
## [19] processx_3.3.1    compiler_3.6.1    cli_1.1.0        
## [22] rvest_0.3.4       xml2_1.2.0        desc_1.2.0       
## [25] labeling_0.3      scales_1.0.0      callr_3.2.0      
## [28] digest_0.6.20     pkgconfig_2.0.2   htmltools_0.3.6  
## [31] sessioninfo_1.1.1 rlang_0.4.0       readxl_1.3.1     
## [34] rstudioapi_0.10   generics_0.0.2    jsonlite_1.6     
## [37] gtools_3.8.1      magrittr_1.5      kableExtra_1.1.0 
## [40] Rcpp_1.0.1        munsell_0.5.0     fansi_0.4.0      
## [43] stringi_1.4.3     yaml_2.2.0        MASS_7.3-51.4    
## [46] pkgbuild_1.0.3    plyr_1.8.4        grid_3.6.1       
## [49] crayon_1.3.4      lattice_0.20-38   haven_2.1.0      
## [52] hms_0.4.2         zeallot_0.1.0     knitr_1.23       
## [55] ps_1.3.0          pillar_1.4.2      igraph_1.2.4.1   
## [58] reshape2_1.4.3    codetools_0.2-16  pkgload_1.0.2    
## [61] glue_1.3.1        evaluate_0.14     remotes_2.0.4    
## [64] modelr_0.1.4      vctrs_0.2.0       testthat_2.1.1   
## [67] cellranger_1.1.0  gtable_0.3.0      assertthat_0.2.1 
## [70] xfun_0.7          broom_0.5.2       viridisLite_0.3.0
## [73] beeswarm_0.2.3    memoise_1.1.0     ellipsis_0.2.0.1