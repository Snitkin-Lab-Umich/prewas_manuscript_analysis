multi_vs_MEAN_model <-
  lm(formula = 100 * NumMultiallelicSite/NumVariantSite ~ MeanSNPDist,
     data = variant_summary)
r2 <- round(summary(multi_vs_MEAN_model)$r.squared, 2)

r_string <- paste("~R^2==~", r2)
x_pos <- max(variant_summary$MeanSNPDist) * 0.85
y_pos <-
  min(variant_summary$NumMultiallelicSite / variant_summary$NumVariantSite)
variant_summary %>%
  ggplot(aes(x = MeanSNPDist, y = NumMultiallelicSite/NumVariantSite)) +
  geom_point(mapping = aes(color = Dataset, size = `Dataset Size (#)`)) +
  theme_bw() +
  geom_smooth(method = "lm",
              se = FALSE,
              color = "black",
              size = default_line_size) +
  scale_color_manual(values = palette) +
  ylab("Multiallelic Sites / All Variant Sites") +
  xlab("Mean Pairwise SNP Distance (BP)") +
  scale_x_continuous(labels = scales::comma) +
  guides(color = FALSE) +
  theme(text = element_text(size = default_font_size)) +
  scale_y_continuous(labels=scales::percent_format(accuracy = 0.1)) +
  annotate("text",
           size = rsquared_font_size,
           x = x_pos,
           y = y_pos,
           label = r_string,
           parse = TRUE) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))
save_as_pdf_eps_png("MEAN_multiallelic_vs_SNP_distance",
                    default_width,
                    default_height)


### Median

multi_vs_MEDIAN_model <-
  lm(formula = 100 * NumMultiallelicSite/NumVariantSite ~ MedianSNPDist,
     data = variant_summary)
r2 <- round(summary(multi_vs_MEDIAN_model)$r.squared, 2)

r_string <- paste("~R^2==~", r2)
x_pos <- max(variant_summary$MedianSNPDist) * 0.85
y_pos <-
  min(variant_summary$NumMultiallelicSite / variant_summary$NumVariantSite)
variant_summary %>%
  ggplot(aes(x = MedianSNPDist, y = NumMultiallelicSite/NumVariantSite)) +
  geom_point(mapping = aes(color = Dataset, size = `Dataset Size (#)`)) +
  theme_bw() +
  geom_smooth(method = "lm",
              se = FALSE,
              color = "black",
              size = default_line_size) +
  scale_color_manual(values = palette) +
  ylab("Multiallelic Sites / All Variant Sites") +
  xlab("Median Pairwise SNP Distance (BP)") +
  scale_x_continuous(labels = scales::comma) +
  guides(color = FALSE) +
  theme(text = element_text(size = default_font_size)) +
  scale_y_continuous(labels=scales::percent_format(accuracy = 0.1)) +
  annotate("text",
           size = rsquared_font_size,
           x = x_pos,
           y = y_pos,
           label = r_string,
           parse = TRUE) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))
save_as_pdf_eps_png("MEDIAN_multiallelic_vs_SNP_distance",
                    default_width,
                    default_height)

