# _targets.R

# setup -------------------------------------------------------------------

devtools::load_all()
library(targets)
library(tarchetypes)

invisible(
  lapply(
    list.files(path = "R", pattern = "\\.R$", full.names = TRUE),
    source
  )
)

conflicted::conflict_prefer("filter", "dplyr")

options(
  tidyverse.quiet = TRUE,
  usethis.quiet = TRUE
)

future::plan(future.callr::callr(workers = future::availableCores() - 1))

# target-specific options
tar_option_set(
  packages = c("tidyverse", "patchwork"),
  # packages = c("tidyverse", "patchwork", "xcms"),
  imports = c("rnaseq.lf.hypoxia.molidustat"),
  format = "qs"
)

# values <- tibble::tibble(
#   polarity = c("positive", "negative"),
#   names = c("pos", "neg")
# )

# list of targets ---------------------------------------------------------

list(

  # dna per cell ------------------------------------------------------------

  tar_target(
    dna_per_cell_file,
    path_to_data("dna-per-cell-number.xlsx"),
    format = "file"
  ),
  tar_target(
    dna_per_cell_raw,
    clean_dna_per_cell(dna_per_cell_file)
  ),
  tar_target(
    dna_per_cell_std,
    make_std_curves(dna_per_cell_raw)
  ),
  tar_target(
    dna_per_cell_clean,
    interp_data(dna_per_cell_raw, dna_per_cell_std)
  ),
  tar_target(
    cells_per_dna,
    calculate_cells_per_dna(dna_per_cell_clean)
  ),
  tar_render(
    dna_per_cell_report,
    path = path_to_reports("dna-per-cell.Rmd"),
    output_dir = system.file("analysis/pdfs", package = "Copeland.2021.hypoxia.flux")
  ),

  # extracellular fluxes ----------------------------------------------------

  tar_target(
    fluxes_meta_files,
    path_to_data("(lf|pasmc)_.*_meta\\.csv"),
    format = "file"
  ),
  tar_target(
    fluxes_meta,
    clean_flux_meta(fluxes_meta_files)
  ),
  tar_target(
    fluxes_data_files,
    path_to_data("(lf|pasmc)_.*_[a-z]_\\d{4}-\\d{2}-\\d{2}\\.xlsx"),
    format = "file"
  ),
  tar_target(
    fluxes_data,
    assemble_flux_data(fluxes_data_files)
  ),
  tar_target(
    conc_raw,
    clean_fluxes(fluxes_data)
  ),
  tar_target(
    conc_std,
    make_std_curves(conc_raw)
  ),
  tar_target(
    conc_std_plots,
    print_plots(conc_std$plots, conc_std$title, "fluxes/01_standard_curves"),
    format = "file"
  ),
  tar_target(
    conc_std_clean_fld,
    make_std_curves(dplyr::filter(conc_raw, !(detector == "fld" & conc > 900)))
  ),
  tar_target(
    conc_std_clean,
    clean_flux_std(conc_raw)
  ),
  tar_target(
    conc_interp,
    interp_data(conc_raw, conc_std_clean)
  ),
  tar_target(
    conc_with_missing,
    fill_missing_fluxes(conc_interp, fluxes_meta)
  ),
  tar_target(
    conc_clean,
    filter_assays(conc_with_missing)
  ),
  tar_target(
    evap_raw,
    assemble_evap_data(fluxes_data)
  ),
  tar_target(
    evap_clean,
    fill_missing_evap(evap_raw, conc_clean)
  ),
  tar_target(
    flux_measurements,
    assemble_flux_measurements(conc_clean, evap_clean)
  ),
  tar_target(
    growth_curves,
    plot_growth_curves(flux_measurements)
  ),
  tar_target(
    growth_curve_plots,
    print_plots(growth_curves$plots, growth_curves$title, "fluxes/02_growth_curves"),
    format = "file"
  ),
  tar_target(
    growth_rates,
    calculate_growth_rates(growth_curves)
  ),
  tar_target(
    degradation_curves,
    plot_degradation_curves(flux_measurements)
  ),
  tar_target(
    degradation_curve_plots,
    print_plots(degradation_curves$plots, degradation_curves$title, "fluxes/03_degradation_curves"),
    format = "file"
  ),
  tar_target(
    degradation_rates,
    calculate_degradation_rates(degradation_curves)
  ),
  tar_target(
    k,
    clean_degradation_rates(degradation_rates)
  ),
  tar_target(
    mass_curves,
    plot_mass_curves(flux_measurements)
  ),
  tar_target(
    mass_curve_plots,
    print_plots(mass_curves$plots, mass_curves$title, "fluxes/04_mass_curves"),
    format = "file"
  ),
  tar_target(
    flux_curves,
    plot_flux_curves(mass_curves, k, growth_rates)
  ),
  tar_target(
    flux_curve_plots,
    print_plots(flux_curves$plots, flux_curves$title, "fluxes/05_flux_curves"),
    format = "file"
  ),
  tar_target(
    fluxes,
    calculate_fluxes(flux_curves)
  ),
  # tar_target(
  #   fluxes_pairwise_annot,
  #   annot_pairwise(fluxes)
  # ),
  tar_render(
    extracellular_fluxes_report,
    path = path_to_reports("extracellular-fluxes.Rmd"),
    output_dir = system.file("analysis/pdfs", package = "Copeland.2021.hypoxia.flux")
  ),

  # q bias correction -------------------------------------------------------

  tar_target(
    qbias_files,
    path_to_data("q-bias-correction"),
    format = "file"
  ),
  tar_target(
    qbias_ratios,
    import_qbias(qbias_files)
  ),
  tar_target(
    pred_ratios,
    calculate_predicted_ratios()
  ),
  tar_target(
    correction_factors,
    calculate_correction_factors(qbias_ratios, pred_ratios)
  ),
  tar_render(
    qbias_correction_factor_report,
    path = path_to_reports("qbias-correction-factors.Rmd"),
    output_dir = system.file("analysis/pdfs", package = "Copeland.2021.hypoxia.flux")
  ),

  # mass isotopomer distributions -------------------------------------------

  tar_target(
    mid_files,
    path_to_data("(a|b)_(fs|sim)_(lf|pasmc)_.*\\.csv"),
    format = "file"
  ),
  tar_target(
    mid_clean,
    clean_mids(mid_files)
  ),
  tar_target(
    mid_correct,
    correct_mid(mid_clean)
  ),
  tar_target(
    mids,
    remove_mid_outliers(mid_correct)
  ),
  tar_target(
    mid_curves,
    plot_mid_curves(mids)
  ),
  tar_target(
    mid_curve_plots,
    print_plots(mid_curves$plots, mid_curves$title, "mids"),
    format = "file"
  ),
  tar_render(
    mid_report,
    path = path_to_reports("mass-isotope-distributions.Rmd"),
    output_dir = system.file("analysis/pdfs", package = "Copeland.2021.hypoxia.flux")
  ),
  tar_target(
    pasmc_m5,
    get_m5_citrate(pruned_mids)
  ),

  # biomass -----------------------------------------------------------------

  tar_target(
    biomass_file,
    path_to_data("cell-mass.csv"),
    format = "file"
  ),
  tar_target(
    biomass_clean,
    clean_biomass(biomass_file)
  ),
  tar_target(
    biomass,
    calculate_biomass(biomass_clean)
  ),
  tar_target(
    biomass_equations,
    calculate_biomass_equations(biomass)
  ),
  tar_target(
    biomass_equations_out,
    write_matlab_input(biomass_equations, coefs, "_biomass.csv"),
    format = "file"
  ),
  tar_render(
    biomass_report,
    path = path_to_reports("biomass.Rmd"),
    output_dir = system.file("analysis/pdfs", package = "Copeland.2021.hypoxia.flux")
  ),

  # matlab input ------------------------------------------------------------

  tar_target(
    reactions_file,
    path_to_reports("modeling/matlab-input/reactions.csv"),
    format = "file"
  ),
  tar_target(
    model_reactions,
    format_reactions(reactions_file)
  ),
  tar_target(
    model_fluxes,
    format_fluxes(growth_rates, fluxes)
  ),
  tar_target(
    model_fluxes_out,
    write_matlab_input(model_fluxes, data, "_fluxes.csv"),
    format = "file"
  ),
  tar_target(
    pruned_mids,
    format_mids(mids)
  ),
  tar_target(
    model_mids,
    summarize_mids(pruned_mids)
  ),
  tar_target(
    model_mids_out,
    write_matlab_input(model_mids, data, "_mids.csv"),
    format = "file"
  ),

  # cell viability ----------------------------------------------------------

  tar_target(
    viability_file,
    path_to_data("cell-viability.csv"),
    format = "file"
  ),
  tar_target(
    viability,
    clean_viability(viability_file)
  ),

  # immunoblots -------------------------------------------------------------

  tar_target(
    blot_files,
    path_to_data("immunoblots"),
    format = "file"
  ),
  tar_target(
    blot_raw,
    read_data(blot_files)
  ),
  tar_target(
    blot_norm,
    normalize_densities(blot_raw)
  ),

  # mrna --------------------------------------------------------------------

  tar_target(
    mrna_files,
    path_to_data("mrna"),
    format = "file"
  ),
  tar_target(
    mrna_raw,
    read_data(mrna_files)
  ),
  tar_target(
    mrna_norm,
    normalize_qpcr(mrna_raw)
  ),

  # model fluxes ------------------------------------------------------------

  tar_target(
    map_flux_files,
    path_to_data("model"),
    format = "file"
  ),
  tar_target(
    map_fluxes,
    clean_model_fluxes(map_flux_files, model_reactions)
  ),
  tar_target(
    map_flux_differences,
    assemble_flux_differences(map_fluxes)
  ),
  tar_render(
    map_flux_difference_report,
    path = path_to_reports("flux-differences.Rmd"),
    output_dir = system.file("analysis/pdfs", package = "Copeland.2021.hypoxia.flux")
  ),

  # NAD assay ---------------------------------------------------------------

  tar_target(
    nad_files,
    path_to_data("nad-assay_.*\\.xlsx"),
    format = "file"
  ),
  tar_target(
    nad_data,
    assemble_flux_data(nad_files)
  ),
  tar_target(
    nad_raw,
    clean_nad(nad_data)
  ),
  tar_target(
    nad_conc_std,
    make_std_curves(nad_raw, fo = ~MASS::rlm(value ~ poly(conc, 2, raw = TRUE), data = .x))
  ),
  tar_target(
    nad_conc_std_plots,
    print_plots(nad_conc_std$plots, nad_conc_std$title, "nad/01_standard_curves"),
    format = "file"
  ),
  tar_target(
    nad_interp,
    interp_data(nad_raw, nad_conc_std)
  ),
  tar_target(
    nad_final,
    finalize_nad(nad_interp, cells_per_dna)
  ),
  tar_target(
    nad_annot,
    annot_nad(nad_final)
  ),

  # RNA-seq -----------------------------------------------------------------

  tar_target(
    dds,
    count_rnaseq()
  ),
  tar_target(
    pca_data,
    vst_rnaseq(dds)
  ),
  tar_target(
    rnaseq_pca,
    plot_rnaseq_pca(pca_data)
  ),
  tar_target(
    rnaseq_different_differences,
    identify_deg(dds, expr((h.dmso - n.dmso) - (n.bay - n.dmso)))
  ),
  tar_target(
    rnaseq_hyp,
    identify_deg(dds, expr((h.dmso - n.dmso)))
  ),
  tar_target(
    rnaseq_bay,
    identify_deg(dds, expr((n.bay - n.dmso)))
  ),
  tar_target(
    rnaseq_hyp_bay,
    identify_deg(dds, expr((h.bay - n.bay)))
  ),
  tar_target(
    rnaseq_volcano,
    plot_rnaseq_volcano(rnaseq_different_differences, xlab = "ΔHypoxia/ΔBAY", gois = c("EPAS1","HDAC9", "P4HA2", "RBM3"))
  ),
  tar_target(
    rnaseq_hyp_volcano,
    plot_rnaseq_volcano(rnaseq_hyp, xlab = "Hypoxia/Normoxia")
  ),
  tar_target(
    rnaseq_bay_volcano,
    plot_rnaseq_volcano(rnaseq_bay, xlab = "BAY/DMSO")
  ),
  tar_target(
    rnaseq_hyp_bay_volcano,
    plot_rnaseq_volcano(rnaseq_hyp_bay, xlab = "Hypoxia/Normoxia")
  ),
  tar_target(
    rnaseq_overlap,
    find_rnaseq_overlap(dds)
  ),
  tar_target(
    rnaseq_venn,
    plot_rnaseq_venn(rnaseq_overlap, "Transcripts")
  ),
  tar_target(
    rnaseq_gsea,
    run_gsea(rnaseq_different_differences)
  ),
  tar_target(
    rnaseq_hyp_gsea,
    run_gsea(rnaseq_hyp)
  ),
  tar_target(
    rnaseq_bay_gsea,
    run_gsea(rnaseq_bay)
  ),
  tar_target(
    rnaseq_hyp_bay_gsea,
    run_gsea(rnaseq_hyp_bay)
  ),
  tar_target(
    rnaseq_goi,
    plot_rnaseq_goi(dds_symbols, c("EPAS1", "P4HA2", "RBM3", "HDAC9"))
  ),
  tar_target(
    rnaseq_gsea_plot,
    plot_gsea(rnaseq_gsea, "HALLMARK", lbls = c("With BAY", "With Hypoxia"), vals = unname(clrs[c(4, 2)]))
  ),
  tar_target(
    rnaseq_hyp_gsea_plot,
    plot_gsea(rnaseq_hyp_gsea, "HALLMARK", lbls = c("Down in Hypoxia", "Up in Hypoxia"), vals = unname(clrs[c(1, 2)]))
  ),
  tar_target(
    rnaseq_bay_gsea_plot,
    plot_gsea(rnaseq_bay_gsea, "HALLMARK", lbls = c("Down in BAY", "Up in BAY"), vals = unname(clrs[c(3, 4)]))
  ),
  tar_target(
    rnaseq_hyp_bay_gsea_plot,
    plot_gsea(rnaseq_hyp_bay_gsea, "HALLMARK", lbls = c("Down in Hypoxia", "Up in Hypoxia"), vals = unname(clrs[c(1, 2)]))
  ),
  tar_target(
    gsea_overlap,
    find_gsea_overlap(rnaseq_hyp_gsea, rnaseq_bay_gsea)
  ),
  tar_target(
    gsea_venn,
    plot_rnaseq_venn(gsea_overlap, "Gene Sets")
  ),
  tar_target(
    unique_symbol_ids,
    get_unique_symbol_ids(dds)
  ),
  tar_target(
    dds_symbols,
    dds_to_symbols(dds, unique_symbol_ids)
  ),
  tar_target(
    rnaseq_tfea,
    run_tfea(rnaseq_different_differences)
  ),
  tar_target(
    rnaseq_tfea_plot,
    plot_tfea(rnaseq_tfea)
  ),
  tar_render(
    rnaseq_report,
    path = path_to_reports("rnaseq.Rmd"),
    output_dir = system.file("analysis/pdfs", package = "Copeland.2021.hypoxia.flux")
  ),

  # metabolomics ------------------------------------------------------------

  tar_target(
    metab_targeted_files,
    path_to_data("lf_05-bay_metabolomics-targeted.xlsx"),
    format = "file"
  ),
  tar_target(
    metab_targeted_raw,
    format_metab_targeted(metab_targeted_files)
  ),
  tar_target(
    metab_targeted_clean,
    remove_missing_metab(metab_targeted_raw) %>%
      correct_drift() %>%
      quality_control() %>%
      impute_missing() %>%
      pqn() %>%
      log_transform()
  ),
  tar_target(
    metab_targeted_pca,
    plot_metab_pca(metab_targeted_clean)
  ),
  tar_target(
    metab_targeted_limma,
    fit_metab_limma(metab_targeted_clean)
  ),
  tar_target(
    metab_different_differences,
    metab_top_table(metab_targeted_clean, metab_targeted_limma, "deltas")
  ),
  tar_target(
    metab_hyp,
    metab_top_table(metab_targeted_clean, metab_targeted_limma, "hyp_on_dmso")
  ),
  tar_target(
    metab_bay,
    metab_top_table(metab_targeted_clean, metab_targeted_limma, "bay_on_norm")
  ),
  tar_target(
    metab_volcano,
    plot_metab_volcano(
      metab_different_differences,
      mois = c("GAP", "2-hydroxyglutarate", "aconitate", "taurine", "hydroxyproline", "GABA"),
      colors = clrs[c(2, 4)],
      xlab = "ΔHypoxia/ΔBAY"
    )
  ),
  tar_target(
    metab_volcano_hyp,
    plot_metab_volcano(
      metab_hyp,
      colors = clrs[2:1],
      xlab = "Hypoxia/Normoxia"
    )
  ),
  tar_target(
    metab_volcano_bay,
    plot_metab_volcano(
      metab_bay,
      colors = clrs[4:3],
      xlab = "BAY/DMSO"
    )
  ),
  tar_target(
    metab_venn,
    plot_metab_venn(metab_hyp, metab_bay)
  ),
  tar_target(
    metab_moi,
    plot_mois(metab_targeted_clean, c("GAP", "2-hydroxyglutarate", "aconitate", "taurine", "hydroxyproline", "GABA"))
  ),
  tar_target(
    metab_msea,
    run_msea(metab_different_differences, metab_pathways)
  ),
  tar_target(
    metab_msea_hyp,
    run_msea(metab_hyp, metab_pathways)
  ),
  tar_target(
    metab_msea_bay,
    run_msea(metab_bay, metab_pathways)
  ),
  tar_target(
    metab_pathways,
    get_metab_pathways()
  ),
  tar_target(
    msea_plot,
    plot_msea(metab_msea, lbls = c("With BAY", "With Hypoxia"), vals = unname(clrs[c(4, 2)]))
  ),
  tar_target(
    msea_hyp_plot,
    plot_msea(metab_msea_hyp, lbls = c("Down in Hypoxia", "Up in Hypoxia"), vals = unname(clrs[c(1, 2)]))
  ),
  tar_target(
    msea_bay_plot,
    plot_msea(metab_msea_bay, lbls = c("Down in BAY", "Up in BAY"), vals = unname(clrs[c(3, 4)]))
  ),
  tar_target(
    leading_edge,
    plot_leading_edge(metab_different_differences, metab_pathways[["(KEGG) Citrate cycle (TCA cycle)"]])
  ),
  tar_render(
    metabolomics_report,
    path = path_to_reports("metabolomics-targeted.Rmd"),
    output_dir = system.file("analysis/pdfs", package = "Copeland.2021.hypoxia.flux")
  ),
  # tar_target(
  #   metab_untargeted_files,
  #   path_to_data(".mzML"),
  #   format = "file"
  # ),
  # tar_target(
  #   metab_untargeted_sample_file,
  #   path_to_data("sample-sheet.csv"),
  #   format = "file"
  # ),
  # tar_target(
  #   metab_untargeted_samples,
  #   format_sample_info(metab_untargeted_sample_file)
  # ),
  # tar_map(
  #   values = values,
  #   names = "names",
  #   tar_target(
  #     raw,
  #     read_metab_raw(metab_untargeted_files, polarity, metab_untargeted_samples)
  #   ),
  #   tar_target(
  #     tics,
  #     get_tics(raw, metab_untargeted_samples)
  #   ),
  #   tar_target(
  #     chromatograms,
  #     plot_chromatograms(tics)
  #   ),
  #   tar_target(
  #     intensities,
  #     plot_intensities(tics)
  #   ),
  #   tar_target(
  #     cwp,
  #     optimize_centwave_params(raw, polarity)
  #   ),
  #   tar_target(
  #     peaks,
  #     findChromPeaks(raw, cwp$best_cwp)
  #   ),
  #   tar_target(
  #     merged,
  #     refineChromPeaks(peaks, param = MergeNeighboringPeaksParam(expandRt = 4, ppm = 2.5))
  #   ),
  #   tar_target(
  #     align_group,
  #     optimize_align_group_params(merged, polarity)
  #   ),
  #   tar_target(
  #     aligned,
  #     adjust_rtime(merged, align_group$best_obi)
  #   ),
  #   tar_target(
  #     grouped,
  #     group_peaks(aligned, align_group$best_density)
  #   ),
  #   tar_target(
  #     filled,
  #     fillChromPeaks(grouped)
  #   )
  # ),

  # MYC ---------------------------------------------------------------------

  tar_target(
    simyc_fluxes,
    combine_fluxes(growth_rates, fluxes, exp = "05-simyc")
  ),
  tar_target(
    simyc_fluxes_annot,
    annot_fluxes_simyc(simyc_fluxes)
  ),
  tar_target(
    simyc_fluxes_growth_plot,
    plot_myc(simyc_fluxes, simyc_fluxes_annot, "growth", "Growth Rate (/h)", x = treatment, fill = oxygen) +
      ggplot2::geom_hline(yintercept = 0, size = 0.25)
  ),
  tar_target(
    simyc_fluxes_lactate_plot,
    plot_myc(simyc_fluxes, simyc_fluxes_annot, "lactate", "Lactate\n(fmol/cell/h)", x = treatment, fill = oxygen)
  ),
  tar_target(
    oemyc_fluxes,
    combine_fluxes(growth_rates, fluxes, exp = "bay-myc")
  ),
  tar_target(
    oemyc_fluxes_annot,
    annot_fluxes_oemyc(oemyc_fluxes)
  ),
  tar_target(
    oemyc_fluxes_growth_plot,
    plot_myc(oemyc_fluxes, oemyc_fluxes_annot, "growth", "Growth Rate (/h)", x = virus, fill = treatment)
  ),
  tar_target(
    oemyc_fluxes_lactate_plot,
    plot_myc(oemyc_fluxes, oemyc_fluxes_annot, "lactate", "Lactate\n(fmol/cell/h)", x = virus, fill = treatment)
  ),

  # M1 ----------------------------------------------------------------------

  tar_target(
    m1a,
    plot_growth_curve(flux_measurements, cell = "lf", exp = c("05", "bay"))
  ),
  tar_target(
    m1b,
    plot_growth_rates(growth_rates, cell = "lf", exp = c("05", "bay"))
  ),
  # tar_target(
  #   m1c,
  #   plot_viability(viability) + ggplot2::coord_cartesian(ylim = c(NA, NA))
  # ),
  tar_target(
    m1c_image,
    path_to_manuscript("figures/images/lf_05_hif1a-ldha-blots.png"),
    format = "file"
  ),
  tar_target(
    m1c,
    plot_blot(m1c_image)
  ),
  tar_target(
    m1d_image,
    path_to_manuscript("figures/images/lf_bay_hif1a-ldha-blots.png"),
    format = "file"
  ),
  tar_target(
    m1d,
    plot_blot(m1d_image)
  ),
  tar_target(
    m1e,
    plot_expression(blot_norm, c("lf_05", "lf_bay"), "hif1a", "HIF-1α protein\n(normalized)")
  ),
  tar_target(
    m1f,
    plot_expression(blot_norm, c("lf_05", "lf_bay"), "ldha", "LDHA protein\n(normalized)")
  ),
  tar_target(
    m1g,
    plot_expression(mrna_norm, c("lf_05", "lf_bay"), "glut1", "GLUT1 mRNA\n(normalized)")
  ),
  tar_target(
    m1h,
    plot_expression(mrna_norm, c("lf_05", "lf_bay"), "ldha", "LDHA mRNA\n(normalized)")
  ),
  tar_target(
    m1i,
    plot_high_fluxes(fluxes, "lf", c("05", "bay"))
  ),
  tar_target(
    m1j,
    plot_low_fluxes(fluxes, "lf", c("05", "bay"))
  ),
  tar_target(
    m1,
    arrange_fluxes(m1a, m1b, m1c, m1d, m1e, m1f, m1g, m1h, m1i, m1j)
  ),
  tar_target(
    m1_figure,
    write_figures(m1, "m1.pdf"),
    format = "file"
  ),

  # S1 ----------------------------------------------------------------------

  tar_target(
    s1a,
    plot_time_lines(viability, y = viability, ylab = "Cell viability (%)", clr = "oxygen")
  ),
  tar_target(
    s1b,
    plot_cells_per_dna(dna_per_cell_clean)
  ),
  tar_target(
    dna_count_hypoxia_file,
    path_to_data("dna-count-hypoxia.csv"),
    format = "file"
  ),
  tar_target(
    dna_count_hypoxia,
    clean_dna_count_hypoxia(dna_count_hypoxia_file)
  ),
  tar_target(
    s1c,
    plot_dna_count_hypoxia(dna_count_hypoxia)
  ),
  tar_target(
    s1d,
    plot_evap_data(evap_clean)
  ),
  tar_target(
    s1e,
    plot_k(degradation_rates, k)
  ),
  tar_target(
    s1,
    arrange_s1(s1a, s1b, s1c, s1d, s1e)
  ),
  tar_target(
    s1_figure,
    write_figures(s1, "s1.pdf"),
    format = "file"
  ),

  # S2 ----------------------------------------------------------------------

  tar_target(
    s2a,
    plot_growth_curve(flux_measurements, cell = "lf", exp = c("02"))
  ),
  tar_target(
    s2b,
    plot_growth_rates(growth_rates, cell = "lf", exp = c("02"))
  ),
  tar_target(
    s2c,
    patchwork::plot_spacer()
  ),
  tar_target(
    s2d_image,
    path_to_manuscript("figures/images/lf_02_hif1a-ldha-blots.png"),
    format = "file"
  ),
  tar_target(
    s2d,
    plot_blot(s2d_image)
  ),
  tar_target(
    s2e,
    plot_expression(blot_norm, c("lf_02"), "hif1a", "HIF-1α protein\n(normalized)")
  ),
  tar_target(
    s2f,
    plot_expression(blot_norm, c("lf_02"), "ldha", "LDHA protein\n(normalized)")
  ),
  tar_target(
    s2g,
    plot_expression(mrna_norm, c("lf_02"), "glut1", "GLUT1 mRNA\n(normalized)")
  ),
  tar_target(
    s2h,
    plot_expression(mrna_norm, c("lf_02"), "ldha", "LDHA mRNA\n(normalized)")
  ),
  tar_target(
    s2i,
    plot_high_fluxes(fluxes, "lf", c("02"))
  ),
  tar_target(
    s2j,
    plot_low_fluxes(fluxes, "lf", c("02"))
  ),
  tar_target(
    s2,
    arrange_fluxes(s2a, s2b, s2c, s2d, s2e, s2f, s2g, s2h, s2i, s2j)
  ),
  tar_target(
    s2_figure,
    write_figures(s2, "s2.pdf"),
    format = "file"
  ),

  # S3 ----------------------------------------------------------------------

  tar_target(
    s3a,
    plot_growth_curve(flux_measurements, cell = "pasmc", exper = c("05"))
  ),
  tar_target(
    s3b,
    plot_growth_rates(growth_rates, cell = "pasmc", exper = c("05"))
  ),
  tar_target(
    s3c,
    patchwork::plot_spacer()
  ),
  tar_target(
    s3d_image,
    path_to_manuscript("figures/images/pasmc_05_hif1a-ldha-blots.png"),
    format = "file"
  ),
  tar_target(
    s3d,
    plot_blot(s3d_image)
  ),
  tar_target(
    s3e,
    plot_expression(blot_norm, c("pasmc_05"), "hif1a", "HIF-1α protein\n(normalized)")
  ),
  tar_target(
    s3f,
    plot_expression(blot_norm, c("pasmc_05"), "ldha", "LDHA protein\n(normalized)")
  ),
  tar_target(
    s3g,
    plot_expression(mrna_norm, c("pasmc_05"), "glut1", "GLUT1 mRNA\n(normalized)")
  ),
  tar_target(
    s3h,
    plot_expression(mrna_norm, c("pasmc_05"), "ldha", "LDHA mRNA\n(normalized)")
  ),
  tar_target(
    s3i,
    plot_high_fluxes(fluxes, "pasmc", c("05"))
  ),
  tar_target(
    s3j,
    plot_low_fluxes(fluxes, "pasmc", c("05"))
  ),
  tar_target(
    s3,
    arrange_fluxes(s3a, s3b, s3c, s3d, s3e, s3f, s3g, s3h, s3i, s3j)
  ),
  tar_target(
    s3_figure,
    write_figures(s3, "s3.pdf"),
    format = "file"
  ),

  # M2 ----------------------------------------------------------------------

  tar_target(
    m2ab,
    plot_labeling_rate(mids)
  ),
  tar_target(
    m2c,
    plot_manuscript_mids(pruned_mids)
  ),
  tar_target(
    m2,
    arrange_m2(m2ab, m2c)
  ),
  tar_target(
    m2_figure,
    write_figures(m2, "m2.pdf"),
    format = "file"
  ),

  # S4 ----------------------------------------------------------------------

  tar_target(
    s4,
    plot_lf_mids(pruned_mids)
  ),
  tar_target(
    s4_figure,
    write_figures(s4, "s4.pdf")
  ),

  # S5 ----------------------------------------------------------------------

  tar_target(
    s5,
    plot_pasmc_mids(pruned_mids)
  ),
  tar_target(
    s5_figure,
    write_figures(s5, "s5.pdf")
  ),

  # S6 ----------------------------------------------------------------------

  tar_target(
    time_course_mids,
    format_time_course_mids(model_mids)
  ),
  tar_target(
    s6a,
    plot_mid_time_course(time_course_mids, "lf", "21%", "None", "plasma")
  ),
  tar_target(
    s6b,
    plot_mid_time_course(time_course_mids, "lf", "0.5%", "None", "viridis")
  ),
  tar_target(
    s6,
    arrange_s6(s6a, s6b)
  ),
  tar_target(
    s6_figure,
    write_figures(s6, "s6.pdf")
  ),

  # S7 ----------------------------------------------------------------------

  tar_target(
    s7a,
    plot_normoxia_network(lf_hypoxia_graph, "LF\nNormoxia")
  ),
  tar_target(
    hypoxia_growth_graph,
    make_graph(map_flux_differences, nodes, cell = "lf", treat = "0.5%", normalizer = "growth")
  ),
  tar_target(
    s7d,
    plot_ratio_network(hypoxia_growth_graph, "Hypoxia/Normoxia\nGrowth Rate Normalized", edges = FALSE)
  ),
  tar_target(
    pasmc_hypoxia_graph,
    make_graph(map_flux_differences, nodes, cell = "pasmc", treat = "21%", normalizer = "none")
  ),
  tar_target(
    s7b,
    plot_normoxia_network(pasmc_hypoxia_graph, "PASMC\nNormoxia")
  ),
  tar_target(
    s7c,
    plot_ratio_network(pasmc_hypoxia_graph, "PASMC\nHypoxia/Normoxia")
  ),
  tar_target(
    s7,
    arrange_s7(s7a, s7b, s7c, s7d)
  ),
  tar_target(
    s7_figure,
    write_figures(s7, "s7.pdf")
  ),

  # M3 ----------------------------------------------------------------------

  tar_target(
    node_file,
    path_to_data("nodes\\.csv"),
    format = "file"
  ),
  tar_target(
    nodes,
    readr::read_csv(node_file)
  ),
  tar_target(
    lf_hypoxia_graph,
    make_graph(map_flux_differences, nodes, cell = "lf", treat = "21%", normalizer = "none")
  ),
  tar_target(
    lf_hypoxia_graph_ratio_plot,
    plot_ratio_network(lf_hypoxia_graph, "Hypoxia/Normoxia")
  ),
  tar_target(
    bay_graph,
    make_graph(map_flux_differences, nodes, cell = "lf", treat = "DMSO", normalizer = "none")
  ),
  tar_target(
    bay_graph_ratio_plot,
    plot_ratio_network(bay_graph, "BAY/DMSO")
  ),
  tar_target(
    m3,
    arrange_m4(lf_hypoxia_graph_ratio_plot, bay_graph_ratio_plot)
  ),
  tar_target(
    m3_figure,
    write_figures(m3, "m3.pdf")
  ),

  # M4 ----------------------------------------------------------------------

  tar_target(
    m4,
    plot_lactate_mids(pruned_mids, "lf")
  ),
  tar_target(
    m4_figure,
    write_figures(m4, "m4.pdf")
  ),

  # M5 ----------------------------------------------------------------------

  tar_target(
    twoby_fluxes,
    analyze_twoby_fluxes(growth_rates, fluxes)
  ),
  tar_target(
    m5a,
    plot_twoby_fluxes(twoby_fluxes$data, twoby_fluxes$annot, "growth", "Growth Rate (/h)")
  ),
  tar_target(
    m5b,
    plot_twoby_fluxes(twoby_fluxes$data, twoby_fluxes$annot, "glucose", "Glucose\n(fmol/cell/h)") + ggplot2::scale_y_reverse()
  ),
  tar_target(
    m5c,
    plot_twoby_fluxes(twoby_fluxes$data, twoby_fluxes$annot, "lactate", "Lactate\n(fmol/cell/h)")
  ),
  tar_target(
    m5g,
    plot_nad(nad_final, nad_annot, "NAD", "NAD\n(nmol/cell)")
  ),
  tar_target(
    m5h,
    plot_nad(nad_final, nad_annot, "NADH", "NADH\n(nmol/cell)")
  ),
  tar_target(
    m5i,
    plot_nad(nad_final, nad_annot, "Ratio", "NADH/NAD ratio")
  ),
  tar_target(
    m5,
    arrange_m5(m5a, m5b, m5c, metab_targeted_pca, metab_volcano, metab_moi, msea_plot, leading_edge, m5g, m5h, m5i)
  ),
  tar_target(
    m5_figure,
    write_figures(m5, "m5.pdf")
  ),

  # S8 ----------------------------------------------------------------------

  tar_target(
    s8,
    arrange_s8(metab_volcano_hyp, metab_volcano_bay, metab_venn, msea_hyp_plot, msea_bay_plot)
  ),
  tar_target(
    s8_figure,
    write_figures(s8, "s8.pdf")
  ),

  # M6 ----------------------------------------------------------------------

  tar_target(
    twoby_densities_annot,
    annot_twoby_densities(blot_norm)
  ),
  tar_target(
    m6,
    arrange_m6(rnaseq_pca, rnaseq_volcano, rnaseq_goi, rnaseq_gsea_plot, rnaseq_tfea_plot)
  ),
  tar_target(
    m6_figure,
    write_figures(m6, "m6.pdf")
  ),

  # S9 ----------------------------------------------------------------------

  tar_target(
    s9,
    arrange_s9(rnaseq_hyp_volcano, rnaseq_bay_volcano, rnaseq_venn, gsea_venn, rnaseq_hyp_gsea_plot, rnaseq_bay_gsea_plot)
  ),
  tar_target(
    s9_figure,
    write_figures(s9, "s9.pdf")
  ),

  # M7 ----------------------------------------------------------------------

  tar_target(
    myc_image_file,
    path_to_manuscript("figures/images/lf_05-bay_myc-blots.png"),
    format = "file"
  ),
  tar_target(
    myc_image,
    plot_blot(myc_image_file, scale = 1, vjust = 0, hjust = 0)
  ),
  tar_target(
    myc_blot_quant,
    plot_twoby_densities(blot_norm, "myc", twoby_densities_annot, "MYC protein\n(normalized)")
  ),
  tar_target(
    simyc_image_file,
    path_to_manuscript("figures/images/lf_05-simyc_myc-blots.png"),
    format = "file"
  ),
  tar_target(
    simyc_image,
    plot_blot(simyc_image_file, scale = 1, vjust = 0, hjust = 0)
  ),
  tar_target(
    oemyc_image_file,
    path_to_manuscript("figures/images/lf_bay-myc_myc-blots.png"),
    format = "file"
  ),
  tar_target(
    oemyc_image,
    plot_blot(oemyc_image_file, scale = 1, vjust = 0, hjust = 0)
  ),
  tar_target(
    m7,
    arrange_m7(myc_image, myc_blot_quant, simyc_image, simyc_fluxes_growth_plot, simyc_fluxes_lactate_plot, oemyc_image, oemyc_fluxes_growth_plot, oemyc_fluxes_lactate_plot)
  ),
  tar_target(
    m7_figure,
    write_figures(m7, "m7.pdf")
  ),

  # resources table ---------------------------------------------------------

  tar_target(
    resources_table,
    create_resources()
  ),

  # flux comparison tables --------------------------------------------------

  tar_target(
    lf_hypoxia_table,
    format_flux_table(map_flux_differences, "lf", "0.5%", " SSR 391.7 [311.2-416.6] (95% CI, 362 DOF)", " SSR 334.3 [311.2-416.6] (95% CI, 362 DOF)")
  ),
  tar_target(
    lf_bay_table,
    format_flux_table(map_flux_differences, "lf", "BAY", " SSR 393.5 [311.2-416.6] (95% CI, 362 DOF)", " SSR 392.4 [308.4-413.4] (95% CI, 359 DOF)")
  ),
  tar_target(
    pasmc_hypoxia_table,
    format_flux_table(map_flux_differences, "pasmc", "0.5%", " SSR 575.6 [499.1-630.6] (95% CI, 563 DOF)", " SSR 521.3 [482.2-611.6] (95% CI, 545 DOF)")
  ),

  # write manuscript --------------------------------------------------------

  tar_target(
    template,
    system.file("manuscript/template.docx", package = "Copeland.2021.hypoxia.flux"),
    format = "file"
  ),
  tar_render(
    manuscript,
    path = path_to_manuscript("manuscript.Rmd"),
    output_dir = path_to_manuscript(""),
    output_format = bookdown::word_document2(
      reference_docx = template,
      df_print = "kable",
      fig_caption = TRUE,
      number_sections = FALSE,
      pandoc_args = c(
        "--lua-filter=scholarly-metadata.lua",
        "--lua-filter=author-info-blocks.lua",
        "--lua-filter=pagebreak.lua"
      )
    )
  ),
  tar_render(
    supplement,
    path = path_to_manuscript("supplement.Rmd"),
    output_dir = path_to_manuscript(""),
    output_format = bookdown::word_document2(
      reference_docx = template,
      df_print = "kable",
      fig_caption = TRUE,
      number_sections = FALSE,
      pandoc_args = c(
        "--lua-filter=scholarly-metadata.lua",
        "--lua-filter=author-info-blocks.lua",
        "--lua-filter=pagebreak.lua",
        "--lua-filter=multiple-bibliographies.lua"
      )
    )
  ),
  NULL
)

