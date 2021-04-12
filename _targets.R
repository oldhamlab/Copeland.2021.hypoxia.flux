# _targets.R

# setup -------------------------------------------------------------------

devtools::load_all()
library(targets)
library(tarchetypes)

src <- list.files("R", pattern = "^_.*\\.R", full.names = TRUE)
invisible(lapply(src, source))

conflicted::conflict_prefer("filter", "dplyr")

options(
  tidyverse.quiet = TRUE,
  usethis.quiet = TRUE
)

future::plan(future::multisession(workers = future::availableCores() - 1))

# target-specific options
tar_option_set(
  packages = c("tidyverse", "patchwork"),
  imports = c("lf.hypoxia.molidustat.rnaseq"),
  format = "qs"
)

# plot setup
clrs <- c(RColorBrewer::brewer.pal(4, "Set1")[1:4], "#08306b")
names(clrs) <- c("21%", "0.5%", "DMSO", "BAY", "0.2%")


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
    plot_degradation_curves(flux_measurement)
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
  tar_target(
    fluxes_pairwise_annot,
    annot_pairwise(fluxes)
  ),
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
    model_mids,
    format_mids(mids)
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
    rnaseq_volcano,
    plot_rnaseq_volcano(dds)
  ),
  tar_target(
    rnaseq_venn,
    plot_rnaseq_venn(dds)
  ),

  # M1 ----------------------------------------------------------------------

  tar_target(
    m1a,
    plot_growth_curve(cell = "lf", exp = "05", clr = "oxygen") + ggplot2::coord_cartesian(ylim = c(0, 200))
  ),
  tar_target(
    m1b,
    plot_growth_rates(cell = "lf", exp = "05", clr = "oxygen")
  ),
  tar_target(
    m1c,
    plot_viability(viability) + ggplot2::coord_cartesian(ylim = c(NA, NA))
  ),
  tar_target(
    m1d_image,
    path_to_manuscript("figures/images/lf_05_hif1a-ldha-blots.png"),
    format = "file"
  ),
  tar_target(
    m1d,
    plot_blot(m1d_image)
  ),
  tar_target(
    m1e,
    plot_densities(blot_norm, "lf_05", "hif1a", "HIF-1α protein\n(normalized)", "oxygen")
  ),
  tar_target(
    m1f,
    plot_densities(blot_norm, "lf_05", "ldha", "LDHA protein\n(normalized)", "oxygen")
  ),
  tar_target(
    m1g,
    plot_mrna(mrna_norm, "lf_05", "glut1", "GLUT1 mRNA\n(normalized)", "oxygen")
  ),
  tar_target(
    m1h,
    plot_mrna(mrna_norm, "lf_05", "ldha", "LDHA mRNA\n(normalized)", "oxygen")
  ),
  tar_target(
    m1i,
    plot_high_fluxes("lf", "05", "oxygen", annot = fluxes_pairwise_annot)
  ),
  tar_target(
    m1j,
    plot_low_fluxes("lf", "05", "oxygen", annot = fluxes_pairwise_annot)
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

  # M2 ----------------------------------------------------------------------

  tar_target(
    m2a,
    plot_growth_curve(cell = "lf", exp = "bay", clr = "treatment") + ggplot2::coord_cartesian(ylim = c(0, 200))
  ),
  tar_target(
    m2b,
    plot_growth_rates(cell = "lf", exp = "bay", clr = "treatment")
  ),
  tar_target(
    m2c,
    patchwork::plot_spacer()
  ),
  tar_target(
    m2d_image,
    path_to_manuscript("figures/images/lf_bay_hif1a-ldha-blots.png"),
    format = "file"
  ),
  tar_target(
    m2d,
    plot_blot(m2d_image)
  ),
  tar_target(
    m2e,
    plot_densities(blot_norm, "lf_bay", "hif1a", "HIF-1α protein\n(normalized)", "treatment")
  ),
  tar_target(
    m2f,
    plot_densities(blot_norm, "lf_bay", "ldha", "LDHA protein\n(normalized)", "treatment")
  ),
  tar_target(
    m2g,
    plot_mrna(mrna_norm, "lf_bay", "glut1", "GLUT1 mRNA\n(normalized)", "treatment")
  ),
  tar_target(
    m2h,
    plot_mrna(mrna_norm, "lf_bay", "ldha", "LDHA mRNA\n(normalized)", "treatment")
  ),
  tar_target(
    m2i,
    plot_high_fluxes("lf", "bay", "treatment", annot = fluxes_pairwise_annot)
  ),
  tar_target(
    m2j,
    plot_low_fluxes("lf", "bay", "treatment", annot = fluxes_pairwise_annot)
  ),
  tar_target(
    m2,
    arrange_fluxes(m2a, m2b, m2c, m2d, m2e, m2f, m2g, m2h, m2i, m2j)
  ),
  tar_target(
    m2_figure,
    write_figures(m2, "m2.pdf"),
    format = "file"
  ),

  # S1 ----------------------------------------------------------------------

  tar_target(
    s1a,
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
    s1b,
    plot_dna_count_hypoxia(dna_count_hypoxia)
  ),
  tar_target(
    s1c,
    plot_evap_data(evap_clean)
  ),
  tar_target(
    s1d,
    plot_k(degradation_rates, k)
  ),
  tar_target(
    s1,
    arrange_s1(s1a, s1b, s1c, s1d)
  ),
  tar_target(
    s1_figure,
    write_figures(s1, "s1.pdf"),
    format = "file"
  ),

  # S2 ----------------------------------------------------------------------

  tar_target(
    s2a,
    plot_growth_curve(cell = "lf", exp = "02", clr = "oxygen") + ggplot2::coord_cartesian(ylim = c(0, 200))
  ),
  tar_target(
    s2b,
    plot_growth_rates(cell = "lf", exp = "02", clr = "oxygen")
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
    plot_densities(blot_norm, "lf_02", "hif1a", "HIF-1α protein\n(normalized)", "oxygen")
  ),
  tar_target(
    s2f,
    plot_densities(blot_norm, "lf_02", "ldha", "LDHA protein\n(normalized)", "oxygen")
  ),
  tar_target(
    s2g,
    plot_mrna(mrna_norm, "lf_02", "glut1", "GLUT1 mRNA\n(normalized)", "oxygen")
  ),
  tar_target(
    s2h,
    plot_mrna(mrna_norm, "lf_02", "ldha", "LDHA mRNA\n(normalized)", "oxygen")
  ),
  tar_target(
    s2i,
    plot_high_fluxes("lf", "02", "oxygen", annot = fluxes_pairwise_annot)
  ),
  tar_target(
    s2j,
    plot_low_fluxes("lf", "02", "oxygen", annot = fluxes_pairwise_annot)
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
    plot_growth_curve(cell = "pasmc", exp = "05", clr = "oxygen") + ggplot2::coord_cartesian(ylim = c(0, 300))
  ),
  tar_target(
    s3b,
    plot_growth_rates(cell = "pasmc", exp = "05", clr = "oxygen")
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
    plot_densities(blot_norm, "pasmc_05", "hif1a", "HIF-1α protein\n(normalized)", "oxygen")
  ),
  tar_target(
    s3f,
    plot_densities(blot_norm, "pasmc_05", "ldha", "LDHA protein\n(normalized)", "oxygen")
  ),
  tar_target(
    s3g,
    plot_mrna(mrna_norm, "pasmc_05", "glut1", "GLUT1 mRNA\n(normalized)", "oxygen")
  ),
  tar_target(
    s3h,
    plot_mrna(mrna_norm, "pasmc_05", "ldha", "LDHA mRNA\n(normalized)", "oxygen")
  ),
  tar_target(
    s3i,
    plot_high_fluxes("pasmc", "05", "oxygen", annot = fluxes_pairwise_annot)
  ),
  tar_target(
    s3j,
    plot_low_fluxes("pasmc", "05", "oxygen", annot = fluxes_pairwise_annot)
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

  # M3 ----------------------------------------------------------------------

  tar_target(
    m3ab,
    plot_labeling_rate(mids)
  ),
  tar_target(
    m3c,
    plot_manuscript_mids(model_mids)
  ),
  tar_target(
    m3,
    arrange_m3(m3ab, m3c)
  ),
  tar_target(
    m3_figure,
    write_figures(m3, "m3.pdf"),
    format = "file"
  ),

  # S4 ----------------------------------------------------------------------

  tar_target(
    s4,
    plot_lf_mids(model_mids)
  ),
  tar_target(
    s4_figure,
    write_figures(s4, "s4.pdf")
  ),

  # S5 ----------------------------------------------------------------------

  tar_target(
    s5,
    plot_pasmc_mids(model_mids)
  ),
  tar_target(
    s5_figure,
    write_figures(s5, "s5.pdf")
  ),

  # M4 ----------------------------------------------------------------------

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
    hypoxia_graph,
    make_graph(map_flux_differences, nodes, treat = "21%", normalizer = "none")
  ),
  tar_target(
    hypoxia_graph_ratio_plot,
    plot_ratio_network(hypoxia_graph, "Hypoxia/Normoxia")
  ),
  tar_target(
    bay_graph,
    make_graph(map_flux_differences, nodes, treat = "DMSO", normalizer = "none")
  ),
  tar_target(
    bay_graph_ratio_plot,
    plot_ratio_network(bay_graph, "BAY/DMSO")
  ),
  tar_target(
    m4c,
    plot_lactate_mids(model_mids, "lf")
  ),
  tar_target(
    m4,
    arrange_m4(hypoxia_graph_ratio_plot, bay_graph_ratio_plot, m4c)
  ),
  tar_target(
    m4_figure,
    write_figures(m4, "m4.pdf")
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
    s6c,
    plot_normoxia_network(hypoxia_graph)
  ),
  tar_target(
    hypoxia_growth_graph,
    make_graph(map_flux_differences, nodes, treat = "0.5%", normalizer = "growth")
  ),
  tar_target(
    s6d,
    plot_ratio_network(hypoxia_growth_graph, "Hypoxia/Normoxia\nGrowth Rate Normalized")
  ),
  tar_target(
    s6,
    arrange_s6(s6a, s6b, s6c, s6d)
  ),
  tar_target(
    s6_figure,
    write_figures(s6, "s6.pdf")
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
  # tar_target(
  #   m5d,
  #   plot_nad(nad_final)
  # ),
  tar_target(
    m5,
    arrange_m5(m5a, m5b, m5c, rnaseq_pca, rnaseq_volcano)
  ),
  tar_target(
    m5_figure,
    write_figures(m5, "m5.pdf")
  )

)

