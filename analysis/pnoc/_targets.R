library(targets)

source("functions.R")

tar_option_set(packages = c("magrittr", "tidyverse"))
options(tidyverse.quiet = TRUE)

# Define cores
future::plan(
  future::multicore,
  workers = future::availableCores()
)

list(
  tar_target(
    root,
    "/beegfs/scratch/bruening_scratch/pklemm/gabaergic_disinhibition_bnst_pnocarc_hfd_hyperphagia"
  ),
  tar_target(
    export,
    paste0(root, "/release/pnoc") %>%
      (function(path) {
        dir.create(path, recursive = TRUE)
        return(path)
      })
  ),
  tar_target(
    export_plots,
    paste0(export, "/plots") %>%
      (function(path) {
        dir.create(path, recursive = TRUE)
        return(path)
      })
  ),
  tar_target(
    path_intervention,
    paste0(root, "/primary_data_v2/intervention")
  ),
  tar_target(
    path_gastric_infusion,
    paste0(root, "/primary_data_v2/gastric_infusion")
  ),
  tar_target(
    path_intervention_positions,
    paste0(root, "/primary_data_v2/intervention/intervention_positions_v2.csv"),
    format = "file"
  ),
  tar_target(
    intervention_positions,
    path_intervention_positions %>%
      readr::read_csv(
        col_types = list(
          mouse_id = readr::col_character(),
          component = readr::col_character()
        )
      )
  ),
  tar_target(
    file_intervention_times,
    paste0(root, "/primary_data_v2/intervention/intervention_times_v2.csv"),
    format = "file"
  ),
  tar_target(
    file_intervention_times_gastric_infusion,
    paste0(root, "/primary_data_v2/intervention_times_gastric_infusion.csv"),
    format = "file"
  ),
  tar_target(
    intervention_times_gastric_infusion,
    file_intervention_times_gastric_infusion %>%
      readr::read_csv(
        col_types = list(
          mouse_id = readr::col_character(),
          condtion = readr::col_character(),
          Cagedfood_Intervention = readr::col_time(format = "%M:%OS")
        )
      ) %>%
      dplyr::mutate(
        Cagedfood_Intervention = Cagedfood_Intervention %>%
          lubridate::as.duration() %>%
          as.numeric()
      )
  ),
  # Get intervention times
  tar_target(
    intervention_times,
    file_intervention_times %>%
      readr::read_csv(
        col_types = list(
          intervention_time = readr::col_time(format = "%M:%OS"),
          mouse_id = readr::col_character()
        )
      ) %>%
      dplyr::mutate(intervention_time = lubridate::as.duration(intervention_time))
  ),
  tar_target(
    neurons_intervention_all,
    load_neurons_intervention(
      path_dat = path_intervention,
      intervention_positions = intervention_positions,
      intervention_times = intervention_times,
      downsample_to_seconds = FALSE,
      verbose = TRUE
    )
  ),
  tar_target(
    neurons_gastric_infusion_all,
    load_neurons_gastric_infusion(
      path_dat = path_gastric_infusion,
      intervention_times = intervention_times_gastric_infusion,
      downsample_to_seconds = FALSE,
      verbose = TRUE
    )
  ),
  tar_target(
    neurons_gastric_infusion,
    neurons_gastric_infusion_all %>%
      # Filter out cagedfood intervention and leave 2 minutes space
      dplyr::filter(second <= (600 + (60 * 7) + 480)) %>%
      add_pattern(
        k = 3,
        pattern_names = c(
          "Up" = "3",
          "Non-responder" = "1",
          "Down" = "2"
        )
      )
  ),
  tar_target(
    gastric_infusion_conditions_narrowed,
    c("ensure_chow", "ensure_HFD")
  ),
  tar_target(
    neurons_gastric_infusion_narrowed,
    neurons_gastric_infusion_all %>%
      # Filter for conditions of interest
      dplyr::filter(condition %in% gastric_infusion_conditions_narrowed) %>%
      dplyr::filter(second <= (600 + (60 * 7) + 480)) %>%
      add_pattern(
        k = 3,
        pattern_names = c(
          "Up" = "3",
          "Non-responder" = "1",
          "Down" = "2"
        )
      )
  ),
  tar_target(
    neurons_gastric_infusion_ensure,
    neurons_gastric_infusion_all %>%
      # Filter for conditions of interest
      dplyr::filter(condition %in% c("ensure_chow", "ensure_HFD")) %>%
      dplyr::filter(second <= (600 + (60 * 7) + 480)) %>%
      add_pattern(
        k = 3,
        pattern_names = c(
          "Up" = "3",
          "Non-responder" = "1",
          "Down" = "2"
        )
      )
  ),
  tar_target(
    limit_post_intervention_time,
    # Set to NULL or "max" to remove
    # 60
    "max"
  ),
  tar_target(
    limit_pre_intervention_time,
    # Set to NULL or "max" to remove
    # 60
    "max"
  ),
  tar_target(
    neurons_intervention,
    neurons_intervention_all %>%
      create_joint_window() %>%
      # Limit Intervention time
      limit_intervention_time(
        limit_post_intervention_time = limit_post_intervention_time,
        limit_pre_intervention_time = limit_pre_intervention_time
      ) %>%
      add_pattern(
        k = 3,
        pattern_names = c(
          "Up" = "3",
          "Non-responder" = "1",
          "Down" = "2"
        )
      )
  ),
  tar_target(
    lineplot_pattern_neurons_intervention_dat,
    neurons_intervention %>%
      dplyr::mutate(second_ceiling = ceiling(second)) %>%
      dplyr::group_by(component_id, second_ceiling, pattern) %>%
      dplyr::summarise(z_score_by_baseline = mean(z_score_by_baseline)) %>%
      dplyr::ungroup() %>%
      dplyr::rename(second = second_ceiling)
  ),
  tar_target(
    lineplot_pattern_neurons_intervention,
    lineplot_pattern_neurons_intervention_dat %>%
      lineplot_meanpattern() +
      ggplot2::geom_vline(
        # xintercept = . %>% dplyr::pull(intervention_time) %>% min() %>% as.numeric(),
        xintercept = 600,
        linetype = "dashed",
        color = "black"
      )
  ),
  tar_target(
    heatmap_intervention_dat,
    neurons_intervention %>%
      dplyr::mutate(second_ceiling = ceiling(second)) %>%
      dplyr::group_by(component_id, second_ceiling, pattern, condition) %>%
      dplyr::summarise(z_score_by_baseline = mean(z_score_by_baseline)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(second_ceiling = as.numeric(second_ceiling))
  ),
  tar_target(
    heatmap_intervention,
    heatmap_intervention_dat %>%
      plot_heatmap()
  ),
  tar_target(
    plot_heatmap_neuron_intervention_paper,
    neurons_intervention %>%
      # dplyr::filter(condition == "refeeding_chow") %>%
      plot_tiled_heatmap() +
      patchwork::plot_annotation(
        title = "Heatmaps of Intervention Groups"
      )
  ),
  tar_target(
    neurons_intervention_pattern_percentage,
    neurons_intervention %>%
      get_pattern_percentage()
  ),
  tar_target(
    plot_neurons_intervention_pattern_percentage,
    neurons_intervention_pattern_percentage %>%
      plot_pattern_percentage() +
      ggplot2::ggtitle("Intervention: Pattern percentage per condition")
  ),
  tar_target(
    mean_intervention_z_scores,
    neurons_intervention %>%
      dplyr::filter(measurement == "intervention") %>%
      dplyr::mutate(second_ceiling = as.numeric(ceiling(second))) %>%
      dplyr::group_by(component_id, second_ceiling, pattern, condition) %>%
      dplyr::summarise(z_score_by_baseline = mean(z_score_by_baseline)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(component_id, pattern, condition) %>%
      dplyr::summarise(mean_z_score_by_baseline = mean(z_score_by_baseline)) %>%
      dplyr::ungroup()
  ),
  tar_target(
    plot_zscores_intervention_per_pattern_cagefood_dat,
    mean_intervention_z_scores %>%
      # Remove two samples that skew the analysis
      dplyr::filter(condition %in% c("cagefood_chow", "cagefood_HFD")) %>%
      dplyr::filter(mean_z_score_by_baseline < 30)
  ),
  tar_target(
    plot_zscores_intervention_per_pattern_cagefood,
    plot_zscores_intervention_per_pattern_cagefood_dat %>%
      dplyr::group_by(pattern) %>%
      dplyr::group_map(function(dat, pattern) {
        current_pattern <- dplyr::pull(pattern)
        ggstatsplot::ggbetweenstats(
          data = dat,
          x = condition,
          y = mean_z_score_by_baseline,
          title = current_pattern
          # type = "nonparametric",
          # p.adjust.method = "none"
        )
      }) %>%
      patchwork::wrap_plots(ncol = 1) +
      patchwork::plot_annotation(
        title = "Mean Z-Score by Baseline for the intervention phase of each condition",
        subtitle = "The down-regulated neurons in refeeding chow seem more down-regulated compared to the rest, but not significantly"
      )
  ),
  tar_target(
    plot_zscores_intervention_per_pattern_refeeding_dat,
    mean_intervention_z_scores %>%
      # Remove two samples that skew the analysis
      dplyr::filter(condition %in% c("refeeding_chow", "refeeding_HFD")) %>%
      dplyr::filter(mean_z_score_by_baseline < 30)
  ),
  tar_target(
    plot_zscores_intervention_per_pattern_refeeding,
    plot_zscores_intervention_per_pattern_refeeding_dat %>%
      dplyr::group_by(pattern) %>%
      dplyr::group_map(function(dat, pattern) {
        current_pattern <- dplyr::pull(pattern)
        ggstatsplot::ggbetweenstats(
          data = dat,
          x = condition,
          y = mean_z_score_by_baseline,
          title = current_pattern
          # pairwise.comparisons = TRUE
          # type = "nonparametric",
          # pairwise.display = "all",
          # p.adjust.method = "none"
        )
      }) %>%
      patchwork::wrap_plots(ncol = 1) +
      patchwork::plot_annotation(
        title = "Mean Z-Score by Baseline for the intervention phase of each condition",
        subtitle = "The down-regulated neurons in refeeding chow seem more down-regulated compared to the rest, but not significantly"
      )
  ),
  tar_target(
    lineplot_pattern_neurons_gastric_infusion_narrowed_dat,
    neurons_gastric_infusion_narrowed %>%
      dplyr::mutate(second_ceiling = ceiling(second)) %>%
      dplyr::group_by(component_id, second_ceiling, pattern) %>%
      dplyr::summarise(z_score_by_baseline = mean(z_score_by_baseline)) %>%
      dplyr::ungroup() %>%
      dplyr::rename(second = second_ceiling)
  ),
  tar_target(
    lineplot_pattern_neurons_gastric_infusion_narrowed,
    lineplot_pattern_neurons_gastric_infusion_narrowed_dat %>%
      lineplot_meanpattern() %>%
      add_gastric_infusion_vlines()
  ),
  tar_target(
    lineplot_pattern_neurons_gastric_infusion_ensure,
    neurons_gastric_infusion_ensure %>%
      dplyr::mutate(second_ceiling = ceiling(second)) %>%
      dplyr::group_by(component_id, second_ceiling, pattern) %>%
      dplyr::summarise(z_score_by_baseline = mean(z_score_by_baseline)) %>%
      dplyr::ungroup() %>%
      dplyr::rename(second = second_ceiling) %>%
      lineplot_meanpattern() %>%
      add_gastric_infusion_vlines()
  ),
  tar_target(
    neurons_gastric_infusion_narrowed_pattern_percentage,
    neurons_gastric_infusion_narrowed %>%
      get_pattern_percentage()
  ),
  tar_target(
    plot_neurons_gastric_infusion_narrowed_pattern_percentage,
    neurons_gastric_infusion_narrowed_pattern_percentage %>%
      plot_pattern_percentage() +
      ggplot2::ggtitle("Gastric Infusion: Pattern percentage per condition")
  ),
  tar_target(
    plot_tiled_heatmap_gastric_infusion_ensure,
    neurons_gastric_infusion_ensure %>%
      plot_tiled_heatmap_gi(add_vlines = FALSE) +
      patchwork::plot_annotation(
        title = "Heatmaps of Gastric Infusion Ensure Groups"
      )
  ),
  tar_target(
    mean_gastric_infusion_narrowed_z_scores_combined,
    neurons_gastric_infusion_narrowed %>%
      dplyr::filter(measurement %in% c("gastric_infusion", "post_gastric_infusion")) %>%
      dplyr::mutate(second_ceiling = as.numeric(ceiling(second))) %>%
      dplyr::group_by(component_id, second_ceiling, pattern, condition) %>%
      dplyr::summarise(z_score_by_baseline = mean(z_score_by_baseline)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(component_id, pattern, condition) %>%
      dplyr::summarise(mean_z_score_by_baseline = mean(z_score_by_baseline)) %>%
      dplyr::ungroup()
  ),
  tar_target(
    plot_zscores_gastric_infusion_narrowed_per_pattern_ensure_combined,
    mean_gastric_infusion_narrowed_z_scores_combined %>%
      dplyr::filter(condition %in% c("ensure_chow", "ensure_HFD")) %>%
      plot_zscores_by_pattern() +
      patchwork::plot_annotation(
        title = "Mean Z-Score by Baseline for the GI + Post-GI phase of each condition",
        subtitle = "Ensure Chow vs Ensure HFD"
      )
  ),
  tar_target(
    lineplot_pattern_neurons_intervention_chow_hfd_dat,
    neurons_intervention %>%
      dplyr::filter(condition %in% c("refeeding_chow", "refeeding_HFD")) %>%
      dplyr::mutate(feeding = ifelse(condition %in% c("refeeding_chow", "refeeding_HFD"), "refeeding", "cagefood")) %>%
      dplyr::mutate(second_ceiling = ceiling(second)) %>%
      dplyr::group_by(component_id, second_ceiling, pattern, feeding) %>%
      dplyr::summarise(z_score_by_baseline = mean(z_score_by_baseline)) %>%
      dplyr::ungroup() %>%
      dplyr::rename(second = second_ceiling)
  ),
  tar_target(
    lineplot_pattern_neurons_intervention_chow_hfd,
    lineplot_pattern_neurons_intervention_chow_hfd_dat %>%
      lineplot_meanpattern() +
      ggplot2::geom_vline(
        # xintercept = . %>% dplyr::pull(intervention_time) %>% min() %>% as.numeric(),
        xintercept = 600,
        linetype = "dashed",
        color = "black"
      ) +
      ggplot2::ggtitle(
        "Mean z-score by baseline separated by cluster"
      ) +
      ggplot2::facet_grid(
        feeding~.
        # scales = "free_y",
        # space = "free_y"
      )
  ),
  tar_target(
    neurons_gastric_infusion_narrowed_normalised_zscore,
    neurons_gastric_infusion_narrowed %>%
      dplyr::group_by(component_id, measurement, condition, measurement_length, pattern) %>%
      dplyr::summarise(
        mean_zscore = mean(z_score)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::distinct()
  ),
  tar_target(
    plot_zscore_gastric_infusion_narrowed_ggbetweenstats,
    neurons_gastric_infusion_narrowed_normalised_zscore %>%
      # Get peak analysis per group
      dplyr::group_by(condition) %>%
      dplyr::group_map(function(dat, condition) {
        current_condition <- dplyr::pull(condition)
        ggstatsplot::ggbetweenstats(
          data = dat,
          x = measurement,
          y = mean_zscore,
          title = current_condition
        )
      }) %>%
      # Wrap everything into one figure
      patchwork::wrap_plots()
  ),
  tar_target(
    plot_neurons_gastric_infusion_narrowed_delta_baseline_vs_postgi_dat,
    neurons_gastric_infusion_narrowed_normalised_zscore %>%
      dplyr::filter(condition %in% c("ensure_chow", "ensure_HFD")) %>%
      dplyr::filter(measurement != "gastric_infusion") %>%
      dplyr::select(component_id, condition, mean_zscore, measurement) %>%
      tidyr::pivot_wider(names_from = measurement, values_from = mean_zscore) %>%
      dplyr::mutate(
        condition =
          condition %>%
          forcats::as_factor() %>%
          forcats::fct_relevel("ensure_chow", "ensure_HFD")
      ) %>%
      dplyr::mutate(
        fc_postgi_vs_baseline = post_gastric_infusion / baseline,
        fc_postgi_vs_baseline_pseudo = fc_postgi_vs_baseline + 2,
        log2_fc_postgi_vs_baseline = log2(fc_postgi_vs_baseline_pseudo),
        delta_postgi_vs_baseline = post_gastric_infusion - baseline
      )
  ),
  tar_target(
  plot_neurons_gastric_infusion_narrowed_delta_baseline_vs_postgi,
    plot_neurons_gastric_infusion_narrowed_delta_baseline_vs_postgi_dat %>%
      ggstatsplot::ggbetweenstats(
        x = condition,
        y = delta_postgi_vs_baseline,
        title = "Delta Post-Gastric Infusion - Baseline"
      )
  ),
  tar_target(
    dat_figure_raw_figure_data_export,
    list(
      "Fig3D" =
        lineplot_pattern_neurons_intervention_chow_hfd_dat %>%
        dplyr::filter(feeding == "refeeding"),
      "Fig3E" =
        neurons_intervention_pattern_percentage %>%
        dplyr::filter(condition %in% c("refeeding_HFD", "refeeding_chow")),
      "Fig3F" =
        heatmap_intervention_dat %>%
        dplyr::filter(condition %in% c("refeeding_HFD", "refeeding_chow")),
      "Fig3G" = 
        plot_zscores_intervention_per_pattern_refeeding_dat,
      "Fig4B" = 
        lineplot_pattern_neurons_gastric_infusion_narrowed_dat,
      "Fig4C" =
        neurons_gastric_infusion_narrowed_pattern_percentage %>%
        dplyr::filter(condition %in% c("ensure_HFD", "ensure_chow")),
      "Fig4D" = 
        neurons_gastric_infusion_narrowed_normalised_zscore %>%
        dplyr::select(component_id, measurement, mean_zscore),
      "Fig4F" = 
        neurons_gastric_infusion_ensure %>%
        dplyr::mutate(second_ceiling = ceiling(second)) %>%
        dplyr::group_by(component_id, second_ceiling, pattern, condition) %>%
        dplyr::summarise(z_score_by_baseline = mean(z_score_by_baseline)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(second_ceiling = as.numeric(second_ceiling)),
      "Fig4E" = 
        plot_neurons_gastric_infusion_narrowed_delta_baseline_vs_postgi_dat,
      "Fig4E_FILTERED" = 
        plot_neurons_gastric_infusion_narrowed_delta_baseline_vs_postgi_dat %>%
        dplyr::select(component_id, condition, delta_postgi_vs_baseline) %>%
        dplyr::arrange(condition),
      "Fig4G" = 
        mean_gastric_infusion_narrowed_z_scores_combined %>%
        dplyr::filter(condition %in% c("ensure_chow", "ensure_HFD")),
      "FigS3B" = 
        lineplot_pattern_neurons_intervention_chow_hfd_dat %>%
        dplyr::filter(feeding == "cagefood"),
      "FigS3C" = 
        neurons_intervention_pattern_percentage %>%
        dplyr::filter(condition %in% c("cagefood_HFD", "cagefood_chow")),
      "FigS3D" = 
        heatmap_intervention_dat %>%
        dplyr::filter(condition %in% c("cagefood_HFD", "cagefood_chow")),
      "FigS3E" = 
        plot_zscores_intervention_per_pattern_cagefood_dat
    )
  ),
  tar_target(
    file_figure_raw_figure_data_export,
    NULL %>%
      (function(null) {
        export_path <-
          file.path(export, "raw_figure_data.xlsx")

        rmyknife::write_xls(
          dat = dat_figure_raw_figure_data_export,
          ExcelFileName = export_path,
          SheetNames = names(dat_figure_raw_figure_data_export)
        )
        return(export_path)
      }),
    format = "file"
  ),
  tar_target(
    sems_paper,
    # Calculate SEM's
    list(
      "Fig3G_Up_Chow" =
        dat_figure_raw_figure_data_export$Fig3G %>%
          dplyr::filter(condition == "refeeding_chow" & pattern == "Up") %>%
          dplyr::pull(mean_z_score_by_baseline) %>%
          std_error(),
      "Fig3G_Up_Hfd" =
        dat_figure_raw_figure_data_export$Fig3G %>%
          dplyr::filter(condition == "refeeding_HFD" & pattern == "Up") %>%
          dplyr::pull(mean_z_score_by_baseline) %>%
          std_error(),
      "Fig3G_NonResponder_Chow" =
        dat_figure_raw_figure_data_export$Fig3G %>%
          dplyr::filter(condition == "refeeding_chow" & pattern == "Non-responder") %>%
          dplyr::pull(mean_z_score_by_baseline) %>%
          std_error(),
      "Fig3G_NonResponder_Hfd" =
        dat_figure_raw_figure_data_export$Fig3G %>%
          dplyr::filter(condition == "refeeding_HFD" & pattern == "Non-responder") %>%
          dplyr::pull(mean_z_score_by_baseline) %>%
          std_error(),
      "Fig3G_Down_Chow" =
        dat_figure_raw_figure_data_export$Fig3G %>%
          dplyr::filter(condition == "refeeding_chow" & pattern == "Down") %>%
          dplyr::pull(mean_z_score_by_baseline) %>%
          std_error(),
      "Fig3G_Down_Hfd" =
        dat_figure_raw_figure_data_export$Fig3G %>%
          dplyr::filter(condition == "refeeding_HFD" & pattern == "Down") %>%
          dplyr::pull(mean_z_score_by_baseline) %>%
          std_error(),
      "Fig4D_Baseline_Chow" =
        neurons_gastric_infusion_narrowed_normalised_zscore %>%
          dplyr::select(component_id, measurement, mean_zscore, condition) %>%
          dplyr::filter(measurement == "baseline" & condition == "ensure_chow") %>%
          dplyr::pull(mean_zscore) %>%
          std_error(),
      "Fig4D_Gastric_Infusion_Chow" =
        neurons_gastric_infusion_narrowed_normalised_zscore %>%
          dplyr::select(component_id, measurement, mean_zscore, condition) %>%
          dplyr::filter(measurement == "gastric_infusion" & condition == "ensure_chow") %>%
          dplyr::pull(mean_zscore) %>%
          std_error(),
      "Fig4D_Post_Gastric_Infusion_Chow" =
        neurons_gastric_infusion_narrowed_normalised_zscore %>%
          dplyr::select(component_id, measurement, mean_zscore, condition) %>%
          dplyr::filter(measurement == "post_gastric_infusion" & condition == "ensure_chow") %>%
          dplyr::pull(mean_zscore) %>%
          std_error(),
      "Fig4D_Baseline_HFD" =
        neurons_gastric_infusion_narrowed_normalised_zscore %>%
          dplyr::select(component_id, measurement, mean_zscore, condition) %>%
          dplyr::filter(measurement == "baseline" & condition == "ensure_HFD") %>%
          dplyr::pull(mean_zscore) %>%
          std_error(),
      "Fig4D_Gastric_Infusion_HFD" =
        neurons_gastric_infusion_narrowed_normalised_zscore %>%
          dplyr::select(component_id, measurement, mean_zscore, condition) %>%
          dplyr::filter(measurement == "gastric_infusion" & condition == "ensure_HFD") %>%
          dplyr::pull(mean_zscore) %>%
          std_error(),
      "Fig4D_Post_Gastric_Infusion_HFD" =
        neurons_gastric_infusion_narrowed_normalised_zscore %>%
          dplyr::select(component_id, measurement, mean_zscore, condition) %>%
          dplyr::filter(measurement == "post_gastric_infusion" & condition == "ensure_HFD") %>%
          dplyr::pull(mean_zscore) %>%
          std_error(),
      "Fig4E_Chow_delta" =
        dat_figure_raw_figure_data_export$Fig4E_FILTERED %>%
          dplyr::filter(condition == "ensure_chow") %>%
          dplyr::pull(delta_postgi_vs_baseline) %>%
          std_error(),
      "Fig4E_HFD_delta" =
        dat_figure_raw_figure_data_export$Fig4E_FILTERED %>%
          dplyr::filter(condition == "ensure_HFD") %>%
          dplyr::pull(delta_postgi_vs_baseline) %>%
          std_error(),
      "Fig4G_Up_Chow" =
        dat_figure_raw_figure_data_export$Fig4G %>%
          dplyr::filter(pattern == "Up" & condition == "ensure_chow") %>%
          dplyr::pull(mean_z_score_by_baseline) %>%
          std_error(),
      "Fig4G_Up_HFD" =
        dat_figure_raw_figure_data_export$Fig4G %>%
          dplyr::filter(pattern == "Up" & condition == "ensure_HFD") %>%
          dplyr::pull(mean_z_score_by_baseline) %>%
          std_error(),
      "Fig4G_Non_responder_Chow" =
        dat_figure_raw_figure_data_export$Fig4G %>%
          dplyr::filter(pattern == "Non-responder" & condition == "ensure_chow") %>%
          dplyr::pull(mean_z_score_by_baseline) %>%
          std_error(),
      "Fig4G_Non_responder_HFD" =
        dat_figure_raw_figure_data_export$Fig4G %>%
          dplyr::filter(pattern == "Non-responder" & condition == "ensure_HFD") %>%
          dplyr::pull(mean_z_score_by_baseline) %>%
          std_error(),
      "Fig4G_Down_Chow" =
        dat_figure_raw_figure_data_export$Fig4G %>%
          dplyr::filter(pattern == "Down" & condition == "ensure_chow") %>%
          dplyr::pull(mean_z_score_by_baseline) %>%
          std_error(),
      "Fig4G_Down_HFD" =
        dat_figure_raw_figure_data_export$Fig4G %>%
          dplyr::filter(pattern == "Down" & condition == "ensure_HFD") %>%
          dplyr::pull(mean_z_score_by_baseline) %>%
          std_error()
    ) %>%
      tibble::as_tibble() %>%
      tidyr::pivot_longer(everything(), names_to = "Plot", values_to = "SEM")
  )
)