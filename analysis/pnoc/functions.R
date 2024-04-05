#' Run PAM/kMedoids
#' @param x data frame
#' @param k Cluster number
run_pam_wrapper <- function(
  x,
  k
) {
  cluster::pam(
    x,
    k = k,
    diss = FALSE,
    medoids = NULL,
    nstart = 100,
    stand = FALSE, cluster.only = FALSE,
    do.swap = TRUE,
    keep.diss = F,
    keep.data = F,
    pamonce = 5,
    trace.lev = 0
  ) %>%
    return()
}

#' Convert to cosine dissimilarity matrix (distance matrix)
#' See https://stats.stackexchange.com/questions/31565/compute-a-cosine-dissimilarity-matrix-in-r
#' @param mat Input matrix
fast_cosine <- function(mat) {
  # normalize by sum of each row (vector)
  sim <- mat / sqrt(rowSums(mat * mat))
  sim <- sim %*% t(sim)
  return(sim)
}

#' Attach pattern column to neuron table
#' @param tbl Table containing columns `component_id`, `z_score_by_baseline`, `second`
#' @param clustering Clustering method used
#' @param k Number of expected cluster
#' @param pattern_names Rename cluster names for better understanding of the plots
add_pattern <- function(
  tbl,
  clustering = "kmedoids",
  k = 3,
  pattern_names = NA,
  name_order = NA
) {
  # tbl <- test
  # Get maximum length of all components
  max_seconds_per_component <-
    tbl %>%
    dplyr::group_by(component_id) %>%
    dplyr::summarise(max_second = max(second)) %>%
    dplyr::ungroup()

  max_second_cutoff <-
    max_seconds_per_component %>%
    dplyr::pull(max_second) %>%
    min()

  max_second_absolute <-
    max_seconds_per_component %>%
    dplyr::pull(max_second) %>%
    max()

  components_above_max_second <-
    max_seconds_per_component %>%
    dplyr::filter(max_second > max_second_cutoff)

  components_within_max_second <-
    max_seconds_per_component %>%
    dplyr::filter(max_second == max_second_cutoff)

  paste0("Filtering for, ", max_second_cutoff, " seconds (", nrow(components_within_max_second), " components of this length) . There are ", nrow(max_seconds_per_component), " components with more traces after that (maximum seconds is ", max_second_absolute, ")") %>%
    message()

  tbl_wide <-
    tbl %>%
    dplyr::filter(second <= max_second_cutoff) %>%
    dplyr::select(component_id, second, z_score_by_baseline) %>%
    tidyr::pivot_wider(names_from = component_id, values_from = z_score_by_baseline)


  if (clustering == "kmeans") {
    clustering_result <-
      t(tbl_wide[-1]) %>%
      # Create distance matrix
      dist(method = "euclidean") %>%
      # Hierarchical clustering
      kmeans(k, nstart = 100, iter.max = 100) %>%
      .$cluster
  } else {
    clustering_result <-
      t(tbl_wide[-1]) %>%
      fast_cosine() %>%
      as.dist() %>%
      run_pam_wrapper(k = k) %>%
      .$clustering
  }

  cluster <-
    clustering_result %>%
    dplyr::bind_rows() %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "component_id", values_to = "pattern") %>%
    dplyr::mutate(pattern = forcats::as_factor(pattern))

  # Rename pattern if name vector provided
  if (sum(is.na(pattern_names)) == 0) {
    cluster <-
      cluster %>%
      dplyr::mutate(
        pattern =
          pattern %>%
          forcats::fct_recode(!!!pattern_names) %>%
          # Change the order to the one provided with the list
          forcats::fct_relevel(names(pattern_names))
      )
  }

  tbl %>%
    dplyr::left_join(cluster, by = "component_id")
}

# https://github.com/stas-g/findPeaks
find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
      z <- i - m + 1
      z <- ifelse(z > 0, z, 1)
      w <- i + m + 1
      w <- ifelse(w < length(x), w, length(x))
      if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}

#' Helper function to test get_tidy_video_intervention
#' @param path_intervention Path to intervention files
#' @examples
#'   tidy_videos <- test_get_tidy_video_intervention()
test_get_tidy_video_intervention <- function(
  path_intervention = "/beegfs/scratch/bruening_scratch/pklemm/2022-10-tamara-pnoc-ca-imaging/primary_data/intervention/10077/2_saline.csv.gz"
) {
  tbl <-
    path_intervention %>%
    readr::read_csv() %>%
    dplyr::select(component, F_denoised, frames) %>%
    dplyr::rename(trace = F_denoised, frame = frames)

  result <- list()

  result$downsampled <-
    get_tidy_video_intervention(
      dat = tbl,
      frames_per_second = 17,
      intervention_second = 840,
      downsample_to_seconds = TRUE
    )
  
  result$not_downsampled <-
    get_tidy_video_intervention(
      dat = tbl,
      frames_per_second = 17,
      intervention_second = 840,
      downsample_to_seconds = FALSE
    )
  return(result)
}

#' Get tidy dataframe of video
#' @param dat Dataframe containing `component`, `trace` and `frame` column
#' @param frames_per_second How many traces/frames make up one second
#' @param intervention_second Second of intervention start
#' @export
#' @return Tidy dataframe with rows component, second, trace_mean_per_second and trace_mean_per_second_by_baseline (z-score but normalized using baseline)
get_tidy_video_intervention <- function(
  dat,
  intervention_second,
  frames_per_second = 17,
  downsample_to_seconds = TRUE
) {
  # frames_per_second <- 17
  # intervention_second <- 840
  # downsample_to_seconds <- TRUE
  # dat <- paste0(path_intervention, "/10077/2_saline.csv.gz") %>% readr::read_csv() %>% dplyr::select(component, F_denoised, frames) %>% dplyr::rename(trace = F_denoised, frame = frames)
  # dat <- glue::glue("{rmd_params$root}/data_raw/heterogeneity/Mouse300/registration_glp1_gip_coagonist/registered_videos_glp1_gip_coagonist_video_0.csv.gz") %>% readr::read_csv()
  temp_dat <-
    dat %>%
    dplyr::mutate(frame = as.double(frame)) %>%
    dplyr::mutate(time = frame / frames_per_second) %>%
    # Create a new variable where miliseconds are always rounded up to the next full second
    dplyr::mutate(second = ceiling(time))

  if (downsample_to_seconds) {
    temp_dat <-
      temp_dat %>%
      # Group component and timestamp together
      dplyr::group_by(component, second) %>%
      # Perfom mean fluorescence per group (per component/second)
      dplyr::summarise(trace_mean_per_second = mean(trace)) %>%
      dplyr::ungroup() %>%
      # Quick hack to make the code below consistent whether we are downsampling to seconds or not
      dplyr::rename(trace = trace_mean_per_second)
  }

  # Calculate z_score based on the baseline
  baseline_metrics <-
    temp_dat %>%
    dplyr::filter(second < intervention_second) %>%
    dplyr::group_by(component) %>%
    dplyr::summarise(
      baseline_sd = sd(trace),
      baseline_mean = mean(trace)
    )
  
  # Add measurement
  temp_dat <-
    temp_dat %>%
    dplyr::mutate(measurement =
      ifelse(second < intervention_second, "baseline", "intervention") %>%
      forcats::as_factor()
    ) %>%
    # Attach measurement_length, which is needed for normalising the peak counts
    dplyr::group_by(measurement) %>%
    dplyr::mutate(
      measurement_length = max(second) - min(second)
    )

  temp_dat <-
    temp_dat %>%
    # Attach these metrics
    dplyr::left_join(
      baseline_metrics,
      by = "component"
    ) %>%
    dplyr::mutate(z_score_by_baseline = (trace - baseline_mean) / baseline_sd) %>%
    dplyr::select(-baseline_sd, -baseline_mean) %>%
    # Calculate z-score
    dplyr::group_by(component) %>%
    dplyr::mutate(z_score = scale(trace)[,1]) %>%
    dplyr::ungroup()
  
  if (downsample_to_seconds) {
    # Create duration from seconds column
    temp_dat <-
      temp_dat %>%
      dplyr::mutate(
        second = second %>%
          lubridate::seconds() %>%
          lubridate::as.duration()
        )
  } else {
    # Create duration from seconds time column
    temp_dat <-
      temp_dat %>%
      dplyr::mutate(
        second = time %>%
          lubridate::seconds() %>%
          lubridate::as.duration()
      ) %>%
      # Remove time
      dplyr::select(-time)
  }

  return(temp_dat)
}

#' Load intervention neurons
#' @param path_intervention Path to intervention files
#' @param intervention_positions Table containing `mouse_id`, `condition`, `component`, `position`
#' @param intervention_times Table containing `mouse_id`, `condition`, `intervention_time`
#' @param downsample_to_seconds Whether to downsample to seconds
#' @param accepting_all_neurons Accept all these neurons per default
#' @param verbose Whether to print messages
load_neurons_intervention <- function(
  path_dat = "/beegfs/scratch/bruening_scratch/pklemm/2022-10-tamara-pnoc-ca-imaging/primary_data/intervention_testing",
  intervention_positions,
  intervention_times,
  downsample_to_seconds = TRUE,
  accepting_all_neurons = c("2_saline", "3_liraglutide", "saline", "liraglutide"),
  no_intervention_time = c("baseline_chow", "baseline_HFD", "baseline"),
  verbose = TRUE
) {
  # path_dat <- "/beegfs/scratch/bruening_scratch/pklemm/2022-10-tamara-pnoc-ca-imaging/primary_data_v2/intervention"
  # intervention_times <- intervention_times_v2
  # intervention_positions <- intervention_positions_v2
  # verbose <- TRUE
  # downsample_to_seconds <- FALSE
  path_dat %>%
    list.dirs(full.names = FALSE, recursive = FALSE) %>%
    # Iterate over all mice
    furrr::future_map(function(current_mouse_id) {
    # purrr::map(function(current_mouse_id) {
      # current_mouse_id <- "11606"
      paste0(path_dat, "/", current_mouse_id) %>%
        list.files(full.names = FALSE, recursive = FALSE) %>%
        # Iterate over all videos associated with a mouse
        furrr::future_map(function(video_path) {
        # purrr::map(function(video_path) {
          # video_path <- "baseline_chow.csv.gz"
          # video_path <- "saline.csv.gz"
          # Get condition from video name path
          current_condition <- stringr::str_remove(video_path, pattern = ".csv.gz$")
          # Get accepted components as vector
          current_accepted_components <-
            intervention_positions %>%
            dplyr::filter(
              condition == current_condition &
              mouse_id == current_mouse_id
            ) %>%
            dplyr::distinct(component) %>%
            dplyr::pull()
          
          if (current_condition %in% no_intervention_time) {
            current_intervention_time <- lubridate::as.duration(9999)
          } else {
            # Get intervention time as period and numeric
            current_intervention_time <-
              intervention_times %>%
              dplyr::filter(
                condition == current_condition &
                mouse_id == current_mouse_id
              ) %>%
              dplyr::pull(intervention_time)
          }
          
          current_intervention_time_numeric <-
            current_intervention_time %>%
            as.numeric()
          
          if (verbose) {
            glue::glue("Mouse {current_mouse_id} - {current_condition} - Intervention time {current_intervention_time} (numeric: {current_intervention_time_numeric})") %>%
              message()
            if (current_condition %in% accepting_all_neurons) {
              "Accept all neurons" %>% message()
            } else {
              glue::glue("{length(current_accepted_components)} accepted components") %>%
                message()
              paste0(capture.output(current_accepted_components), collapse = "\n") %>%
                message()
            }
          }

          # Load video
          full_path <- paste0(path_dat, "/", current_mouse_id, "/", video_path)
          full_path %>%
            readr::read_csv(
              col_types = list(
                component = readr::col_character()
              )
            ) %>%
            dplyr::select(component, F_denoised, frames) %>%
            dplyr::rename(
              trace = F_denoised,
              frame = frames
            ) %>%
            # Filter for accepted components
            dplyr::filter(
              (component %in% current_accepted_components) |
              (current_condition %in% accepting_all_neurons)
            ) %>%
            # Get tidy video table
            get_tidy_video_intervention(
              intervention_second = current_intervention_time_numeric,
              downsample_to_seconds = downsample_to_seconds
            ) %>%
            # Attach metadata
            dplyr::mutate(
              component_id = paste0(component, "_", current_mouse_id, "_", current_condition),
              mouse_id = current_mouse_id,
              intervention_time = current_intervention_time,
              condition = current_condition
            ) %>%
            # Attach is_peak
            dplyr::group_by(component_id) %>%
            # Get a tibble for each group
            dplyr::group_map(function(tbl, group_id) {
              # Get vector of traces
              trace_vector <-
                tbl %>%
                dplyr::pull(trace)
              # Get peaks
              peaks <- find_peaks(trace_vector, 2)
              tbl %>%
                # peaks contains the indices of the peaks, so to use them we need the row numbers
                dplyr::mutate(id = dplyr::row_number()) %>%
                dplyr::mutate(is_peak = id %in% peaks) %>%
                # Re-attach grouping variables since they are removed in group_map
                dplyr::bind_cols(group_id) %>%
                # Remove id column again since it is not needed
                dplyr::select(-id)
            }) %>%
            dplyr::bind_rows() %>%
            dplyr::ungroup()
        })
    }) %>%
    # Concatinate everything into one table
    dplyr::bind_rows() %>%
    # Attach positions
    dplyr::left_join(
      intervention_positions,
      by = c("mouse_id", "condition", "component")
    ) %>%
    # Replace NA positions which were accepted by default
    dplyr::mutate(
      position = ifelse(
        is.na(position),
        "accepted_by_default",
        position
      )
    )
}

#' Visualize mean pattern
#' @param dat dataframe containing second, z_score_by_baseline and pattern
#' @param intervention_time Add dashed intervention time marker
#' @return ggplot2 object
lineplot_meanpattern <- function(dat, intervention_time = NULL) {
  result <-
    ggplot2::ggplot(
      data = dat,
      mapping = ggplot2::aes(
        x = second,
        y = z_score_by_baseline
      )
    ) +
    ggplot2::ylab("Z-Score by baseline") +
    ggplot2::xlab("Time (S)") +
    ggplot2::geom_smooth(
      mapping = ggplot2::aes(color = pattern),
      method = "gam",
      formula = y ~ s(x, bs = "cs"),
      se = TRUE,
      size = 0.5
    )
  
  pattern_count <-
    dat %>%
    dplyr::distinct(pattern) %>%
    nrow()
  
  # Add custom coloring for n = 3
  if (pattern_count == 3) {
    result <-
      result +
      # Add custom coloring
      # ggplot2::scale_color_manual(values = c("red", "grey", "blue", "green", "yellow", "orange", "purple"))
      ggplot2::scale_color_manual(values = c("red", "grey", "blue"))
  }
  
  if (is.numeric(intervention_time)) {
    result <-
      result +
      ggplot2::geom_vline(xintercept = intervention_time, linetype = "dotted")
  }

  return(result)
}

plot_heatmap <- function(tbl) {
  ggplot2::ggplot(
    data = tbl,
    mapping = ggplot2::aes(
      x = second_ceiling,
      y = component_id,
      fill = z_score_by_baseline
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::facet_grid(
      pattern + condition~.,
      scales = "free_y",
      space = "free_y"
    ) +
    # Add lims to fix outlier
    # https://github.com/tidyverse/ggplot2/issues/866
    ggplot2::scale_fill_gradient2(
      limits = c(-4, 4),
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      oob = scales::squish
    ) +
    ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
      # aspect.ratio = 1
    ) +
    ggplot2::theme(strip.text.y.right = ggplot2::element_text(angle = 0))
}

#' Create a joint window for all components so we can compare them
#' @param dat A data frame with columns `component_id`, `second`, `trace` and `intervention_time`
#' @param verbose A logical indicating whether to print messages
#' @return dat with joint analysis window
#' @export
create_joint_window <- function(
  tbl,
  verbose = TRUE
) {
  # tbl <- neurons_intervention
  # verbose <- TRUE

  max_seconds_per_component <-
    tbl %>%
    dplyr::group_by(component_id) %>%
    dplyr::summarise(max_second = max(second)) %>%
    dplyr::ungroup()

  min_length <-
    max_seconds_per_component %>%
    dplyr::pull(max_second) %>%
    min()

  max_length <-
    max_seconds_per_component %>%
    dplyr::pull(max_second) %>%
    max()

  intervention_times_per_component <-
    tbl %>%
    dplyr::group_by(component_id) %>%
    dplyr::summarise(intervention = min(intervention_time)) %>%
    dplyr::ungroup()

  min_intervention <-
    intervention_times_per_component %>%
    dplyr::pull(intervention) %>%
    min()

  max_intervention <-
    intervention_times_per_component %>%
    dplyr::pull(intervention) %>%
    max()

  if (verbose) {
    glue::glue("Using minimum intervention time of {min_intervention} (max in other components {max_intervention}) and minimum length of {round(min_length)} (max in other components {round(max_length)})")
  }

  # Get intervention frame
  second_to_frame <- function(sec) {
    as.numeric(sec) * 17
  }

  # We recalculate the windowed seconds based on the frames to not run into rounding errors
  tbl %>%
    # Ungroup just to be sure
    dplyr::ungroup() %>%
    # Have a consistent time frame around the intervention
    dplyr::mutate(frame_windowed = frame - second_to_frame(intervention_time) + second_to_frame(min_intervention)) %>%
    dplyr::filter(frame_windowed > 0 & frame_windowed < second_to_frame(min_length)) %>%
    dplyr::mutate(
      frame_original = frame,
      frame = frame_windowed,
      second_original = second,
      second = (frame_windowed / 17) %>%
        lubridate::seconds() %>%
        lubridate::as.duration()
    ) %>%
    dplyr::select(frame, second, frame_original, second_original, dplyr::everything()) %>%
    dplyr::mutate(intervention_time = lubridate::as.duration(min_intervention)) %>%
    # Update measurement length to reflect the new lengths
    dplyr::ungroup() %>%
    dplyr::group_by(component_id, condition, measurement) %>%
    dplyr::mutate(
      min_second = min(second),
      max_second = max(second),
      measurement_length = max_second - min_second
    ) %>%
    dplyr::ungroup()

  # tbl %>%
  #   # Have a consistent time frame around the intervention
  #   dplyr::mutate(second_windowed = second - intervention_time + min_intervention) %>%
  #   dplyr::filter(second_windowed >= 0 & second_windowed < min_length) %>%
  #   dplyr::mutate(
  #     second_original = second,
  #     second = second_windowed
  #   ) %>%
  #   dplyr::select(-second_windowed) %>%
  #   dplyr::select(second, second_original, dplyr::everything()) %>%
  #   dplyr::mutate(intervention_time = lubridate::as.duration(min_intervention))

}

add_pattern_based_on_peaks <- function(
  tbl,
  k = 3,
  clustering = "kmeans",
  pattern_names = NA
) {
  tbl_wide <-
    tbl %>%
    dplyr::select(component_id, second_ceiling, peaks) %>%
    tidyr::pivot_wider(names_from = component_id, values_from = peaks)

  if (clustering == "kmeans") {
    clustering_result <-
      t(tbl_wide[-1]) %>%
      # Create distance matrix
      dist(method = "euclidean") %>%
      # Hierarchical clustering
      kmeans(k, nstart = 100, iter.max = 100) %>%
      .$cluster
  } else {
    clustering_result <-
      t(tbl_wide[-1]) %>%
      fast_cosine() %>%
      as.dist() %>%
      run_pam_wrapper(k = k) %>%
      .$clustering
  }

  cluster <-
    clustering_result %>%
    dplyr::bind_rows() %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "component_id", values_to = "pattern") %>%
    dplyr::mutate(pattern = forcats::as_factor(pattern))

  # Rename pattern if name vector provided
  if (sum(is.na(pattern_names)) == 0) {
    cluster <-
      cluster %>%
      dplyr::mutate(
        pattern =
          pattern %>%
          forcats::fct_recode(!!!pattern_names) %>%
          # Change the order to the one provided with the list
          forcats::fct_relevel(names(pattern_names))
      )
  }
  
  tbl %>%
    dplyr::left_join(cluster, by = "component_id")
}

#' Load neurons gastric infusion
#' @param path_dat Path to gastric infusion files
#' @param intervention_times Table containing `Mouse`, `Condition`, `Intervention_Time`, `Cagedfood_Intervention`
#' @param downsample_to_seconds Whether to downsample to seconds
#' @param verbose Whether to print messages
load_neurons_gastric_infusion <- function(
  path_dat = "/beegfs/scratch/bruening_scratch/pklemm/2022-10-tamara-pnoc-ca-imaging/primary_data/intervention_testing",
  intervention_times,
  downsample_to_seconds = FALSE,
  verbose = TRUE
) {
  # path_dat <- "/beegfs/scratch/bruening_scratch/pklemm/2022-10-tamara-pnoc-ca-imaging/primary_data_v2/gastric_infusion"
  path_dat %>%
    list.dirs(full.names = FALSE, recursive = FALSE) %>%
    # Iterate over all mice
    furrr::future_map(function(current_mouse_id) {
      # purrr::map(function(current_mouse_id) {
      # current_mouse_id <- "mouse300"
      paste0(path_dat, "/", current_mouse_id) %>%
        list.files(
          full.names = FALSE,
          recursive = FALSE,
          pattern = "*.csv.gz$"
        ) %>%
        # Iterate over all videos associated with a mouse
        furrr::future_map(function(video_path) {
          # purrr::map(function(video_path) {
          # video_path <- "cagedfood_16hfast.csv.gz"
          # Get condition from video name path
          current_condition <- stringr::str_remove(video_path, pattern = ".csv.gz$")
          # Get intervention time as period and numeric
          current_intervention_time <-
            600
            # intervention_times %>%
            # dplyr::filter(
            #   Condition == current_condition &
            #   Mouse == current_mouse_id
            # ) %>%
            # dplyr::pull(Intervention_Time)
          
          # Gastric infusion is 7 times 60 seconds
          current_gastric_infusion_time_end <-
            current_intervention_time + (60 * 7)
          
          current_cagedfood_intervention_time <-
            intervention_times %>%
            dplyr::filter(
              Condition == current_condition &
              Mouse == current_mouse_id
            ) %>%
            dplyr::pull(Cagedfood_Intervention)
          
          # Tamara noted the intervention time as relative to the start of the post-intervention period
          current_cagedfood_intervention_time <-
            current_cagedfood_intervention_time + current_gastric_infusion_time_end

          if (verbose) {
            glue::glue("Mouse {current_mouse_id} - {current_condition} - Intervention time: {current_intervention_time} - Cagedfood Intervention: {current_cagedfood_intervention_time}") %>%
              message()
          }

          # Load video
          full_path <- paste0(path_dat, "/", current_mouse_id, "/", video_path)
          full_path %>%
            readr::read_csv(
              col_types = list(
                component = readr::col_character()
              )
            ) %>%
            dplyr::select(component, F_denoised, frames) %>%
            dplyr::rename(
              trace = F_denoised,
              frame = frames
            ) %>%
            # Get tidy video table
            get_tidy_video_intervention(
              intervention_second = current_intervention_time,
              downsample_to_seconds = downsample_to_seconds
            ) %>%
            # Attach metadata
            dplyr::mutate(
              component_id = paste0(component, "_", current_mouse_id, "_", current_condition),
              mouse_id = current_mouse_id,
              intervention_time = current_intervention_time,
              gastric_infusion_time_end = current_gastric_infusion_time_end,
              cagedfood_intervention_time = current_cagedfood_intervention_time,
              condition = current_condition
            ) %>%
            # Attach is_peak
            dplyr::group_by(component_id) %>%
            # Get a tibble for each group
            dplyr::group_map(function(tbl, group_id) {
              # Get vector of traces
              trace_vector <-
                tbl %>%
                dplyr::pull(trace)
              # Get peaks
              peaks <- find_peaks(trace_vector, 2)
              tbl %>%
                # peaks contains the indices of the peaks, so to use them we need the row numbers
                dplyr::mutate(id = dplyr::row_number()) %>%
                dplyr::mutate(is_peak = id %in% peaks) %>%
                # Re-attach grouping variables since they are removed in group_map
                dplyr::bind_cols(group_id) %>%
                # Remove id column again since it is not needed
                dplyr::select(-id)
            }) %>%
            dplyr::bind_rows() %>%
            dplyr::ungroup()
        })
    }) %>%
    # Concatinate everything into one table
    dplyr::bind_rows() %>%
    # Fix measurement column
    dplyr::mutate(measurement =
      ifelse(second < intervention_time, "baseline",
      ifelse(second >= intervention_time & second < gastric_infusion_time_end, "gastric_infusion",
      ifelse(second >= gastric_infusion_time_end & second < cagedfood_intervention_time, "post_gastric_infusion",
      ifelse(second >= cagedfood_intervention_time, "cagedfood_intervention", "should_not_exist")))) %>%
      forcats::as_factor() %>%
      forcats::fct_relevel("baseline", "gastric_infusion", "post_gastric_infusion", "cagedfood_intervention")
    ) %>%
    # Fix measurement_length column
    dplyr::group_by(component_id) %>%
    dplyr::mutate(measurement_length =
      ifelse(second < intervention_time, intervention_time,
      ifelse(second >= intervention_time & second < gastric_infusion_time_end, gastric_infusion_time_end - intervention_time,
      ifelse(second >= gastric_infusion_time_end & second < cagedfood_intervention_time, cagedfood_intervention_time - gastric_infusion_time_end,
      ifelse(second >= cagedfood_intervention_time, max(second) - cagedfood_intervention_time, NA))))
    ) %>%
    dplyr::ungroup()
    # # Reorder conditions
    # dplyr::mutate(
    #   condition = condition %>%
    #     forcats::as_factor() %>%
    #     forcats::fct_relevel("water", "ensure", "lipids", "peptides", "glc25", "glc50", "glc125", "glc625")
    # )
}


add_gastric_infusion_vlines <- function(ggplot_obj) {
  ggplot_obj +
    ggplot2::geom_vline(
      xintercept = 600,
      linetype = "dashed",
      color = "black"
    ) +
    ggplot2::geom_vline(
      xintercept = 660,
      linetype = "dashed",
      color = "grey"
    ) +
    ggplot2::geom_vline(
      xintercept = 720,
      linetype = "dashed",
      color = "grey"
    ) +
    ggplot2::geom_vline(
      xintercept = 780,
      linetype = "dashed",
      color = "grey"
    ) +
    ggplot2::geom_vline(
      xintercept = 840,
      linetype = "dashed",
      color = "grey"
    ) +
    ggplot2::geom_vline(
      xintercept = 900,
      linetype = "dashed",
      color = "grey"
    ) +
    ggplot2::geom_vline(
      xintercept = 960,
      linetype = "dashed",
      color = "grey"
    ) +
    ggplot2::geom_vline(
      xintercept = 1020,
      linetype = "dashed",
      color = "black"
    )
}

get_pattern_percentage <- function(tbl) {
  tbl %>%
    dplyr::select(condition, pattern, component_id) %>%
    dplyr::distinct() %>%
    dplyr::group_by(condition) %>%
    dplyr::add_count(name = "condition_neuron_count") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(condition, pattern) %>%
    dplyr::add_count(name = "pattern_condition_neuron_count") %>%
    dplyr::select(condition_neuron_count, pattern_condition_neuron_count) %>%
    dplyr::distinct() %>%
    dplyr::ungroup() %>%
    dplyr::mutate(percentage = (pattern_condition_neuron_count / condition_neuron_count) * 100)
}

plot_pattern_percentage <- function(tbl) {
  tbl %>%
    ggplot2::ggplot(
      mapping = ggplot2::aes(
        x = condition,
        y = percentage,
        fill = pattern
      )
    ) +
    ggplot2::geom_bar(
      stat = "identity",
      position = "stack"
    ) +
    ggplot2::scale_fill_manual(values = c("red", "grey", "blue")) +
    ggplot2::coord_flip()
}

plot_tiled_heatmap <- function(tbl) {
  tbl %>%
    dplyr::mutate(second_ceiling = ceiling(second)) %>%
    dplyr::group_by(component_id, second_ceiling, pattern, condition) %>%
    dplyr::summarise(z_score_by_baseline = mean(z_score_by_baseline)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(second_ceiling = as.numeric(second_ceiling)) %>%
    dplyr::group_by(condition) %>%
    dplyr::group_map(function(dat, condition) {
      current_condition <- dplyr::pull(condition)
      dat %>%
        dplyr::mutate(condition = current_condition) %>%
        plot_heatmap() +
        ggplot2::facet_grid(
          pattern~.,
          scales = "free_y",
          space = "free_y"
        ) +
        # Remove z_score_by_baseline from legend
        ggplot2::theme(legend.position = "none") +
        ggplot2::ggtitle(current_condition)
    }) %>%
    # Wrap everything into one figure
    patchwork::wrap_plots()
}

plot_tiled_heatmap_gi <- function(tbl, add_vlines = FALSE) {
  tbl %>%
    dplyr::mutate(second_ceiling = ceiling(second)) %>%
    dplyr::group_by(component_id, second_ceiling, pattern, condition) %>%
    dplyr::summarise(z_score_by_baseline = mean(z_score_by_baseline)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(second_ceiling = as.numeric(second_ceiling)) %>%
    dplyr::group_by(condition) %>%
    dplyr::group_map(function(dat, condition) {
      current_condition <- dplyr::pull(condition)
      heat_plot <-
        dat %>%
        dplyr::mutate(condition = current_condition) %>%
        plot_heatmap() +
        ggplot2::facet_grid(
          pattern~.,
          scales = "free_y",
          space = "free_y"
        ) +
        # Remove z_score_by_baseline from legend
        ggplot2::theme(legend.position = "none") +
        ggplot2::ggtitle(current_condition)
      
      if (add_vlines) {
        heat_plot %>%
          add_gastric_infusion_vlines()
      } else {
        heat_plot
      }
    }) %>%
    # Wrap everything into one figure
    patchwork::wrap_plots()
}

plot_zscores_by_pattern <- function(tbl) {
  tbl %>%
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
    patchwork::wrap_plots(ncol = 1)
}

#' Limit Intervention Time
#'
#' This function filters a data frame to limit the time range around an intervention point.
#'
#' @param dat A data frame containing the data to be filtered.
#' @param limit_post_intervention_time Numeric value specifying the limit for post-intervention time.
#' @param limit_pre_intervention_time Numeric value specifying the limit for pre-intervention time.
#'
#' @return A filtered data frame based on the specified time limits.
#'
#' @import dplyr
#' @importFrom dplyr filter mutate
#'
#' @export
limit_intervention_time <- function(
  dat,
  limit_post_intervention_time,
  limit_pre_intervention_time
) {
  # Check if we have a limit
  intervention_time_set <- is.numeric(limit_post_intervention_time)
  pre_intervention_time_set <- is.numeric(limit_pre_intervention_time)
  dat_return <- dat
  # Limit post-intervention time
  if (intervention_time_set) {
    dat_return <-
      dat_return %>%
      # Filter post-intervention time
      dplyr::filter(second <= intervention_time + limit_post_intervention_time) %>%
      # Adjust measurement length
      dplyr::mutate(measurement_length = ifelse(
        measurement == "baseline",
        measurement_length,
        limit_post_intervention_time
      ))
  }
  # Limit pre-intervention time
  if (pre_intervention_time_set) {
    dat_return <-
      dat_return %>%
      # Filter pre-intervention time
      dplyr::filter(second >= intervention_time - limit_pre_intervention_time) %>%
      # Adjust measurement length
      dplyr::mutate(measurement_length = ifelse(
        measurement == "baseline",
        limit_pre_intervention_time,
        measurement_length
      ))
  }
  return(dat_return)
}

#' Calculate the standard error of a numeric vector
#'
#' This function calculates the standard error of a numeric vector by dividing
#' the standard deviation by the square root of the vector length.
#'
#' @param x A numeric vector.
#' @return The standard error of the input vector.
#' @examples
#' std_error(c(1, 2, 3, 4, 5))
#' @export
#' @importFrom stats sd
#' @importFrom base sqrt
std_error <- function(x) {
  sd(x) / sqrt(length(x))
}
