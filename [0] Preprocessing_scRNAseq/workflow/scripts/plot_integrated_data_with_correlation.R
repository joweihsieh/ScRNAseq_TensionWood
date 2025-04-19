suppressPackageStartupMessages({
    library(magrittr)
    library(dplyr)
    library(tools)
    library(RColorBrewer)
    ori_par <- par(no.readonly = TRUE)
})


# Define functions =============================================================
#' Output figure
output_png_figure <- function(
    plotting_function,
    # x, y, col, main = "",
    output_figure = FALSE,
    output_path = "temp.png",
    output_without_margin = FALSE,
    ...
) {
    if (output_figure) {
        png(output_path,
            pointsize = 10, res = 300,
            width = 20, height = 15, units = "cm")
    }

    if (output_without_margin) {
        par(mai = c(0, 0, 0, 0))
    } else {
        par(mai = ori_par$mai)
    }

    plotting_function(
        output_figure = output_figure,
        output_path = output_path,
        output_without_margin = output_without_margin,
        ...
    )

    par(mai = ori_par$mai)

    if (output_figure) dev.off()
}

#' Plot seurat UMAP colored by correlation
plot_with_correlation <- function(
    MS_plotting_df = merge_plotting_df,
    cor_colname,
    output_without_margin,
    xlim, ylim,
    sorted_order,
    ...
) {
    x <- MS_plotting_df$UMAP.1
    y <- MS_plotting_df$UMAP.2
    cor_vector <- MS_plotting_df[, cor_colname] # "Correlation_PtrFiber1_1"

    upper_bound <- quantile(cor_vector, 0.95)
    lower_bound <- upper_bound * 0.8 # quantile(cor_vector, 0.66)
    lower <- 0.25
    minivalue <- min(cor_vector)
    #maxivalue <- max(cor_vector)
    maxivalue <- 0.55

    for (i in 1:length(cor_vector)){
        if (cor_vector[i] < lower){
            cor_vector[i] = 0
        }
    }


    color_index <-
        #((cor_vector - minivalue) / (maxivalue - minivalue)) %>%
        ((cor_vector - lower) / (maxivalue - lower)) %>%
        pmax(0) %>%
        pmin(1) %>%
        multiply_by(500) %>%
        round(digits = 0) %>%
        add(1)


    # color_tick <- rev(brewer.pal(11, "Spectral"))
    # color_pool <- colorRampPalette(color_tick[7:11])(501)
    color_tick <- c("#EEF2F9", "#C44233")
    color_pool <- colorRampPalette(color_tick)(501)

    if (sorted_order) {
        plot_order <- order(cor_vector)
    } else {
        plot_order <- sample(length(x))
    }
    plot(x[plot_order], y[plot_order],
        col = color_pool[color_index[plot_order]],
        pch = 20, cex = 0.5,
        xlab = "UMAP_1", ylab = "UMAP_2",
        main = ifelse(
            output_without_margin, "",
            paste0(
                "Min:", round(minivalue, 2), " ",
                #"Min:", round(min(cor_vector), 2), " ",
                "LB:", round(lower_bound, 2), " ",
                "UB:", round(upper_bound, 2), " ",
                "Max:", round(max(cor_vector), 2)
            )
        ),
        xlim = xlim, ylim = ylim,
        axes = !output_without_margin, las = 1
    )

    par(xpd = TRUE)
    current_xlim <- par()$usr[1:2]
    current_ylim <- par()$usr[3:4]
    legend_x_range <- c(sum(current_xlim * c(0.25, 0.75)), current_xlim[2])
    legend_x_vector <-
        seq(
            from = legend_x_range[1],
            to = legend_x_range[2],
            length.out = length(color_pool) + 1
        )
    legend_y_range <-
        c(sum(current_ylim * c(-0.1, 1.1)), sum(current_ylim * c(-0.15, 1.15)))
    for (i in seq_along(color_pool)) {
        lines(
            legend_x_vector[c(i, i + 1)],
            rep(legend_y_range[1], 2),
            col = color_pool[i],
            lwd = 20, lend = "square"
        )
    }
    text(
        legend_x_range[1], legend_y_range[2],
        #sprintf("%.2f", lower_bound), cex = 1
        sprintf("%.2f", lower), cex = 1

    )
    text(
        legend_x_range[2], legend_y_range[2],
        #sprintf("%.2f", upper_bound), cex = 1
        sprintf("%.2f", maxivalue), cex = 1

    )
    par(xpd = FALSE)
}


# Set parameters ===============================================================
#' Get input parameters from command line
input_MS_plotting_csv <- snakemake@input$MS_plotting_csv
input_SS_MQD_plotting_csv <- snakemake@input$SS_MQD_plotting_csv
SC_sample <- snakemake@wildcards$SC_sample
output_figure_folder <- snakemake@output$figure_folder


# input_MS_plotting_csv <-
#     file.path(
#         "results/Multi_species_analysis/all_plotting_tables_addSS",
#         "plotting_PalTenXBatch2_seed_42_md_0.3_nn_30.csv"
#     )
# input_SS_MQD_plotting_csv <-
#     file.path(
#         "results/Single_species_analysis/all_plotting_tables_cor_MQD",
#         "plotting_{SC_sample}_cor_{MQD_group}.csv"
#     )
# output_figure_folder <-
#     file.path(
#         "results/Multi_species_analysis/UMAP_with_correlation",
#         "integration_PalTenXBatch2_seed_42_md_0.3_nn_30",
#         "{SC_sample}_cor_{MQD_group}"
#     )
# SC_sample <- c("TenX_PalBarkChen2021", "TenX_PalWoodChen2021")
# MQD_group <- c("PtrCambium1", "PtrFiber1", "PtrVessel1", "PtrRay1")


# Implementation ===============================================================
#' Create output directory
if (!dir.exists(output_figure_folder)) {
    dir.create(output_figure_folder, recursive = TRUE)
}

#' Input plotting information
MS_plotting_df <- read.csv(input_MS_plotting_csv)
SS_MQD_plotting_df <- read.csv(input_SS_MQD_plotting_csv)
is_barcode_column <- grepl("Barcode", colnames(SS_MQD_plotting_df))
is_correlation_column <- grepl("Correlation_", colnames(SS_MQD_plotting_df))

#' Merge correlation into plotting data frame
selected_MS_plotting_df <- filter(MS_plotting_df, Sample == SC_sample)
if (nrow(selected_MS_plotting_df) == 0) {
    stop("No barcodes in SC_sample are in MS_plotting_csv.")
}
merge_plotting_df <- full_join(
    selected_MS_plotting_df,
    SS_MQD_plotting_df[, is_barcode_column | is_correlation_column],
    by = "Barcode"
)

#' Get xlim and ylim
plot(MS_plotting_df$UMAP.1, MS_plotting_df$UMAP.2)
plot_xlim <- par()$usr[1:2]
plot_ylim <- par()$usr[3:4]

#' Output figures
correlation_columns <- colnames(SS_MQD_plotting_df)[is_correlation_column]
for (selected_column in correlation_columns) {
    output_png_figure(
        plotting_function = plot_with_correlation,
        MS_plotting_df = merge_plotting_df,
        cor_colname = selected_column,
        output_figure = TRUE,
        output_path =
            paste0(output_figure_folder, "/", selected_column, ".png"),
        output_without_margin = FALSE,
        xlim = plot_xlim, ylim = plot_ylim,
        sorted_order = FALSE
    )
    output_png_figure(
        plotting_function = plot_with_correlation,
        MS_plotting_df = merge_plotting_df,
        cor_colname = selected_column,
        output_figure = TRUE,
        output_path =
            paste0(output_figure_folder, "/", selected_column, "_Clear.png"),
        output_without_margin = TRUE,
        xlim = plot_xlim, ylim = plot_ylim,
        sorted_order = FALSE
    )
    output_png_figure(
        plotting_function = plot_with_correlation,
        MS_plotting_df = merge_plotting_df,
        cor_colname = selected_column,
        output_figure = TRUE,
        output_path =
            paste0(output_figure_folder, "/", selected_column, "_Sorted.png"),
        output_without_margin = FALSE,
        xlim = plot_xlim, ylim = plot_ylim,
        sorted_order = TRUE
    )
    output_png_figure(
        plotting_function = plot_with_correlation,
        MS_plotting_df = merge_plotting_df,
        cor_colname = selected_column,
        output_figure = TRUE,
        output_path =
            paste0(output_figure_folder, "/", selected_column, "_Sorted_Clear.png"),
        output_without_margin = TRUE,
        xlim = plot_xlim, ylim = plot_ylim,
        sorted_order = TRUE
    )
}
