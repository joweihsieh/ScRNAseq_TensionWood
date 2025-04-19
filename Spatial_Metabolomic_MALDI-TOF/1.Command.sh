
############################### Bio1
library(readxl)

input_dir <- "20241205_stem+CHCA_image"
output_dir <- paste0(input_dir,"/output")
script_path <- "2.Spatial_Heatmap_ExtractValue_Final.R"
excel_file <- "20241205_stem+CHCA_image.xlsx"

setwd(input_dir)

data <- read_excel(excel_file)


png_files <- list.files(
  input_dir, 
  pattern = "^(Tension|Normal).*\\.png$",
  full.names = FALSE
)

generate_command <- function(png_file, x_grids, y_grids) {
  id <- sub("\\.png$", "", png_file)
  paste(
    "Rscript", script_path,
    paste0('"', input_dir, '"'),
    paste0('"', output_dir, '"'),
    paste0('"', png_file, '"'),
    x_grids, y_grids, paste0('"', id, '"')
  )
}

commands <- character()

for (i in seq_along(png_files)) {
  file_names_without_png <- sub("\\.png$", "", png_files[i])

  x_grids <- data[data$file_name == file_names_without_png, "X_grids"]
  y_grids <- data[data$file_name == file_names_without_png, "Y_grids"]

  if (length(x_grids) > 0 && length(y_grids) > 0) {
    commands <- c(commands, generate_command(png_files[i], x_grids, y_grids))
  } else {
    cat("Warning:", file_names_without_png, "not found\n")
  }
}

output_file <- "generated_commands.sh"
writeLines(commands, con = output_file)

cat(commands, sep = "\n")

############################### Bio2
library(readxl)

input_dir <- "20241205_stem+CHCA_image_2"
output_dir <- paste0(input_dir,"/output")
script_path <- "2.Spatial_Heatmap_ExtractValue_Final_Rotate.R"
excel_file <- "20241205_stem+CHCA_image.xlsx"

setwd(input_dir)

data <- read_excel(excel_file)



png_files <- list.files(
  input_dir, 
  pattern = "^(Tension|Normal).*\\.png$", 
  full.names = FALSE
)

generate_command <- function(png_file, x_grids, y_grids) {
  id <- sub("\\.png$", "", png_file)
  paste(
    "Rscript", script_path,
    paste0('"', input_dir, '"'),
    paste0('"', output_dir, '"'),
    paste0('"', png_file, '"'),
    x_grids, y_grids, paste0('"', id, '"')
  )
}

commands <- character()

for (i in seq_along(png_files)) {
  file_names_without_png <- sub("\\.png$", "", png_files[i])

  x_grids <- data[data$file_name == file_names_without_png, "X_grids"]
  y_grids <- data[data$file_name == file_names_without_png, "Y_grids"]

  if (length(x_grids) > 0 && length(y_grids) > 0) {
    commands <- c(commands, generate_command(png_files[i], x_grids, y_grids))
  } else {
    cat("Warning:", file_names_without_png, "not found\n")
  }
}

output_file <- "generated_commands.sh"
writeLines(commands, con = output_file)

cat(commands, sep = "\n")


############################### Bio3
library(readxl)

input_dir <- "20241205_stem+CHCA_image_3"
output_dir <- paste0(input_dir,"/output")
script_path <- "2.Spatial_Heatmap_ExtractValue_Final_Rotate_2.R"
excel_file <- "20241205_stem+CHCA_image.xlsx"

setwd(input_dir)

data <- read_excel(excel_file)



png_files <- list.files(
  input_dir, 
  pattern = "^(Tension|Normal).*\\.png$", 
  full.names = FALSE
)

generate_command <- function(png_file, x_grids, y_grids) {
  id <- sub("\\.png$", "", png_file)
  paste(
    "Rscript", script_path,
    paste0('"', input_dir, '"'),
    paste0('"', output_dir, '"'),
    paste0('"', png_file, '"'),
    x_grids, y_grids, paste0('"', id, '"')
  )
}

commands <- character()

for (i in seq_along(png_files)) {
  file_names_without_png <- sub("\\.png$", "", png_files[i])

  x_grids <- data[data$file_name == file_names_without_png, "X_grids"]
  y_grids <- data[data$file_name == file_names_without_png, "Y_grids"]

  if (length(x_grids) > 0 && length(y_grids) > 0) {
    commands <- c(commands, generate_command(png_files[i], x_grids, y_grids))
  } else {
    cat("Warning:", file_names_without_png, "not found\n")
  }
}

output_file <- "generated_commands.sh"
writeLines(commands, con = output_file)

cat(commands, sep = "\n")




############################### run command line in bash 
bash 20241205_stem+CHCA_image/generated_commands.sh
bash 20241205_stem+CHCA_image_2/generated_commands.sh
bash 20241205_stem+CHCA_image_3/generated_commands.sh




