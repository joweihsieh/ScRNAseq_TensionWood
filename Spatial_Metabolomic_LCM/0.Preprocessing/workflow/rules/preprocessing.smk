rule pick_ms1_peaks_from_raw:
    input:
        raw_folder = lambda wildcards: RAW_DATA_DIR / config["RAW_PLATE_FOLDER"][wildcards.plate]
    output:
        ms1_folder = directory(INTERIM_DATA_DIR / "ms1_peak_picking/{plate}")
    log:
        "logs/" + DATETIME + "_pick_ms1_peaks_from_raw_{plate}.log"
    shell:
        """
        mkdir --parents "{output.ms1_folder}"
        ls "{input.raw_folder}" | grep ".raw$" | sed "s/^/\\/data\\//" \
            > "{output.ms1_folder}/input_file_list.txt"
        docker run \
            --rm \
            -e WINEDEBUG=-all \
            -v "{input.raw_folder}":/data:rw \
            -v "{output.ms1_folder}":/output:rw \
            proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses \
                wine msconvert \
                    --singleThreaded \
                    --filelist /output/input_file_list.txt \
                    --outdir /output \
                    --filter "peakPicking vendor" \
                    --filter "msLevel 1" \
                    --filter "scanTime [0,300]" \
                    --mzML \
        > {log} 2>&1
        """

rule average_ms1_from_mzML:
    input:
        ms1_peak_picking_folder = INTERIM_DATA_DIR / "ms1_peak_picking/{plate}"
    output:
        ms1_averaging_folder = directory(INTERIM_DATA_DIR / "ms1_averaging/{plate}")
    params:
        integration_type = "mean"
    log:
        "logs/" + DATETIME + "_average_ms1_from_mzML_{plate}.log"
    shell:
        """
        mkdir --parents "{output.ms1_averaging_folder}"
        for file in "{input.ms1_peak_picking_folder}"/*.mzML; do
            output_file=$(basename "$file")
            output_file="${{output_file%.mzML}}.tsv"
            typer tension_wood.dataset run mzml2tsv \
                --input-path "$file" \
                --output-path "{output.ms1_averaging_folder}/$output_file" \
                --integration-type {params.integration_type}
        done > {log} 2>&1
        """

rule collect_ms1_averaging_tsvs:
    input:
        ms1_averaging_folders =
            lambda wildcards: [
                INTERIM_DATA_DIR / "ms1_averaging" / x
                for x in config["PLATE_GROUP_FOLDERS"][wildcards.plates]
            ]
    output:
        ms1_tsv_folder = directory(INTERIM_DATA_DIR / "ms1_tsv/{plates}")
    shell:
        """
        mkdir --parents "{output.ms1_tsv_folder}"
        inputs="{input.ms1_averaging_folders}"
        for input in $inputs
        do
            cp $input/*.tsv {output.ms1_tsv_folder}
        done
        """

rule bin_ms1_peaks:
    input:
        ms1_tsv_folder = INTERIM_DATA_DIR / "ms1_tsv/{plates}"
    output:
        ms1_binning_folder = directory(INTERIM_DATA_DIR / "ms1_binning/{plates}")
    params:
        bin_width = 0.0003,
        bin_width_unit = "m/z",
        ms_lower_bound = 100,
        ms_upper_bound = 1500
    log:
        "logs/" + DATETIME + "_bin_ms1_peaks_{plates}.log"
    shell:
        """
        mkdir --parents "{output.ms1_binning_folder}"
        for file in "{input.ms1_tsv_folder}"/*.tsv; do
            output_file=$(basename "$file")
            typer tension_wood.dataset run bindata \
                --input-path "$file" \
                --output-path "{output.ms1_binning_folder}/$output_file" \
                --bin-width {params.bin_width} \
                --bin-width-unit {params.bin_width_unit} \
                --ms-lower-bound {params.ms_lower_bound} \
                --ms-upper-bound {params.ms_upper_bound}
        done > {log} 2>&1
        """

rule combine_ms1_binning_tsv:
    input:
        ms1_binning_folder = INTERIM_DATA_DIR / "ms1_binning/{plates}"
    output:
        ms1_combined_tsv = INTERIM_DATA_DIR / "ms1_combined_{plates}.tsv"
    log:
        "logs/" + DATETIME + "_combine_ms1_binning_tsv_{plates}.log"
    shell:
        """
        python tension_wood/dataset.py mergetsv \
            --input-tsv-folder "{input.ms1_binning_folder}" \
            --output-tsv "{output.ms1_combined_tsv}" \
        > {log} 2>&1
        # Run with python to avoid PicklingError: Can't pickle <function read_dataframe_and_mass at 0x7f6c3b631da0>: it's not the same object as tension_wood.dataset.read_dataframe_and_mass
        """