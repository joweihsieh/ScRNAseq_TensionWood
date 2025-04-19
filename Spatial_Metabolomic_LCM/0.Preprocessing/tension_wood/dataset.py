# %%
import numpy as np
import pandas as pd
import pyopenms as oms
import multiprocessing as mp
import itertools
# from functools import partial

from pathlib import Path
import typer
from loguru import logger
# from tqdm import tqdm

# %%
app = typer.Typer(add_completion=False)

@app.command()
def mzml2tsv(input_path = Path('.'), output_path = Path('.'), integration_type = 'mean'):
    logger.info(f'Read file: {input_path}')
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)
    spectra = exp.getSpectra()
    spectra_ms1 = [s for s in spectra if s.getMSLevel() == 1]

    positions = []
    intensities = []
    polarities = []
    polarity_map = {
        oms.IonSource.Polarity.POSITIVE: 'positive',
        oms.IonSource.Polarity.NEGATIVE: 'negative',
        oms.IonSource.Polarity.POLNULL: 'polnull'
    }
    scan_number = len(spectra_ms1)
    for i in range(scan_number):
        positions.extend(spectra_ms1[i].get_peaks()[0])
        intensities.extend(spectra_ms1[i].get_peaks()[1])
        polarities.extend(
            [polarity_map.get(spectra_ms1[i].getInstrumentSettings().getPolarity(), None)] *
            len(spectra_ms1[i].get_peaks()[0])
        )
    data = pd.DataFrame(
        {
            'position': positions,
            'intensity': intensities,
            'mode': polarities
        }
    )
    data['position'] = np.round(data['position'], 4)
    if integration_type == 'mean':
        data['intensity'] = data['intensity'] / scan_number
    data = data.groupby(['mode', 'position']).agg({'intensity': 'sum'}).reset_index()

    logger.info(f'Write file: {output_path}')
    data.to_csv(output_path, sep='\t', index=False, header=False)

def create_bins(lower_bound, upper_bound, bin_width):
    bins = []
    added_bound = lower_bound
    while added_bound < upper_bound:
        bins.append(added_bound)
        added_bound = added_bound + added_bound * bin_width / 1e6
    bins.append(added_bound)
    return bins

def get_bins_center(bins):
    n = len(bins)
    bins = np.array(bins)
    bin_center = (bins[:n-1] + bins[1:]) / 2
    return bin_center

@app.command()
def bindata(
    input_path = Path('.'),
    output_path = Path('.'),
    BIN_WIDTH: float = 0.0003,
    BIN_WIDTH_UNIT: str = 'm/z',
    MS_LOWER_BOUND: float = 100,
    MS_UPPER_BOUND: float = 1500
):
    logger.info(f'Read file: {input_path}')
    data = pd.read_table(input_path, header=None, names=['mode', 'position', 'intensity'])
    data = data[(data['position'] >= MS_LOWER_BOUND) & (data['position'] <= MS_UPPER_BOUND)]

    if BIN_WIDTH_UNIT == 'ppm':
        bins = create_bins(MS_LOWER_BOUND, MS_UPPER_BOUND, BIN_WIDTH)
    elif BIN_WIDTH_UNIT == 'm/z':
        bins = np.arange(MS_LOWER_BOUND - BIN_WIDTH / 2, MS_UPPER_BOUND + BIN_WIDTH, BIN_WIDTH)
    bins_center = get_bins_center(bins)
    data['bin_index'] = np.digitize(data['position'], bins, right=False) - 1

    binned_data = data.groupby(['mode', 'bin_index']).agg({'position': 'mean', 'intensity': 'sum'}).reset_index()
    bin_center = bins_center[binned_data['bin_index']]
    binned_data['bin_center'] = [f'{i:.4f}' for i in bin_center]

    output = binned_data[['mode', 'bin_center', 'intensity']]
    logger.info(f'Write file: {output_path}')
    output.to_csv(output_path, sep='\t', index=False, header=False)


def read_tsv_file(file_path):
    """
    Reads a tsv file and renames the columns for merging.
    """
    try:
        df = pd.read_table(
            file_path,
            sep='\t',
            header=None,
            names=['mode', 'mass_binned', f'intensity_sum.{file_path.stem}']
        )
        return df
    except FileNotFoundError:
        logger.error(f'File not found: {file_path}')
        return None
    except pd.errors.EmptyDataError:
        logger.error(f'Empty tsv file: {file_path}')
        return None
    except Exception as e:
        logger.error(f'Error reading file {file_path}: {e}')
        return None

def concat_dataframes(base_df, df_chunk):
    """
    Function to join a list of dataframes with the base_df.
    """
    df_chunk = [df.set_index(['mode', 'mass_binned']) for df in df_chunk]
    df_chunk = [base_df.join(df, how='left') for df in df_chunk]
    result = pd.concat(df_chunk, axis=1)
    return result

def combine_dataframes(dataframes, mass_set, N_PROCESSORS=16):
    """
    Combines a list of DataFrames on the 'mass_binned' column using outer merge.
    """
    if not dataframes:
        logger.error('No valid dataframes to combine.')
        return None

    try:
        mass_list = sorted(list(mass_set))
        index = pd.MultiIndex.from_tuples(mass_list, names=['mode', 'mass_binned'])
        logger.info(f'Created MultiIndex with {len(index)} (mode, mass_binned) combinations.')
        base_df = pd.DataFrame(index=index)

        logger.info(f'Concat {len(dataframes)} dataframes.')
        chunk_size = len(dataframes) // N_PROCESSORS
        chunks = [dataframes[i:i + chunk_size] for i in range(0, len(dataframes), chunk_size)]
        del dataframes
        logger.info(f'Concat {len(chunks)} dataframe chunks.')
        with mp.Pool(processes=N_PROCESSORS) as pool:
            results = pool.starmap(concat_dataframes, zip(itertools.repeat(base_df, len(chunks)), chunks))
        del chunks
        logger.info(f'Concat {len(results)} combination results.')
        final_result = pd.concat(results, axis=1)
        final_result.reset_index(inplace=True)
        del results
        return final_result
    except Exception as e:
        logger.error(f'Error combining DataFrames: {e}')
        return None

def read_dataframe_and_mass(tsv_path):
    df = read_tsv_file(tsv_path)
    if df is not None:
        return(df, zip(df['mode'], df['mass_binned']))

def combine_tsvs(input_tsvs, N_PROCESSORS=16):
    """
    Combines multiple tsv files into a single DataFrame based on 'mode' and 'mass_binned'.
    """
    logger.info('Input tsv files')
    with mp.Pool(processes=N_PROCESSORS) as pool:
        results = pool.map(read_dataframe_and_mass, input_tsvs)
    dataframes = []
    mass_set = set()
    for content in results:
        dataframes.append(content[0])
        mass_set.update(content[1])

    combined_df = combine_dataframes(dataframes, mass_set)
    if combined_df is not None:
        logger.info(f'Successfully combined {len(dataframes)} tsv files.')
    return combined_df

@app.command()
def mergetsv(input_tsv_folder = Path('.'), output_tsv = Path('.')):
    logger.info(f'Processing {input_tsv_folder}')
    input_tsvs = sorted(Path(input_tsv_folder).glob('*.tsv'))

    combined_df = combine_tsvs(input_tsvs)

    logger.info(f'Output to {output_tsv}')
    combined_df.to_csv(output_tsv, sep='\t', index=False, header=True)

if __name__ == '__main__':
    app()
