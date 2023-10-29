from utils import copy_file_cmd, unzip_file_cmd, write_json_file, load_json_file
import os
import datetime
import matplotlib.pyplot as plt
import seaborn as sns


def initialize_run(raw_data_folder:str,plate_naming_scheme:dict): 
     
    """
    Arguments:
    - raw_data_folder: path to the raw data folder in the run folder
    - plate_naming_scheme: dictionary with plate name as key and path to raw data file as value 
    
    
    Actions
    - Unzips the raw data if it is zipped 
    - Copies the raw data files to the run folder raw data folder
    
    Returns:
    Dictionary with plate raw data file name as key and path as value
    
    """

    
    #new plate dictionary with plate name as key and path as value
    
    raw_run_plates = {}
    
    for plate in plate_naming_scheme.keys():
        
        plate_name = plate
        plate_path = plate_naming_scheme[plate]
        
        # check if file exists
        if not os.path.exists(plate_path):
            raise ValueError(f"File {plate_path} does not exist")
        
        # check if file is gzipped
        if not plate_path.endswith(".gz"):
            copy_file_cmd(plate_path, raw_data_folder, f"{plate_name}.fastq")
        
        # if file is gzipped, unzip it and copy it to raw data folder
        if plate_path.endswith(".gz"):
            unzip_file_cmd(plate_path, raw_data_folder, f"{plate_name}.fastq")
        
        # add plate to run plates dictionary
        raw_run_plates[plate_name] = os.path.join(raw_data_folder,f"{plate_name}.fastq")
        
    return raw_run_plates
        
    
  
def initialize_run_json_file(run_name:str,output_directory:str,json_file_path:str,plate_naming_scheme:dict,run_plates:dict,multiprocessing_cores:int=1):
    
    """
    Arguments:
    - run_name: name of the run
    - output_directory: path to the output directory of the run
    - json_file_path: path to the json file
    - plate_naming_scheme: dictionary with plate name as key and path to raw data file as value
    - run_plates: dictionary with plate name as key and path to run-level raw data file as value
    
    Actions:
    - Creates a dictionary with the run information and the plates information
    - Writes the dictionary to a json file
    
    Returns:
    None
    """
    
    
    run_dictionary = { }
    
    run_dictionary["run_info"]= {
    
    "run_name" : run_name,
    "run_date" : datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    "output_directory" : output_directory,
    "json_file_path" : json_file_path,
    "plates" : [plate for plate in plate_naming_scheme.keys()],
    "raw_plate_fastq_paths" : [plate_naming_scheme[plate] for plate in plate_naming_scheme.keys()],
    "multiprocessing_cores" : multiprocessing_cores,
    "barcodes_csv": "path/to/barcodes/csv/file",
    "freebarcodes_search_error": 1,
    "freebarcodes_file": "path/to/freebarcodes/file",
    'alignment_algorithm' : 'algorithm',
    'consensus_threshold':'',
    'consensus_ambiguous_base':'',
    'run_raw_consensus_csv_file_path' : 'path/to/raw/consensus/csv/file',
    'run_trimmed_reconstructed_consensus_csv_file_path' : 'path/to/reconstructed/consensus/csv/file',
    'visualization_folder_output_path' : 'path/to/visualization/output/file',
    
    }
    
    for plate in plate_naming_scheme.keys():
        
        run_dictionary[plate]={
    
        'fastq_source_file' : plate_naming_scheme[plate],
        'fastq_destination_file' : run_plates[plate],
        'raw_file_stats' : 
            {
            'num_sequences': 000000, 
            'mean_length': 00000000, 
            'median_length': 000000, 
            'min_length': 000000, 
            'max_length': 00000, 
            'length_std': 0000000,
            },
        'length_filtering': (10,100),
        'NGS_minimum_quality_score' : 10,
        'post_fastp_stats' :
            {
            'num_sequences': 000000, 
            'mean_length': 00000000, 
            'median_length': 000000, 
            'min_length': 000000, 
            'max_length': 00000, 
            'length_std': 0000000,
            },
            
        'fastp_folder_output_path' : 'path/to/fastp/output/file',
        'cnst_fwd_read' : '',
        'freebarcodes_decoded_file_path' : 'path/to/freebarcodes/output/file',
        'demultiplexed_folder_output_path' : 'path/to/demultiplexed/fastq/file',
        'alignment_folder_output_path' : 'path/to/alignment/output/file',
        'raw_consensus_csv_file_path' : 'path/to/raw/consensus/csv/file',
        'five_prime_trim_seq' : '',
        'three_prime_trim_seq' : '',
        'five_prime_reconstruct_seq' : '',
        'three_prime_reconstruct_seq' : '',
        'plate_raw_consensus_csv_file_path' : 'path/to/reconstructed/consensus/csv/file',
        'plate_trimmed_reconstructed_consensus_csv_file_path' : 'path/to/reconstructed/consensus/csv/file',
        }
    
    # write to json file
    write_json_file(run_dictionary,json_file_path)
    
    return None
    
  
    
def fastq_sequence_length_histogram(plate,lengths, stats_dict, output_directory_abs_path:str): 
    
    """
    Arguments:
    - plate: name of the plate
    - lengths: list of lengths of sequences
    - stats_dict: dictionary with stats of the plate
    - output_directory_abs_path: absolute path to the output directory
    
    Actions:
    - Plots histogram of sequence lengths
    - Prints stats dictionary of the plate on the histogram
    - Saves histogram to output directory
    - Prints stats dictionary of the plate to console
    
    Returns:
    None
    """

    
    #get the legths of sequences, along with some statistics
    num_sequences = stats_dict["num_sequences"]
    mean_length = stats_dict["mean_length"]
    median_length = stats_dict["median_length"]
    min_length = stats_dict["min_length"]
    max_length = stats_dict["max_length"]
    std_length = stats_dict["length_std"]
    
    #plot histogram
    sns.set(rc={"figure.dpi": 300, "savefig.dpi": 300})
    sns.set_style("white")
    sns.set_context("paper")
    sns.set_palette("colorblind")
    sns.histplot(lengths, bins=100, kde=False)
    plt.xlabel("Sequence length (bp)")
    plt.ylabel("Count")
    plt.title(f"Sequence length histogram for {plate}")
    
    #add legend with number of sequences, mean, median, min, max std
    plt_text = f"number of sequences: {num_sequences}\nmean len: {mean_length:.0f}\nmedian len: {median_length:.0f}\nmin len: {min_length}\nmax len: {max_length}\nlen std: {std_length:.0f}"
    plt.text(0.95, 0.95, plt_text, verticalalignment='top', horizontalalignment='right', transform=plt.gca().transAxes, fontsize=12)
    
    plt.savefig(os.path.join(output_directory_abs_path,f"{plate}_sequence_length_histogram.png"))
    plt.show()
    print(f"Sequence length histogram for {plate} saved to {output_directory_abs_path}")
    
    return None
    
    
  

def raw_stats_update_json_file(json_file_path:str,plate:str,plate_raw_stats_dictionary:dict):
    """
    Arguments:
    - json_file_path: path to the json file
    - plate: name of the plate
    - plate_raw_stats_dictionary: dictionary with stats of the plate
    
    Actions:
    - Updates the json file with the stats of the plate
    
    Returns:
    None
    """
    
    
        
    run_dictionary=load_json_file(json_file_path)
    
    
        
    run_dictionary[plate]['raw_file_stats']['num_sequences'] = plate_raw_stats_dictionary['num_sequences']
    run_dictionary[plate]['raw_file_stats']['mean_length'] = plate_raw_stats_dictionary['mean_length']
    run_dictionary[plate]['raw_file_stats']['median_length'] = plate_raw_stats_dictionary['median_length']
    run_dictionary[plate]['raw_file_stats']['min_length'] = plate_raw_stats_dictionary['min_length']
    run_dictionary[plate]['raw_file_stats']['max_length'] = plate_raw_stats_dictionary['max_length']
    run_dictionary[plate]['raw_file_stats']['length_std'] = plate_raw_stats_dictionary['length_std']
    
    
    write_json_file(run_dictionary,json_file_path)
    return None
    