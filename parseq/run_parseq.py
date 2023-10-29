
from utils import create_folder, fastq_file_read_stats, load_json_file
from initialize import initialize_run, initialize_run_json_file, fastq_sequence_length_histogram, raw_stats_update_json_file
from raw_file_processing import pre_fastp_update_json_file, run_fastp_on_plates, post_fastp_stats_update_json_file
from demultiplex import pre_demultiplex_json_update, create_freebarcodes_file, freebarcodes_decode, pooled_separate_freebarcodes_decoded_file_into_wells
from align import align_plates
from consensus import get_consensus, pre_trimming_reconstruction_json_update, trim_and_reconstruct_consesnsus
from visualizations import create_visualizations_raw_consensus, create_visualizations_trimmed_reconstructed_consensus
import os.path as path



def set_up_run (run_name:str,output_directory:str,plate_naming_scheme:dict,multiprocessing_cores:int=1):

    """
    Sets up the run folder that contains the raw data and the output data
    Unzips the raw data if it is zipped
    Copies the raw data to the run folder
    Gets stats about the raw fastq files
    creates and updates json with relevant information
    """
    
    #create run folder
    run_folder = path.join(output_directory, run_name)
    create_folder(run_folder)
    
    #create raw data folder   
    raw_data_folder = path.join(run_folder,"raw_data")
    create_folder(raw_data_folder)
    
    #initialize run and get run plates
    run_plates = initialize_run(raw_data_folder, plate_naming_scheme)
    
    
    # create run json file
    json_file_path = path.join(output_directory,run_name,run_name+".json")
    initialize_run_json_file(run_name,output_directory, json_file_path, plate_naming_scheme, run_plates,multiprocessing_cores)
   
    
    
    # create lengths histograms and get stats for each plate
    fastq_length_histograms_path = os.path.join(output_directory, run_name, "fastq_length_histograms")
    create_folder(fastq_length_histograms_path)
    
    for item in run_plates.items():
        
        plate = item[0]
        
        fastq_path = run_plates[plate]
        
        #Generate fastq length histogram
        raw_plate_lengths, raw_plate_stats_dictionary = fastq_file_read_stats(fastq_path)
        fastq_sequence_length_histogram(plate, raw_plate_lengths, raw_plate_stats_dictionary, fastq_length_histograms_path)
        
        # update json file
        raw_stats_update_json_file(json_file_path, plate, raw_plate_stats_dictionary)
    
        # print stats dictionary for each plate
        print(f"Plate {plate} length dictionary: \n{raw_plate_stats_dictionary}")
    

    return json_file_path
    
    
    
    
def fastp_process(run_json_file_path:str, length_filtering_dict:dict,ngs_read_quality_filtering_dict:dict):
    
    run_dictionary = load_json_file(run_json_file_path)
    
    output_directory = run_dictionary["run_info"]["output_directory"]
    run_name = run_dictionary["run_info"]["run_name"]
    
    fastp_output_path = path.join(output_directory,run_name,"fastp_output")
    create_folder(fastp_output_path)
    
    pre_fastp_update_json_file(run_json_file_path,fastp_output_path,length_filtering_dict,ngs_read_quality_filtering_dict)
    
    run_fastp_on_plates(run_json_file_path)
    post_fastp_stats_update_json_file(run_json_file_path)
    
    
    return None



def demultiplex_plates(run_json_file_path:str,barcodes_csv_path:str, barcode_search_error:int,fwd_read_constant_sequence_dict:dict):
    
    run_dictionary = load_json_file(run_json_file_path) 
    output_directory = run_dictionary["run_info"]["output_directory"]
    run_name = run_dictionary["run_info"]["run_name"]  
    
      
    # update json file
    pre_demultiplex_json_update(run_json_file_path, barcodes_csv_path, barcode_search_error, fwd_read_constant_sequence_dict)
    
    
    # create free barcodes output folder, free barcodes file and run freebarcodes decoding
    # free barcodes file  is a text file with a list of all the barcodes formatted for freebarcodes
    freebarcodes_output_path = path.join(output_directory,run_name,"freebarcodes_output")
    create_folder(freebarcodes_output_path)
    create_freebarcodes_file(run_json_file_path,freebarcodes_output_path)
    freebarcodes_decode(run_json_file_path, freebarcodes_output_path)
    
    
    # create demultiplexed output folder and run demultiplexing on each plate decoding output
    demultiplexed_files_path = path.join(output_directory,run_name,"demultiplexed")
    create_folder(demultiplexed_files_path)
    pooled_separate_freebarcodes_decoded_file_into_wells(run_json_file_path, demultiplexed_files_path)
    
    return None


def align_wells_and_get_consensus (run_json_file_path:str, alignment_algorithm:str, consensus_threshold:float,ambiguous_base:str):
    
    
    run_dictionary = load_json_file(run_json_file_path)
    output_directory = run_dictionary["run_info"]["output_directory"]
    run_name = run_dictionary["run_info"]["run_name"]
    
    #create alignment folder
    alignment_folder_path = path.join(output_directory,run_name,"aligned")
    create_folder(alignment_folder_path)
    #align plates
    align_plates(run_json_file_path,alignment_folder_path,alignment_algorithm)
    
    #create consensus folder
    consensus_folder_path = path.join(output_directory,run_name,"consensus")
    create_folder(consensus_folder_path)
    #get consensus
    get_consensus(run_json_file_path,consensus_folder_path,consensus_threshold,ambiguous_base)
    
    return None
    
    
        
    




def trim_and_reconstruct(run_json_file_path:str,trim_dictionary:dict,reconstruct_dictionary:dict):
    
   # load json
    run_dictionary = load_json_file(run_json_file_path)
    output_directory = run_dictionary["run_info"]["output_directory"]
    run_name = run_dictionary["run_info"]["run_name"]
    
    #path to output folder (same comsensus folder as align_wells_and_get_consensus)
    consensus_folder_path = path.join(output_directory,run_name,"consensus")
    
    pre_trimming_reconstruction_json_update(run_json_file_path,trim_dictionary,reconstruct_dictionary)
    trim_and_reconstruct_consesnsus(run_json_file_path,consensus_folder_path)
    
    return None
    
    
    

def create_visualizations(run_json_file_path:str, raw_consensus_visualizations:bool,trimmed_reconstructed_visualizations:bool):
    
    run_dictionary = load_json_file(run_json_file_path)
    output_directory = run_dictionary["run_info"]["output_directory"]
    run_name = run_dictionary["run_info"]["run_name"]
    
    visualizations_directory = path.join(output_directory,run_name,"visualizations")
    create_folder(visualizations_directory)
    
    if raw_consensus_visualizations:
        create_visualizations_raw_consensus(run_json_file_path, visualizations_directory)
    
    if trimmed_reconstructed_visualizations:
        create_visualizations_trimmed_reconstructed_consensus(run_json_file_path, visualizations_directory)    
    
    return None
      