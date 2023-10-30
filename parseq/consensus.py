from .utils import load_json_file, write_json_file, create_folder, get_read_count_from_fasta
from Bio import AlignIO
from Bio.Align import AlignInfo
import pandas as pd
import os
from multiprocessing import Pool
import edlib





def get_consensus(run_json_file_path:str,consensus_folder_path:str,consensus_threshold:float,ambiguous_base:str):
    
    """
    Arguments:
    - run_json_file_path: path to the json file
    - consensus_folder_path: path to the consensus output folder
    - consensus_threshold: threshold for consensus
    - ambiguous_base: ambiguous base for consensus
    
    Actions:
    - Creates consensus csv file for each plate
    - Creates run level consensus csv file
    - updates json file with the consensus folder path and consensus threshold
    - Runs multiprocessing on the plate level
    
    Returns:
    None
    
    """
    
    # load json file
    run_dictionary = load_json_file(run_json_file_path)
    
    # load variables from json file
    multiprocessing_cores = run_dictionary["run_info"]["multiprocessing_cores"]
    run_name = run_dictionary["run_info"]["run_name"]
    output_directory = run_dictionary["run_info"]["output_directory"]
    
    #create raw consensus folder
    raw_consensus_folder_path = os.path.join(consensus_folder_path,'raw_consensus')
    create_folder(raw_consensus_folder_path)
    
    #update dictionary with consensus_threshold and ambiguous_base
    
    run_dictionary["run_info"]["consensus_threshold"] = consensus_threshold
    run_dictionary["run_info"]["consensus_ambiguous_base"] = ambiguous_base
    
    
    # get variables for multiprocessing    
    plates = []
    input_folders = []
    output_files = []
    consensus_thresholds = []
    ambiguous_bases = []
    
    
    for item in run_dictionary.items(): # loop over plates to get and set variables from the json file
        
        if item[0] == 'run_info':
            pass
        else:
            plate = item[0]
            input_folder = run_dictionary[plate]['alignment_folder_output_path']
            output_file = os.path.join(raw_consensus_folder_path,plate+'_raw_consesus.csv')
            run_dictionary[plate]['plate_raw_consensus_csv_file_path'] = output_file
            
            plates.append(plate)
            input_folders.append(input_folder)
            output_files.append(output_file)
            consensus_thresholds.append(consensus_threshold)
            ambiguous_bases.append(ambiguous_base)

    # update json file
    write_json_file(run_dictionary,run_json_file_path)
    
    # zip input to get_raw_plate_consensus for multiprocessing
    zipped_input = zip(plates, input_folders, output_files, consensus_thresholds, ambiguous_bases)
    
    print("Creating raw consensus csv files for each plate")
    
    # run multiprocessing
    with Pool(multiprocessing_cores) as pool:
        pool.starmap(get_raw_plate_consensus, zipped_input)
    
    
    # create run level raw consensus csv file
    raw_run_df = pd.DataFrame(columns=['plate','well','raw_consensus','raw_consensus_length','raw_amb_char_count', 'read_count'])
    
    # loop over plates to concat the raw consensus csv files
    for plate in plates:
        plate_df = pd.read_csv(run_dictionary[plate]['plate_raw_consensus_csv_file_path'])
        raw_run_df = pd.concat([raw_run_df,plate_df])
    
    # save run level raw consensus csv file
    run_consensus_csv_file_path = os.path.join(output_directory,run_name,run_name+'_raw_consesus.csv')
    raw_run_df.to_csv(run_consensus_csv_file_path,index=False)
    run_dictionary["run_info"]["run_raw_consensus_csv_file_path"] = run_consensus_csv_file_path
    # update json file
    write_json_file(run_dictionary,run_json_file_path)
    
    return None
    
    
             
            

def get_raw_plate_consensus(plate:str,input_folder:str,output_file:str,threshold:float=0.6,ambiguous_base:str='N'):
    
    """
    Arguments:
    - plate: plate name
    - input_folder: path to the input folder with the aligned well files for each plate
    - output_file: path to the output file for the raw consensus csv file for each plate
    - threshold: threshold for consensus
    - ambiguous_base: ambiguous base for consensus
    
    Actions:
    - Creates raw consensus csv file for each plate
    
    Returns:
    None
    """
    
    
    raw_consensus_df = pd.DataFrame(columns=['plate','well','raw_consensus','raw_consensus_length','raw_amb_char_count', 'read_count'])
    
    
    
    wells=[]
    raw_consensus_seqs=[]
    raw_consensus_lengths=[]
    raw_amb_char_counts=[]
    read_counts=[]
    
    
    # list all fasta files in the input folder
    fasta_files = [file for file in os.listdir(input_folder) if file.endswith(".fasta")]
    
    
    # loop over fasta files
    for fasta_file in fasta_files:
        
        well = fasta_file.split('.')[0]
        wells.append(well)
        
        well_path = os.path.join(input_folder,fasta_file)
        alignment= AlignIO.read(well_path, "fasta")
        summary_align = AlignInfo.SummaryInfo(alignment)
        consensus = (summary_align.gap_consensus (threshold=threshold, ambiguous=ambiguous_base, require_multiple=False)).replace('-','').upper()
        
        raw_consensus_seqs.append(consensus)
        raw_consensus_lengths.append(len(consensus))
        raw_amb_char_counts.append(consensus.count(ambiguous_base))
        
        read_counts.append(get_read_count_from_fasta(well_path))
    
    raw_consensus_df['well'] = wells
    raw_consensus_df['raw_consensus'] = raw_consensus_seqs
    raw_consensus_df['raw_consensus_length'] = raw_consensus_lengths
    raw_consensus_df['raw_amb_char_count'] = raw_amb_char_counts
    raw_consensus_df['read_count'] = read_counts
    raw_consensus_df['plate'] = plate
    
    # sort by well
    raw_consensus_df = raw_consensus_df.sort_values(by=['well'])
    
    # save to csv file
    raw_consensus_df.to_csv(output_file,index=False)
    print(f"Raw consensus csv file created for {plate}")
    print(f"Output file: {output_file}")
    
    return None
        
        
   

############################## Trim  & Reconstruct ##############################


def pre_trimming_reconstruction_json_update(run_json_file_path:str,trim_dictionary:dict,reconstruct_dictionary:dict):
    
    """
    Arguments:
    - run_json_file_path: path to the json file
    - trim_dictionary: dictionary with the trim sequences for each plate
    - reconstruct_dictionary: dictionary with the reconstruct sequences for each plate
    
    Actions:
    - Updates the json file with the trim and reconstruct sequences for each plate
    
    Returns:
    None
    """
    
    
    #load json file
    run_dictionary = load_json_file(run_json_file_path)
    
    for item in run_dictionary.items():
        
        plate = item[0]
        if plate == 'run_info':
            pass
        else:
            
            if plate in trim_dictionary.keys():
                if 'five_prime' in trim_dictionary[plate].keys():
                    run_dictionary[plate]['five_prime_trim_seq'] = trim_dictionary[plate]['five_prime']
                    
                if 'three_prime' in trim_dictionary[plate].keys():
                    run_dictionary[plate]['three_prime_trim_seq'] = trim_dictionary[plate]['three_prime']
            
            if plate in reconstruct_dictionary.keys():
                if 'five_prime' in reconstruct_dictionary[plate].keys():
                    run_dictionary[plate]['five_prime_reconstruct_seq'] = reconstruct_dictionary[plate]['five_prime']
                    
                if 'three_prime' in reconstruct_dictionary[plate].keys():
                    run_dictionary[plate]['three_prime_reconstruct_seq'] = reconstruct_dictionary[plate]['three_prime']
     
    
     
    
    #save to json file
    write_json_file(run_dictionary,run_json_file_path)
            
    return None


def trim_and_reconstruct_consesnsus(run_json_file_path,consensus_folder_path:str):
    
    """
    Arguments:
    - run_json_file_path: path to the json file
    - consensus_folder_path: path to the consensus output folder
    
    Actions:
    - Trims and reconstructs the consensus for each plate in the run
    - Writes the trimmed reconstructed consensus csv file for each plate
    - Writes the run level trimmed reconstructed consensus csv file
    - Updates the json file with the trimmed reconstructed consensus csv file for each plate
    - Updates the json file with the run level trimmed reconstructed consensus csv file
    
    Returns:
    None
    """
    
    #load json file
    run_dictionary = load_json_file(run_json_file_path)
    
    trimmed_reconstructed_folder_path = os.path.join(consensus_folder_path,'trimmed_reconstructed_consensus')
    create_folder(trimmed_reconstructed_folder_path)
    
    ambiguous_base = run_dictionary["run_info"]["consensus_ambiguous_base"]
    run_name = run_dictionary["run_info"]["run_name"]
    output_directory = run_dictionary["run_info"]["output_directory"]
    
    
    for item in run_dictionary.items():
        
        plate = item[0]
        
        if plate == 'run_info':
            pass
        else:
            # load raw consensus csv file
            raw_consensus_csv_file_path = run_dictionary[plate]['plate_raw_consensus_csv_file_path']
            trimmed_reconstructed_file_path = os.path.join(trimmed_reconstructed_folder_path,plate+'_trimmed_reconstructed_consensus.csv')
            # trim and reconstruct
            five_prime_trim_seq = run_dictionary[plate]['five_prime_trim_seq']
            three_prime_trim_seq = run_dictionary[plate]['three_prime_trim_seq']
            five_prime_reconstruct_seq = run_dictionary[plate]['five_prime_reconstruct_seq']
            three_prime_reconstruct_seq = run_dictionary[plate]['three_prime_reconstruct_seq']
            
            
            # run trimming and reconstruction
            plate_trim_and_reconstruct(raw_consensus_csv_file_path, trimmed_reconstructed_file_path, five_prime_trim_seq, three_prime_trim_seq, five_prime_reconstruct_seq, three_prime_reconstruct_seq,ambiguous_base)
            
            #update dictionary
            run_dictionary[plate]['plate_trimmed_reconstructed_consensus_csv_file_path'] = trimmed_reconstructed_file_path
            
    
    # save to json file
    write_json_file(run_dictionary,run_json_file_path)
    
    # create run level raw consensus csv file
    processed_run_df = pd.DataFrame
    
    # loop over plates to concat the trimmed reconstructed consensus csv files
    for item in run_dictionary.items():
        plate = item[0]
        if plate == 'run_info':
            pass
        else:
            plate_df = pd.read_csv(run_dictionary[plate]['plate_trimmed_reconstructed_consensus_csv_file_path'])
            if processed_run_df.empty:
                processed_run_df = plate_df
            else:
                processed_run_df = pd.concat([processed_run_df,plate_df])
    
    # save run level trimmed reconstructed consensus csv file
    run_consensus_csv_file_path = os.path.join(output_directory,run_name,run_name+'_trimmed_reconstructed_consesus.csv')
    processed_run_df.to_csv(run_consensus_csv_file_path,index=False)
    run_dictionary["run_info"]["run_trimmed_reconstructed_consensus_csv_file_path"] = run_consensus_csv_file_path
    
    print(f"Trimming and reconstruction completed for consensus sequences of {run_name}")
    print(f" Aggregate csv file: {run_consensus_csv_file_path}")
    
    # update json file
    write_json_file(run_dictionary,run_json_file_path) 
            
            
            
    return None



def plate_trim_and_reconstruct(raw_consensus_csv_file_path:str, trimmed_reconstructed_file_path:str, five_prime_trim_seq, three_prime_trim_seq, five_prime_reconstruct_seq, three_prime_reconstruct_seq,ambiguous_base):
    
    """
    Arguments:
    - raw_consensus_csv_file_path: path to the raw consensus csv file
    - trimmed_reconstructed_file_path: path to the trimmed reconstructed consensus csv file
    - five_prime_trim_seq: five prime trim sequence
    - three_prime_trim_seq: three prime trim sequence
    - five_prime_reconstruct_seq: five prime reconstruct sequence
    - three_prime_reconstruct_seq: three prime reconstruct sequence
    - ambiguous_base: ambiguous base for consensus
    
    Actions:
    - Trims and reconstructs the consensus for each plate in the run
    - Writes the trimmed reconstructed consensus csv file for each plate
    
    Returns:
    None
    """
    
    
    # load raw consensus csv file
    raw_consensus_df = pd.read_csv(raw_consensus_csv_file_path)
    
    processed_df = raw_consensus_df.copy()
    
    # add columns 'trimmed_seq'  'reconstructed_seq'  'post_reconstruction_length' 'post_reconstruction_amb_char_count' to the processed_df
    processed_df['trimmed_consensus'] = ''
    processed_df['reconstructed_consensus'] = ''
    processed_df['post_reconstruction_length'] = ''
    processed_df['post_reconstruction_amb_char_count'] = ''
    
    
    for index, row in processed_df.iterrows():
        
        consensus = row['raw_consensus']
        
        # trim the forward read
        start = 0
        end = len(consensus)
        # trim the 5' end
        if five_prime_trim_seq:   
            alignment_5_prime = edlib.align(five_prime_trim_seq, consensus, mode="HW", task="locations")
            start = alignment_5_prime['locations'][0][1] + 1
            
            
        # trim the 3' end
        if three_prime_trim_seq:
            alignment_3_prime = edlib.align(three_prime_trim_seq, consensus, mode="HW", task="locations")
            end = alignment_3_prime['locations'][0][0]
            
        trimmed_seq = consensus[start:end]
        processed_df.at[index,'trimmed_consensus'] = trimmed_seq
      
        # reconstruct the forward read
        reconstructed_seq = processed_df.at[index,'trimmed_consensus']
        if five_prime_reconstruct_seq:
            reconstructed_seq = five_prime_reconstruct_seq + reconstructed_seq
        
        if three_prime_reconstruct_seq:
            reconstructed_seq = reconstructed_seq + three_prime_reconstruct_seq
        
        processed_df.at[index,'reconstructed_consensus'] = reconstructed_seq
        
        processed_df.at[index,'post_reconstruction_length'] = len(reconstructed_seq)
        processed_df.at[index,'post_reconstruction_amb_char_count'] = reconstructed_seq.count(ambiguous_base)
        
        # sort by well
        processed_df = processed_df.sort_values(by=['well'])
        
        # save to csv file
        processed_df.to_csv(trimmed_reconstructed_file_path,index=False)
        
    print (f"Trimming and reconstruction of consensus sequences completed for {trimmed_reconstructed_file_path}")
    print (f"Output csv file: {trimmed_reconstructed_file_path}")
        
    return None
    
