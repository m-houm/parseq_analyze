from utils import load_json_file, write_json_file, create_folder, run_fastp_cmd, fastq_file_read_stats
import os



def pre_fastp_update_json_file(json_file_path:str,fastp_output_path:str,length_filtering_dict:dict, ngs_read_quality_filtering_dict:dict):
    """
    Arguments:
    - json_file_path: path to the json file
    - fastp_output_path: path to the fastp output folder
    - length_filtering_dict: dictionary with plate name as key and length filtering as value
    - ngs_read_quality_filtering_dict: dictionary with plate name as key and NGS minimum phred quality score as value
    
    Actions:
    - Updates the json file with the fastp output path, length filtering, and NGS minimum phred quality score
    
    Returns:
    None
    """
    
    #load json file
    run_dictionary = load_json_file(json_file_path)
    
    for item in run_dictionary.items():

        plate = item[0]
        if plate == 'run_info': #skip run_info
            pass
        else:
           
            # add fastp output path
            run_dictionary[plate]['fastp_folder_output_path'] = os.path.join(fastp_output_path,plate)
            create_folder(run_dictionary[plate]['fastp_folder_output_path'])
            
            # check for read quality filtering
            if plate in ngs_read_quality_filtering_dict.keys():    
                #if yes, then use the defined values
                run_dictionary[plate]['NGS_minimum_quality_score'] = ngs_read_quality_filtering_dict[plate]['min_quality_score']
            else:
                pass
            
            #check if length filtering is defined for this plate
            if plate in length_filtering_dict.keys():    
                #if yes, then use the defined values
                min_length = length_filtering_dict[plate]['min']
                max_length = length_filtering_dict[plate]['max']
                run_dictionary[plate]['length_filtering'] = (min_length,max_length)  
            else:  
                #if not, then use mean length +/- 2*std as default
                min_length = run_dictionary[plate]['raw_file_stats']['mean_length'] - 2*run_dictionary[plate]['raw_file_stats']['length_std']
                max_length = run_dictionary[plate]['raw_file_stats']['mean_length'] + 2*run_dictionary[plate]['raw_file_stats']['length_std']
                run_dictionary[plate]['length_filtering'] = (min_length,max_length)
            
            
            
    #write to json file
    write_json_file(run_dictionary,json_file_path)
    
    return None
    


def run_fastp_on_plates(run_json_file_path:str):
    
    
    """ 
    Runs fastp for on raw sequencing files for each plate in the run.
    Takes all required parameters from the json file.
    Outputs: 
    - fastq file filtered for quality and length.
    - json file with the stats.
    - html file with the stats.
    
    """
    
    #load json file
    run_dictionary = load_json_file(run_json_file_path)
    
    
    for item in run_dictionary.items():
        
        plate = item[0]
        if plate == 'run_info':
            pass
        
        else:   
            fastp_source_file = run_dictionary[plate]['fastq_destination_file'] 
            fastp_output_path = run_dictionary[plate]['fastp_folder_output_path']
            #fastp_plate_output_path =os.path.join(fastp_output_path,plate)
            length_filtering = run_dictionary[plate]['length_filtering']
            NGS_read_quality_filtering = run_dictionary[plate]['NGS_minimum_quality_score']
            
            # run fastp command
            run_fastp_cmd (plate,fastp_source_file, fastp_output_path, length_filtering, NGS_read_quality_filtering)
            print(f"fastp completed for {plate}")
    
    return None


    

def post_fastp_stats_update_json_file(run_json_file_path:str):
    """
    Arguments:
    - json_file_path: path to the json file
    
    Actions:
    - Updates the json file with the post fastp stats of the plate fastq file
    
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
            input_file_abs_path = os.path.join(run_dictionary[plate]['fastp_folder_output_path'],plate+'.fastq')
            #get the legths of sequences, along with some statistics
            _, stats_dictionary = fastq_file_read_stats(input_file_abs_path)
            
            run_dictionary[plate]['post_fastp_stats']['num_sequences'] = stats_dictionary['num_sequences']
            run_dictionary[plate]['post_fastp_stats']['mean_length'] = stats_dictionary['mean_length']
            run_dictionary[plate]['post_fastp_stats']['median_length'] = stats_dictionary['median_length']
            run_dictionary[plate]['post_fastp_stats']['min_length'] = stats_dictionary['min_length']
            run_dictionary[plate]['post_fastp_stats']['max_length'] = stats_dictionary['max_length']
            run_dictionary[plate]['post_fastp_stats']['length_std'] = stats_dictionary['length_std']
            
    # write to json file
    write_json_file(run_dictionary,run_json_file_path)
        
    return None