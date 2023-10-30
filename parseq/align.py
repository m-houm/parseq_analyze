from .utils import load_json_file, write_json_file, create_folder, globalize
import os
from Bio.Align.Applications import MuscleCommandline
from Bio.Align.Applications import MafftCommandline
from multiprocessing import Pool






def align_plates(run_json_file_path:str,alignment_output_path:str,alignment_algorithm:str='mafft'):
    
    """
    Arguments:
    - run_json_file_path: path to the json file
    - alignment_output_path: path to the alignment output folder
    - alignment_algorithm: alignment algorithm to use. Options are 'mafft' or 'muscle'
    
    Actions:
    - Runs alignment for each plate in the run
    - Writes the alignment folder output path to the json file for each plate
    
    Returns:
    None
    """
    
    
    # load json file
    run_dictionary = load_json_file(run_json_file_path)
    
    # check alignment algorithm is one of the two options
    assert alignment_algorithm in ['mafft','muscle'], "Alignment algorithm must be one of 'mafft' or 'muscle'"
    # update dictionary
    run_dictionary["run_info"]["alignment_algorithm"] = alignment_algorithm
    
    # get multiprocessing cores
    multiprocessing_cores = run_dictionary["run_info"]["multiprocessing_cores"]
    
    
    # loop over plates to set read and write alginment paths to json, and start alignment
    for item in run_dictionary.items():
        
        
        
        if item[0] == 'run_info':
            pass
        else:
            plate = item[0]
            input_folder = run_dictionary[plate]['demultiplexed_folder_output_path']
            output_folder = os.path.join(alignment_output_path,plate)
            # update dictionary
            run_dictionary[plate]['alignment_folder_output_path'] = output_folder
            create_folder(output_folder)
            
            if alignment_algorithm == 'mafft': # run mafft alignment
                print(f"Aligning {plate} with mafft")
                mafft_alignment_on_plate(input_folder, output_folder, multiprocessing_cores)
            
            elif alignment_algorithm == 'muscle': # run muscle alignment
                print(f"Aligning {plate} with muscle")
                muscle_alignment_on_plate(input_folder, output_folder, multiprocessing_cores)
                
            
            print(f"Alignment completed for file {plate}")
            print(f"Output folder: {output_folder}")
    
    # update json file        
    write_json_file(run_dictionary,run_json_file_path)
    
    return None


def muscle_alignment_on_plate(input_folder:str, output_folder:str, multiprocessing_cores:int):
    
    """
    Arguments:
    - input_folder: path to the input folder with the demultiplexed well files for each plate
    - output_folder: path to the output folder for the aligned well files for each plate
    - multiprocessing_cores: number of cores to use for multiprocessing
    
    Actions:
    - Runs muscle alignment for each well in the plate
    - Writes the aligned well files to the output folder
    - Runs multiprocessing on the well level (Multi-well alignment concurrently)
    
    Returns:
    None
    """
    
    
    @globalize # globalize the function to be able to use it in multiprocessing
    def muscle_wrapper (in_file:str,out_file:str):
        """muscle wrapper for multiprocessing"""
        print('Aligning: '+str(in_file))
        muscle_cline = MuscleCommandline( input=in_file, out=out_file, gapopen=0.5, gapextend=0.5)
        muscle_cline()
        return
    
    # get all fasta file names in the input folder 
    fasta_file_names = [file for file in os.listdir(input_folder) if file.endswith(".fasta")]
    fasta_file_names = sorted(fasta_file_names)
    input_fasta_file_paths = [os.path.join(input_folder,file) for file in fasta_file_names]
    output_fasta_file_paths = [os.path.join(output_folder,file) for file in fasta_file_names]
    
    # zip input to muscle_wrapper for multiprocessing 
    zipped_input = zip(input_fasta_file_paths, output_fasta_file_paths)
    
    # run multiprocessing
    with Pool(multiprocessing_cores) as pool:
        pool.starmap(muscle_wrapper, zipped_input)
        
    
    return None




def mafft_alignment_on_plate(input_folder:str, output_folder:str, multiprocessing_cores:int):
    
    """
    Arguments:
    - input_folder: path to the input folder with the demultiplexed well files for each plate
    - output_folder: path to the output folder for the aligned well files for each plate
    - multiprocessing_cores: number of cores to use for multiprocessing
    
    Actions:
    - Runs mafft alignment for each well in the plate
    - Writes the aligned well files to the output folder
    - Runs multiprocessing on the well level (Multi-well alignment concurrently)
    
    Returns:
    None
    """

    @globalize # globalize the function to be able to use it in multiprocessing
    def mafft_wrapper(in_file:str,out_file:str):
        """mafft wrapper for multiprocessing"""
        mafft_cline = MafftCommandline(input=in_file, auto=True)
        stdout, _ = mafft_cline()
        with open(out_file, "w") as handle:
            handle.write(stdout)
    
    # get all fasta file names in the input folder 
    fasta_file_names = [file for file in os.listdir(input_folder) if file.endswith(".fasta")]
    input_fasta_file_paths = [os.path.join(input_folder,file) for file in fasta_file_names]
    output_fasta_file_paths = [os.path.join(output_folder,file) for file in fasta_file_names]
    
    # zip input to muscle_wrapper for multiprocessing 
    zipped_input = zip(input_fasta_file_paths, output_fasta_file_paths)
    
    # run multiprocessing
    with Pool(multiprocessing_cores) as pool:
        pool.starmap(mafft_wrapper, zipped_input)
        
    
    return None