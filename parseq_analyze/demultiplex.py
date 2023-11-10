from .utils import load_json_file, write_json_file, create_folder, globalize
import pandas as pd
import os
from multiprocessing import Pool
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import edlib
import subprocess

def check_freebarcodes_installed():
    """
    Checks if freebarcodes is installed
    """
    try:
        subprocess.run(['freebarcodes', '--version'], check=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
    except subprocess.CalledProcessError:
        print("Freebarcodes is not installed or not accessible from the system's PATH.. Please install freebarcodes v2.1 from https://github.com/hawkjo/freebarcodes")
    
    return None
    

def pre_demultiplex_json_update(run_json_file_path:str, barcodes_csv_path:str, barcode_search_error:int, fwd_read_constant_sequence_dict:dict):
    
    # load json file
    run_dictionary = load_json_file(run_json_file_path)
        
    for item in run_dictionary.items():
        
        plate = item[0]
        if plate == 'run_info':
            pass
        else:
            run_dictionary[plate]["cnst_fwd_read"] = fwd_read_constant_sequence_dict[plate]
            
    run_dictionary["run_info"]["barcodes_csv"] = barcodes_csv_path
    run_dictionary["run_info"]["freebarcodes_search_error"] = barcode_search_error
    
    #save to json file
    write_json_file(run_dictionary,run_json_file_path)
    
    return None


      
    
def create_freebarcodes_file(run_json_file_path:str,freebarcodes_output_path:str):
    
    """
    Arguments:
    - run_json_file_path: path to the json file
    
    Actions:
    - Creates freebarcodes file from the barcodes csv file
    - Writes the freebarcodes file to the output directory
    
    Returns:
    None
    """
    
    # load json
    run_dictionary = load_json_file(run_json_file_path)
    
    # get barcodes csv path
    barcodes_csv_path = run_dictionary["run_info"]["barcodes_csv"]
    freebarcodes_error = run_dictionary["run_info"]["freebarcodes_search_error"]

    # import barcodes csv file as df
    df = pd.read_csv(barcodes_csv_path,header=0)
    df = df[['well', 'barcode']]
    
    # check if all barcodes are the same length
    barcode_len = len(df['barcode'][0])
    for index, row in df.iterrows():
        if len(row['barcode']) != barcode_len:
            raise ValueError(f"Barcode at index {index} is not the same length as the first barcode. All barcodes have to be same length")
    
    # create freebarcodes file path and add it to json file
    output_file=os.path.join(freebarcodes_output_path,f"freebarcodes{barcode_len}-{freebarcodes_error}.txt")
    run_dictionary["run_info"]["freebarcodes_file"] = output_file
    write_json_file(run_dictionary,run_json_file_path)
    
    # write freebarcodes text file 
    with open(output_file, 'w') as f:
        for index, row in df.iterrows():
            f.write(f"{row['barcode']}\n")
    
    print(f"Freebarcodes file created at {output_file}")
    return output_file



def freebarcodes_decode(run_json_file_path, decode_output_path:str):
    
    """
    Arguments:
    - run_json_file_path: path to the json file
    - decode_output_path: path to the output folder
    
    Actions:
    - Runs freebarcodes decode command for each plate in the run
    - Writes the freebarcodes folder output path to the json file
    
    Returns:
    None
    """
    
    #load json file
    run_dictionary = load_json_file(run_json_file_path)
    free_barcodes_file = run_dictionary["run_info"]["freebarcodes_file"]
    
    #loop over plates
    for item in run_dictionary.items():
        
        plate = item[0]
        if plate == 'run_info': #skip run_info
            pass
        else:
            
            input_fastq_path = os.path.join(run_dictionary[plate]['fastp_folder_output_path'],plate+'.fastq')
            run_dictionary[plate]["freebarcodes_decoded_file_path"]=os.path.join(decode_output_path,plate+"_decoded.txt")

            freebarcodes_decode_cmd(input_fastq_path, free_barcodes_file, decode_output_path)
            print(f"Freebarcodes decode completed for file {plate}")
            print(f"Output folder: {decode_output_path}")
    
    # write to json file
    write_json_file(run_dictionary,run_json_file_path)
        
    return None
             

def freebarcodes_decode_cmd(fastq_path:str, free_barcodes_file:str, decode_output_path:str):
    
    """
    Arguments:
    - fastq_path: path to the fastq file
    - free_barcodes_file: path to the freebarcodes file
    - decode_output_path: path to the output folder
    
    Actions:
    - Runs freebarcodes decode command
    
    Returns:
    None
    """
    # check if input file exists
    if not os.path.exists(fastq_path):
        raise ValueError(f"File {fastq_path} does not exist")
    
    # check if output folder exists
    if not os.path.exists(decode_output_path):
        raise ValueError(f"Folder {decode_output_path} does not exist")
    
    cmd = f"freebarcodes decode {free_barcodes_file} {fastq_path} --output-dir={decode_output_path}".split()
    subprocess.run(cmd)
    
    return None


def separate_freebarcodes_decoded_file_into_wells(run_json_file_path:str,demultiplexed_folder_path:str):
    
    """
    Arguments:
    - run_json_file_path: path to the json file
    - demultiplexed_folder_path: path to the demultiplexed folder
    
    Actions:
    - Separates freebarcodes_decoded_file_into_wells for each plate in the run
    - Outputs separated well files to the demultiplexed folder for each of the plates
    - Outputs a csv file with the number of reads per well for each plate
    - Outputs a csv file with the demultiplexed fwd reads for each plate
    - Writes the demultiplexed folder output path to the json file
    - Runs multiprocessing on the plate level

    Returns:
    None
    """
    
    
    @globalize
    def pooled_process(plate, freebarcodes_decoded_path, cnst_fwd_read, barcodes_csv_path,demultiplexed_output_path,demultiplexed_folder_path):
   
        # import decoded file into df
        df = pd.read_csv(freebarcodes_decoded_path, sep='\t', header=None)
        df.columns = ['sequence_name','barcode','sequence']

        # write demultiplexed files
        print(f"Separating freebarcodes decoded file into wells for {plate}")
        create_folder(demultiplexed_output_path)
       
        df, reads_df= write_demultiplexed_files(df, barcodes_csv_path,cnst_fwd_read,demultiplexed_output_path)

        # save df and reads_df to plate folder as csv
        df.to_csv(os.path.join(demultiplexed_folder_path,f"{plate}_demultiplexed.csv"),index=False)
        reads_df.to_csv(os.path.join(demultiplexed_folder_path,f"{plate}_well_number_of_reads.csv"),index=False)        


        print(f"Freebarcodes decoded file separated into wells for {plate}")
        print(f"Output fasta folder: {demultiplexed_output_path}")
        print(f" Output csv files: {demultiplexed_folder_path}")
        
        return None
    
    
    #load json file
    run_dictionary = load_json_file(run_json_file_path)  
    barcodes_csv_path = run_dictionary["run_info"]["barcodes_csv"]
    multiprocessing_cores = run_dictionary["run_info"]["multiprocessing_cores"]
    
    plates = []
    freebarcodes_decoded_paths = []
    cnst_fwd_reads = []
    barcodes_csv_paths = []
    demultiplexed_output_paths = []
    demultiplexed_folder_paths = []
    
    #loop over plates to get and set variables from the json file
    for item in run_dictionary.items():
        
        
        plate = item[0]
        if plate == 'run_info':
            pass
        else:
            plates.append(plate)
            freebarcodes_decoded_paths.append(run_dictionary[plate]["freebarcodes_decoded_file_path"])
            cnst_fwd_reads.append(run_dictionary[plate]["cnst_fwd_read"])
            barcodes_csv_paths.append(barcodes_csv_path)
            
            plate_demultiplexed_output_path = os.path.join(demultiplexed_folder_path,plate)
            run_dictionary[plate]["demultiplexed_folder_output_path"] = plate_demultiplexed_output_path
            demultiplexed_output_paths.append(plate_demultiplexed_output_path)
            demultiplexed_folder_paths.append(demultiplexed_folder_path)
    
    # write to json file (demultiplexed output paths)
    write_json_file(run_dictionary,run_json_file_path)
    
    # zip input to pooled_process for multiprocessing       
    zipped_input = zip(plates, freebarcodes_decoded_paths, cnst_fwd_reads, barcodes_csv_paths,demultiplexed_output_paths,demultiplexed_folder_paths)
    
    # run multiprocessing
    with Pool(multiprocessing_cores) as pool:
        pool.starmap(pooled_process, zipped_input)
        
        
    return None
    



def write_demultiplexed_files(df:pd.DataFrame, barcodes_csv_path:str,cnst_fwd_read:str,output_folder:str):
    
    """
    Arguments:
    - df: dataframe with the freebarcodes decoded file
    - barcodes_csv_path: path to the barcodes csv file
    - cnst_fwd_read: constant forward read
    - output_folder: path to the plate output folder for writing demultiplexed well files
    
    Actions:
    - Assigns the well to each sequence
    - reverses the sequence if it is reverse complement
    - Writes the demultiplexed well files to the output folder
    
    Returns:
    - df: dataframe with the well demultiplexed fwd reads for each plate
    - reads_df: dataframe with the number of reads per well for each plate
    """
    
    
        
    # create 384 plates
    wells = [row+column for row in ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P'] for column in ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24']]
    
    #import barcodes csv file
    barcodes_df = pd.read_csv(barcodes_csv_path,header=0)
    barcodes_df = barcodes_df[['well', 'barcode']]
    
    # add well,fwd/rev, and fwd_read to column to df
    df['well'] = ''
    df['fwd/rev'] = ''
    df['fwd_read'] = ''
    
    # match the barcode to the well
    
    for index, row in df.iterrows():
        barcode = row['barcode']
        well = barcodes_df[barcodes_df['barcode']==barcode]['well'].values[0]
        df.at[index,'well'] = well
        
        read=row['sequence']
        read_rc = Seq(read).reverse_complement()
        
        if edlib.align(cnst_fwd_read, read, mode="HW", task="locations")['editDistance'] < edlib.align(cnst_fwd_read, read_rc,mode="HW", task="locations")['editDistance']:
            df.at[index,'fwd/rev'] = 'fwd'
            df.at[index,'fwd_read'] = read
        else:
            df.at[index,'fwd/rev'] = 'rev'
            df.at[index,'fwd_read'] = read_rc
    
    #write to fasta file

    reads_df =pd.DataFrame(columns=['well','number_of_reads'])# dataframe for the number of reads per well
    for well in wells:
        well_df = df[df['well']==well]
        saved_records = []
        for index, row in well_df.iterrows():
            saved_records.append(SeqRecord(Seq(row['fwd_read']), id=row['sequence_name'], description=""))
        num_reads = len(saved_records)
        reads_df = pd.concat([reads_df,pd.DataFrame([[well,num_reads]],columns=['well','number_of_reads'])])
        if num_reads > 0:
            SeqIO.write(saved_records, os.path.join(output_folder,f"{well}.fasta"), "fasta")
    
    # sort df by well
    df = df.sort_values(by=['well'])

    return df, reads_df