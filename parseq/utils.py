import json
import os
import subprocess
import sys
import uuid
from Bio import SeqIO
import numpy as np


def create_folder(folder_path:str):
    """Creates Folder if it does not exist"""
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    return None


def load_json_file(json_path:str):
    """Loads json file into dictionary and returns dictionary"""
    with open(json_path, 'r') as fp:
        run_dictionary = json.load(fp)
    return run_dictionary


def write_json_file(run_dictionary:dict,json_path:str):
    """Writes dictionary to json file"""
    with open(json_path, 'w') as fp:
        json.dump(run_dictionary, fp, indent=4)
    return None


def copy_file_cmd(source_file_path:str, destination_folder:str, destination_file_name:str):
    """
    Copies file from source to destination
    Allows user to change the name of the file in the destination folder
    """
    destination_file_path = os.path.join(destination_folder,destination_file_name)
    copy_cmd=(f"cp {source_file_path} {destination_file_path}")
    subprocess.run(copy_cmd,shell=True)
    print(f"Copied {source_file_path} to {destination_file_path}")
    return None


def unzip_file_cmd(source_file_path:str, destination_folder:str, destination_file_name:str):
    """
    Unzips file from source to destination
    Allows user to change the name of the file in the destination folder
    """
    destination_file_path = os.path.join(destination_folder,destination_file_name)
    unzip_cmd=(f"gunzip {source_file_path} -c > {destination_file_path}")
    subprocess.run(unzip_cmd,shell=True)
    print(f"Unzipped {source_file_path} to {destination_file_path}")
    return None


def globalize(func):
    """
    Decorator that makes a function global
    Useful for multiprocessing nested functions
    """
    def result(*args, **kwargs):
        return func(*args, **kwargs)
    result.__name__ = result.__qualname__ = uuid.uuid4().hex
    setattr(sys.modules[result.__module__], result.__name__, result)
    
    return result


def fastq_file_read_stats(input_file_abs_path:str):
    
    """
    Arguments:
    - input_file_abs_path: absolute path to the fastq file
    
    Actions:
    - Calculates the length of each sequence in the fastq file
    - Calculates the mean, median, min, max, and std of the lengths
    
    Returns:
    - lengths: list of lengths of sequences
    - stats_dict: dictionary with stats of the plate
    """
    #get the legths of sequences, along with some statistics
    lengths = [len(rec) for rec in SeqIO.parse(input_file_abs_path, "fastq")]
    num_sequences = len(lengths)
    mean_length = np.mean(lengths)
    median_length = np.median(lengths)
    min_length = min(lengths)
    max_length = max(lengths)
    std_length = np.std(lengths)
    stats_dict = {
        
        "num_sequences":num_sequences,
        "mean_length":mean_length, 
        "median_length":median_length, 
        "min_length":min_length, 
        "max_length":max_length, 
        "length_std":std_length
    }
    
    return lengths,stats_dict



def run_fastp_cmd (plate:str,fastp_source_file:str, fastp_output_path:str, length_filtering:tuple, NGS_read_quality_filtering:int):
    
    """
    Arguments:
    - fastp_source_file: path to the fastq file
    - fastp_output_path: path to the fastp output folder
    - length_filtering: tuple with min and max length
    - NGS_read_quality_filtering: NGS minimum phred quality score
    
    Actions:
    - Runs fastp command
    
    Returns:
    None
    """
    # check if input file exists
    if not os.path.exists(fastp_source_file):
        raise ValueError(f"File {fastp_source_file} does not exist")
    
    # check if output folder exists
    if not os.path.exists(fastp_output_path):
        raise ValueError(f"Folder {fastp_output_path} does not exist")
    
    plate_output_path = os.path.join(fastp_output_path,plate)
        
    cmd=(f"fastp -A -i {fastp_source_file} --length_required {int(length_filtering[0])} --length_limit {int(length_filtering[1])} -o {plate_output_path}.fastq  -j {plate_output_path}.json -h {plate_output_path}.html -q {NGS_read_quality_filtering}").split()
    print(cmd)
    subprocess.run(cmd)
    
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


def get_read_count_from_fasta(fasta_path:str):
    
    """
    gets the number of reads in a fasta file
    """ 
    with open(fasta_path, "r") as handle:
        read_count = handle.read().count(">")    
        
    return read_count
