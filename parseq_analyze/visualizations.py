# Description: This file contains functions to create visualizations for raw consensus and trimmed reconstructed consensus data.

from .utils import load_json_file, write_json_file, create_folder
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings


def fastq_sequence_length_histogram(plate,lengths, stats_dict, output_directory_abs_path:str,show_plot:bool=False): 
    
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
    
    plot_path = os.path.join(output_directory_abs_path,f"{plate}_sequence_length_histogram.png")
    plt.savefig(plot_path)
    
    if show_plot:
        plt.show()
    print(f"Sequence length histogram for {plate} saved to {output_directory_abs_path}")
    
    plt.close()
    return plot_path


def create_visualizations_raw_consensus(json_file_path:str, visualizations_directory:str):
    """
    Arguments:
    - json_file_path: path to json file
    - visualizations_directory: path to visualizations directory
    
    Actions:
    - Creates visualizations for raw consensus data
    - Run level visualizations:
        - well average read count heatmap
        - well average raw consensus seq length heatmap
        - well average raw consensus ambiguous count heatmap
        - plate average read count bar chart
        - plate average raw consensus seq length bar chart
        - plate average raw consensus ambiguous count bar chart
        
    - Plate level visualizations:
        - well read count heatmap
        - well read count histogram
        - well raw consensus seq length heatmap
        - well raw consensus seq length histogram
        - well raw consensus ambiguous count heatmap
        - well raw consensus ambiguous count histogram
        
    
    Returns:
    None
    """
    
    # load json file
    run_dictionary = load_json_file(json_file_path)
    run_dictionary["run_info"]["visualization_folder_output_path"] = visualizations_directory
    run_name = run_dictionary["run_info"]["run_name"]
    raw_consensus_csv_path = run_dictionary["run_info"]["run_raw_consensus_csv_file_path"]
    write_json_file(run_dictionary, json_file_path)
    
    # check if csv file exists:
    if not os.path.isfile(raw_consensus_csv_path):
        raise ValueError(f"File {raw_consensus_csv_path} does not exist. Please run alignment and consensus first.")
    
    # load csv file
    run_raw_consensus_df = pd.read_csv(raw_consensus_csv_path)

    
    # create visualizations directory
    raw_consensus_visualizations_directory = os.path.join(visualizations_directory,"raw_consensus_visualizations")
    create_folder(raw_consensus_visualizations_directory)
    
    ##### create run level visualizations #####
    
    # get well level averages for run
    run_df = run_raw_consensus_df.groupby(['well']).mean(numeric_only=True)
    run_df.reset_index(inplace=True)
    
    # well level average counts
    plotdf = run_df[["well","read_count"]]
    # set column name to value
    plotdf = plotdf.rename(columns={"read_count":"value"})
    # create average read count heatmap
    name = "well_average_read_count"
    create_384_plate_heatmap(plotdf, raw_consensus_visualizations_directory, run_name , name)
    
    # well level average length
    plotdf = run_df[["well","raw_consensus_length"]]
    # set column name to value
    plotdf = plotdf.rename(columns={"raw_consensus_length":"value"})
    # create average length heatmap
    name = "well_average_raw_consensus_seq_length"
    create_384_plate_heatmap(plotdf, raw_consensus_visualizations_directory, run_name , name)
    
    # well level average ambiguous count
    plotdf = run_df[["well","raw_amb_char_count"]]
    # set column name to value
    plotdf = plotdf.rename(columns={"raw_amb_char_count":"value"})
    # create average ambiguous count heatmap
    name = "well_average_raw_consensus_amb_count"
    create_384_plate_heatmap(plotdf, raw_consensus_visualizations_directory, run_name , name)
    
    
    
    # get plate level averages for run
    run_df = run_raw_consensus_df.groupby(['plate']).mean(numeric_only=True)
    run_df.reset_index(inplace=True)
    
    
    # plate level average read count
    plotdf = run_df[["plate","read_count"]]
    # set column name to value
    plotdf = plotdf.rename(columns={"read_count":"value"})
    # create plate average read count bar chart
    name = "plate_average_read_count"
    create_bar_chart(plotdf, raw_consensus_visualizations_directory, run_name , name)
    
    # plate level average seq length
    plotdf = run_df[["plate","raw_consensus_length"]]
    # set column name to value
    plotdf = plotdf.rename(columns={"raw_consensus_length":"value"})
    # create plate average seq length bar chart
    name = "plate_average_raw_consensus_seq_length"
    create_bar_chart(plotdf, raw_consensus_visualizations_directory, run_name , name)
    
    # plate level average ambiguous count
    plotdf = run_df[["plate","raw_amb_char_count"]]
    # set column name to value
    plotdf = plotdf.rename(columns={"raw_amb_char_count":"value"})
    # create plate average ambiguous count bar chart
    name = "plate_average_raw_consensus_ambiguous_count"
    create_bar_chart(plotdf, raw_consensus_visualizations_directory, run_name , name)
    
    print("Run level visualizations for raw consensus sequences created.")
    print("Output directory: ", raw_consensus_visualizations_directory)
    
    
    ##### create plate level visualizations #####    
    
    for item in run_dictionary.items():
        
        plate = item[0]
        if plate == "run_info": # skip run_info
            pass
        else:
            plate_df = run_raw_consensus_df[run_raw_consensus_df["plate"] == plate]
            
            # create visualizations directory for plate
            plate_visualizations_directory = os.path.join(raw_consensus_visualizations_directory,plate)
            create_folder(plate_visualizations_directory)
            
            #read_count heatmap and histogram
            # take the columns "well" and "read_count"
            plotdf = plate_df[["well","read_count"]]
            # set column name to value
            plotdf = plotdf.rename(columns={"read_count":"value"})
            # create ambiguous count heatmap
            name = "read_count"
            create_384_plate_heatmap(plotdf, plate_visualizations_directory, plate, name)
            create_histogram(plotdf, plate_visualizations_directory, plate, name)
            
            
            
            #raw_consensus_length heatmap
            # take the columns "well" and "raw_amb_char_count"
            plotdf = plate_df[["well","raw_consensus_length"]]
            # set column name to value
            plotdf = plotdf.rename(columns={"raw_consensus_length":"value"})
            # create ambiguous count heatmap
            name = "raw_consensus_length"
            create_384_plate_heatmap(plotdf, plate_visualizations_directory, plate, name)
            create_histogram(plotdf, plate_visualizations_directory, plate, name)
            
            #raw_amb_char_count heatmap
            # take the columns "well" and "raw_amb_char_count"
            plotdf = plate_df[["well","raw_amb_char_count"]]
            # set column name to value
            plotdf = plotdf.rename(columns={"raw_amb_char_count":"value"})
            # create ambiguous count heatmap
            name = "raw_consensus_amb_count"
            create_384_plate_heatmap(plotdf, plate_visualizations_directory, plate, name)
            create_histogram(plotdf, plate_visualizations_directory, plate, name)
    
    print("Plate level visualizations for raw consensus sequences created.")
    print("Output directory: ", raw_consensus_visualizations_directory)       
    return None
            
            


def create_visualizations_trimmed_reconstructed_consensus(json_file_path:str, visualizations_directory:str):
    
    """
    Arguments:
    - json_file_path: path to json file
    - visualizations_directory: path to visualizations directory
    
    Actions:
    - Creates visualizations for trimmed reconstructed consensus data
    - Run level visualizations:
        - well average read count heatmap
        - well average trimmed reconstructed consensus seq length heatmap
        - well average trimmed reconstructed consensus ambiguous count heatmap
        - plate average read count bar chart
        - plate average trimmed reconstructed consensus seq length bar chart
        - plate average trimmed reconstructed consensus ambiguous count bar chart
        
    - Plate level visualizations:
        - well read count heatmap
        - well read count histogram
        - well trimmed reconstructed consensus seq length heatmap
        - well trimmed reconstructed consensus seq length histogram
        - well trimmed reconstructed consensus ambiguous count heatmap
        - well trimmed reconstructed consensus ambiguous count histogram
        
    
    Returns:
    None
    """
    
    
    
    # load json file
    run_dictionary = load_json_file(json_file_path)
    run_dictionary["run_info"]["visualization_folder_output_path"] = visualizations_directory
    run_name = run_dictionary["run_info"]["run_name"]
    recon_consensus_csv_path = run_dictionary["run_info"]["run_trimmed_reconstructed_consensus_csv_file_path"]
    write_json_file(run_dictionary, json_file_path)
    
    
    # check if csv file exists:
    if not os.path.isfile(recon_consensus_csv_path):
        raise ValueError(f"File {recon_consensus_csv_path} does not exist. Please run trimming and reconstruction first.")
    
    # load csv file
    run_recon_consensus_df = pd.read_csv(recon_consensus_csv_path)
    
    # create visualizations directory
    recon_consensus_visualizations_directory = os.path.join(visualizations_directory,"trimmed_recon_consensus_visualizations")
    create_folder(recon_consensus_visualizations_directory)
    
    
    ##### create run level visualizations #####   
    
    # get well level averages for run
    run_df = run_recon_consensus_df.groupby(['well']).mean(numeric_only=True)
    run_df.reset_index(inplace=True)
    
    # well level average counts
    plotdf = run_df[["well","read_count"]]
    # set column name to value
    plotdf = plotdf.rename(columns={"read_count":"value"})
    # create average read count heatmap
    name = "well_average_read_count"
    create_384_plate_heatmap(plotdf, recon_consensus_visualizations_directory, run_name , name)
    
    # well level average length
    plotdf = run_df[["well","post_reconstruction_length"]]
    # set column name to value
    plotdf = plotdf.rename(columns={"post_reconstruction_length":"value"})
    # create average length heatmap
    name = "well_average_post_reconstruction_length"
    create_384_plate_heatmap(plotdf, recon_consensus_visualizations_directory, run_name , name)
    
    # well level average ambiguous count
    plotdf = run_df[["well","post_reconstruction_amb_char_count"]]
    # set column name to value
    plotdf = plotdf.rename(columns={"post_reconstruction_amb_char_count":"value"})
    # create average ambiguous count heatmap
    name = "well_average_post_reconstruction_consensus_amb_count"
    create_384_plate_heatmap(plotdf, recon_consensus_visualizations_directory, run_name , name)
    
    
    # get plate level averages for run
    run_df = run_recon_consensus_df.groupby(['plate']).mean(numeric_only=True)
    run_df.reset_index(inplace=True)
    
    
    # plate level average read count
    plotdf = run_df[["plate","read_count"]]
    # set column name to value
    plotdf = plotdf.rename(columns={"read_count":"value"})
    # create plate average read count bar chart
    name = "plate_average_read_count"
    create_bar_chart(plotdf, recon_consensus_visualizations_directory, run_name , name)
    
    # plate level average seq length
    plotdf = run_df[["plate","post_reconstruction_length"]]
    # set column name to value
    plotdf = plotdf.rename(columns={"post_reconstruction_length":"value"})
    # create plate average seq length bar chart
    name = "plate_average_post_reconstruction_consensus_seq_length"
    create_bar_chart(plotdf, recon_consensus_visualizations_directory, run_name , name)
    
    # plate level average ambiguous count
    plotdf = run_df[["plate","post_reconstruction_amb_char_count"]]
    # set column name to value
    plotdf = plotdf.rename(columns={"post_reconstruction_amb_char_count":"value"})
    # create plate average ambiguous count bar chart
    name = "plate_average_post_reconstruction_consensus_amb_count"
    create_bar_chart(plotdf, recon_consensus_visualizations_directory, run_name , name)
    
    print("Run level visualizations for trimmed reconstructed consensus sequences created.")
    print("Output directory: ", recon_consensus_visualizations_directory)
    
    
    
    ##### create plate level visualizations #####   
    
    for item in run_dictionary.items():
        
        plate = item[0]
        if plate == "run_info":
            pass
        else:
            plate_df = run_recon_consensus_df[run_recon_consensus_df["plate"] == plate]
            
            # create visualizations directory for plate
            plate_visualizations_directory = os.path.join(recon_consensus_visualizations_directory,plate)
            create_folder(plate_visualizations_directory)
            
            #read_count heatmap and histogram
            # take the columns "well" and "read_count"
            plotdf = plate_df[["well","read_count"]]
            # set column name to value
            plotdf = plotdf.rename(columns={"read_count":"value"})
            # create ambiguous count heatmap
            name = "read_count"
            create_384_plate_heatmap(plotdf, plate_visualizations_directory, plate, name)
            create_histogram(plotdf, plate_visualizations_directory, plate, name)
            
            
            
            #raw_consensus_length heatmap
            # take the columns "well" and "raw_amb_char_count"
            plotdf = plate_df[["well","raw_consensus_length"]]
            # set column name to value
            plotdf = plotdf.rename(columns={"raw_consensus_length":"value"})
            # create ambiguous count heatmap
            name = "post_reconstruction_consensus_seq_length"
            create_384_plate_heatmap(plotdf, plate_visualizations_directory, plate, name)
            create_histogram(plotdf, plate_visualizations_directory, plate, name)
            
            #raw_amb_char_count heatmap
            # take the columns "well" and "raw_amb_char_count"
            plotdf = plate_df[["well","raw_amb_char_count"]]
            # set column name to value
            plotdf = plotdf.rename(columns={"raw_amb_char_count":"value"})
            # create ambiguous count heatmap
            name = "post_reconstruction_consensus_amb_count"
            create_384_plate_heatmap(plotdf, plate_visualizations_directory, plate, name)
            create_histogram(plotdf, plate_visualizations_directory, plate, name)
    
    print("Plate level visualizations for trimmed reconstructed consensus sequences created.")
    print("Output directory: ", recon_consensus_visualizations_directory)
           
    return None

          




def create_384_plate_heatmap (input_dataframe:pd.DataFrame,path_to_output_directory:str,entity:str, name:str):

    """
    Arguments:
    - input_dataframe: dataframe containing the data to be plotted (must contain columns 'well' and 'value')
    - path_to_output_directory: path to output directory
    - entity: name of the entity to be plotted
    - name: name of the plot
    
    Actions:
    - Creates a 384 plate heatmap of the input dataframe
    - Saves the heatmap in the output directory
    
    Returns:
    None
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
    
        df = input_dataframe.copy()
        # Sort the dataframe by 'well'
        df['row'] = df['well'].str[0]
        df['col'] = df['well'].str[1:].astype(int)
        df = df.sort_values(by=['row', 'col'])

        # Reshape the dataframe to wide format
        heatmap_data = df.pivot(columns='col', index='row', values='value') 



        # plot the heatmap
        plt.figure(figsize=(20, 10))
        cmap = sns.diverging_palette(220, 20, as_cmap=True)

        vmin = df['value'].min()
        vmax = df['value'].max()
        sns.heatmap(heatmap_data, annot=True, fmt="3.0f",cmap=cmap, vmin=vmin, vmax=vmax)
        plt.title(f"{entity}_{name}_384-well_heatmap") 
        plt.savefig(os.path.join(path_to_output_directory, f"{entity}_{name}_384_well_heatmap.png"))
        plt.close()
    
    return None



def create_histogram(input_dataframe:pd.DataFrame, path_to_output_directory:str, entity:str, name:str):
    
    """
    Arguments:
    - input_dataframe: dataframe containing the data to be plotted (must contain column  'value')
    - path_to_output_directory: path to output directory
    - entity: name of the entity to be plotted
    - name: name of the plot
    
    Actions:
    - Creates a histogram of the input dataframe
    - Saves the histogram in the output directory
    
    Returns:
    None
    """
   
   
    df=input_dataframe.copy()
   
    sns.histplot(df['value'], kde=True, bins=30)
    plt.title(f"{entity}_{name}_histogram")
    plt.xlabel("Value")
    plt.ylabel("Count")
    plt.savefig(os.path.join(path_to_output_directory, f"{entity}_{name}_histogram.png"))
    plt.close()
    return None


def create_bar_chart(input_dataframe, path_to_output_directory:str, entity:str, name:str):
    
    """
    Arguments:
    - input_dataframe: dataframe containing the data to be plotted (must contain columns 'plate' and 'value')
    - path_to_output_directory: path to output directory
    - entity: name of the entity to be plotted
    - name: name of the plot
    
    Actions:
    - Creates a bar chart of the input dataframe
    - Saves the bar chart in the output directory
    
    Returns:
    None
    """
   
   
    df=input_dataframe.copy()

    df.plot.bar(x='plate', y='value', rot=0, legend=False)
    plt.ylabel('Values')
    plt.title(f"{entity}_{name}_bar_chart")
    plt.xticks(rotation=45, fontsize=8)
    plt.savefig(os.path.join(path_to_output_directory, f"{entity}_{name}_bar_chart.png"))
    plt.close()
    return None
    