# plot_alignment.py
# Ashley Mae Conard
# Last Mod. 6/23/2019
# Purpose: Plot alignment scores of input method

import matplotlib as mpl
mpl.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.lines import Line2D
import seaborn as sns
import glob, sys, os, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
sns.set_style('whitegrid')
sns.set_context("notebook", font_scale=2.5)

# Checking input arguments
if len(sys.argv) != 6 :
    print("Usage: python plot_alignment.py \n 1) METHOD_COMPARE (0-no, 1-yes) \n 2) /FULL/PATH/TO/ALIGNMENT/ (e.g. ../timeor/results/preprocess/ALIGNMENT_METHOD_FOLDER) \n 3) PAIRED_OR_NOT (0-single, 1-paired) \n 4) METHOD_NAME_1 (e.g. HISAT2) \n 5) METHOD_NAME_2 (e.g. Bowtie2 or NA)")
    sys.exit (1)    

# Arguments
METHOD_COMP = int(sys.argv[1]) # change to 1 if plotting 1 alignment tool vs. another (e.g. Bowtie2 vs. HISAT2)
INPUT_DIR = sys.argv[2]
PAIRED = int(sys.argv[3]) # change to 1 if paired-end sequencing
METHOD_NAME = str(sys.argv[4])
METHOD2_NAME = str(sys.argv[5])

def extract_alignment_val(align_file):
    # Openning alignment file and read all lines
    with open(align_file) as f:
        lines = f.readlines()

    # If paired-end, check different lines for conc and overall.
    if PAIRED:
        # Getting overall alignment line
        overall_align_line = lines[14:15:1]
    else:
        overall_align_line = lines[5:6:1]

    # Getting overall alignment
    if overall_align_line:
        val=overall_align_line[0].split("%")[0]    
    else:
        sys.exit("ERROR: no overall alignment score found.")
        
    return(float(val))


def create_comparison_alignment_dataframe():

    # Checking to see if summary alignment file exists
    if glob.glob(INPUT_DIR+"/*/*/*/*summaryfile.txt"):
        alignment_file_type = "*summaryfile.txt"
    else:
        sys.exit("Error: no alignment file provided (i.e. summaryfile.txt)")

    # Generating list of overall alignment values for HISAT2 and Bowtie2 across all replicates
    h_list = []
    b_list = []
    for alignment_file in glob.glob(INPUT_DIR+"/*/*/*/"+alignment_file_type):    
        print "Reading in: ", alignment_file
        if "bowtie2" in alignment_file:
            b_list.append(extract_alignment_val(alignment_file))
        elif "hisat2" in alignment_file:
            h_list.append(extract_alignment_val(alignment_file))
    return(h_list, b_list)

def create_alignment_dataframe():
    once = []
    overall = []
    reps = []

    # Checking to see if summary alignment file exists
    if glob.glob(INPUT_DIR+"/*/*/*summaryfile.txt"):
        alignment_file_type = "*summaryfile.txt"
    else:
        sys.exit("Error: no alignment file provided (i.e. summaryfile.txt)")

    for alignment_file in glob.glob(INPUT_DIR+"/*/*/"+alignment_file_type):
        print "Reading in: ", alignment_file

        # Openning alignment file and read all lines
        with open(alignment_file) as f:
            lines = f.readlines()

        # If paired-end, check different lines for conc and overall.
        if PAIRED:
            # get only concordant once line and overall alignment line
            once_line = lines[3:4:1]
            overall_align_line = lines[14:15:1]
        else:
            once_line = lines[3:4:1]
            overall_align_line = lines[5:6:1]

        # Getting concordant alignment exactly once
        if once_line:
            in_parenthesis=re.search('\(([^)]+)', once_line[0]).group(1).strip("%")
            once.append(float(in_parenthesis))

        # Getting overall alignment
        if overall_align_line:
            val=overall_align_line[0].split("%")[0]
            overall.append(float(val))
            
            # Get name of replicate and sample files
            name_rep = alignment_file.split("/")[-2]
            reps.append(name_rep)

    # Creating dataframe from alignment information
    align_list = [once, overall]
    df1 = pd.DataFrame(columns=reps, data=align_list)
    return(df1)

def main(argv):
    if not METHOD_COMP and METHOD2_NAME=="NA":

        # Getting alignment dataframe from alignment method 
        df = create_alignment_dataframe()
        
        # Setting dataframe title depending on paired or single-end sequencing
        if PAIRED:
            df['d']=["Concordant reads map exactly once","Overall alignment rate"]
        else:
            df['d']=["Reads map exactly once","Overall alignment rate"]
        df=df.reset_index()
        df.index = df['d']
        del df['d']
        df2= df
        print(df2)

        # Plotting concordant read mapping and overall alignment score for method and save to file
        ax = plt.figure(figsize=(17.7, 8.27))
        colors = ["blue","grey","lightcoral","green","purple","red", "gold", "olive","peru","aqua",
                  "palegreen","mediumorchid","wheat","darkkhaki","m","saddlebrown","lightblue","teal",
                  "orange","plum","deeppink","crimson","yellowgreen","palevioletred","indigo","mediumaquamarine",
                 "blueviolet","y","indianred","lightskyblue","blue","grey","lightcoral","green","purple","red", "gold", "olive","peru","aqua",
                  "palegreen","mediumorchid","wheat","darkkhaki","m","saddlebrown","lightblue","teal",
                  "orange","plum","deeppink","crimson","yellowgreen","palevioletred","indigo","mediumaquamarine",
                 "blueviolet","y","indianred","lightskyblue"]
        columns = list(df2.columns)

        # Generating stripplot
        for j in range(1, len(df2.columns)):
            color = colors[j-1]
            g = sns.stripplot(x=df2.index, y=columns[j], color=color, jitter=0.3, size=14, dodge=True, data=df2, linewidth=2)

        # Adding colors, legend, title, and plot parameters
        elements = [Line2D([0], [0], color=colors[i]) for i in range(len(df2.columns)-1)]
        ax.legend(handles=elements, labels=list(df2.columns)[1:], bbox_to_anchor=(.60,.70))#, prop={'size': 7})
        plt.ylabel("percent")
        plt.ylim(0,100)
        plt.xlabel("measures")
        title = plt.title("Mapping Metrics for "+METHOD_NAME)
            
        # Saving pdf of alignment plot
        plt.tight_layout()
	plt.savefig(INPUT_DIR+"/"+METHOD_NAME+"_plot_alignment.svg")
        return(plt.show())
    
    elif METHOD_COMP:
        hlist, blist = create_comparison_alignment_dataframe()
        # Create dataframe from input
        df = pd.DataFrame(
            {'idx': range(len(hlist)),'HISAT2': hlist,'Bowtie2': blist})

        # Make 1 col for tophat2 and hisat2
        df1 = pd.melt(df, id_vars=['idx'], var_name='methods', value_name= 'counts')
        df1.head()

        fig, ax = plt.subplots(figsize=(11.7, 8.27))
        ax.grid(color='lightgrey', axis='y', linestyle='-', linewidth=2)
        sns.stripplot(x="methods", y="counts", ax = ax, data=df1, jitter=True, linewidth=1)
        title = plt.title("Overall Alignment Scores Between "+METHOD_NAME+" and "+ METHOD2_NAME)
        
        # Calculating the median overall alignment from each methods' alignments
        hmedian = np.median(np.array(hlist))
        bmedian = np.median(np.array(blist))
        print "HISAT2's median is", hmedian
        print "Bowtie2's median is", bmedian 
        if hmedian >= bmedian:
            print("TIMEOR suggests using HISAT2 alignment results.")
        else:        
            print("TIMEOR suggests using Bowtie2 alignment results.")
        
            
        # Saving pdf of alignment plot
        plt.savefig(INPUT_DIR+"/"+METHOD_NAME+"_"+METHOD2_NAME+"_alignment_comparison.pdf")
	#return(plt.show())

if __name__ == "__main__":
    main(sys.argv[1:])
    
