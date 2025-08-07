import pandas as pd
import os

# Directory containing the individual toxicity files
directory = '/home/sanket/admet-ai-reproduction/All toxicity files'

# Get a list of all the csv files in the directory
all_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.csv') and 'pkcsm_predictions' not in f]

# Create a list to hold the dataframes
df_list = []

# Loop through the files and read them into dataframes
for file in all_files:
    df = pd.read_csv(file, sep='\t')
    df_list.append(df)

# Concatenate all the dataframes into a single dataframe
combined_df = pd.concat(df_list, ignore_index=True)

# Save the combined dataframe to a new csv file
combined_df.to_csv('/home/sanket/admet-ai-reproduction/data/pkcsm_predictions_combined.csv', index=False)

print("Successfully combined all toxicity files into data/pkcsm_predictions_combined.csv")
