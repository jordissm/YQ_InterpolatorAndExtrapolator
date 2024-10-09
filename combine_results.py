import pandas as pd
import os

def combine_results(output_dir):
    all_files = [os.path.join(output_dir, f) for f in os.listdir(output_dir) if os.path.isfile(os.path.join(output_dir, f))]

    for file in all_files:
        print(file)
        df = pd.read_csv(file, sep='\s+')
        YQ_valuesOfInterest = [0.417, 0.458, 0.5]
        df = df[df['YQ'].isin(YQ_valuesOfInterest)]
        df = df.groupby('YQ').agg({'Mean_value': 'mean', 'Min_value': 'min','Max_value': 'max'}).reset_index()
        #    Uncomment the following line if you want to calculate
        #    the central value as the average of the min and max values
        # df['Central_value'] = (df['Min_value'] + df['Max_value']) / 2.0
        #    Uncomment the following line if you want to calculate the
        #    central value as the mean value
        df['Central_value'] = df['Mean_value']
        df['Delta-'] = df['Central_value'] - df['Min_value']
        df['Delta+'] = df['Max_value'] - df['Central_value']
        df = df[['YQ', 'Central_value', 'Delta-', 'Delta+']]
        print(df)

    return 0

if __name__ == "__main__":
    output_directory = "output"
    combined_results_df = combine_results(output_directory)

    # Optionally, save the combined results to a new file
    # combined_results_df.to_csv("combined_results.csv", index=False)