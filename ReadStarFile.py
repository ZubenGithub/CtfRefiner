#!/data/zuben/miniconda2/bin/python

# This will be for generic star file reading
# Currently uses pyem
# I should make my own one -- later

## Tasks
## 1. argument parser to read star file
## 2. Find Defocus statistics PER micrograph
## 3. REMOVE bad particles based on defocus  (need to work out what that means)
## 4. SAVE exaclty same star file (same header etc) but just missing those bad particles


from pyem import star
import re

star_file = "run_data_job203_10000lines.star"

def ReadStarFile(star_file):
    """
    Reads a star file and returns a pandas dataframe.
        NB: I need to test the different programs that it can read.
    Args:
        Star file
    Returns:
        Pandas Dataframe
    """
    dataframe = star.parse_star(star_file)
    columns = dataframe.columns
    new_columns = []
    # Change column names to general column names (without numbers)
    # Just parses the header and finds the names
    for column in columns:
        no_number = column.split()[0]
        new_column_name = re.sub('rln','',no_number) #NB: removed _
        new_columns.append(new_column_name)

    # Renames the columns with the neater and nicer names
    dataframe.columns = new_columns
    return(dataframe)


df= ReadStarFile(star_file)
print(df.columns)
