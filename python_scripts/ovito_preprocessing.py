#  //  Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
#  //  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
#  //  https://www.leibniz-hki.de/en/applied-systems-biology.html
#  //  HKI-Center for Systems Biology of Infection
#  //  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
#  //  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
#  //
#  //  This code is licensed under BSD 2-Clause
#  //  See the LICENSE file provided with this code for the full license.

import numpy as np
import os
import polars as pl

import sys

path = sys.argv[1]

if not os.path.exists(path+"/xyz_files"):
    os.mkdir(path+"/xyz_files")

# loop over all csv files
for file in os.listdir(path):
    if (file == "agent-statistics.csv"):
        print("Creating '.xyz' files for agent-statistics")
        df2 = pl.read_csv(path + "/agent-statistics.csv", sep=";", dtypes={"time":pl.Float64})
        nb_agent = df2.filter(pl.col("time")==0).shape[0]
        count = 0
        for i in list(np.unique(df2["time"])):
            with open(path+"/xyz_files/agent_visualisation_" + str(count)+".xyz", 'w') as f:
                f.write(str(nb_agent))
                f.write('\n')
                f.write('\n')
                df2.filter(pl.col("time") == i).select(["x", "y", "z", "radius", "agent"]).to_pandas().to_csv(f, header=False, sep='\t', index = False)
                count +=1

    elif (file == "molecules.csv"):
        print("Creating '.xyz' files for molecules")
        df = pl.read_csv(path + "/molecules.csv", sep=";", dtypes={"time":pl.Float64, "value":pl.Float64, "x":pl.Float64, "y":pl.Float64, "z":pl.Float64})
        grid_size = df.filter(pl.col("time") == 0).shape[0]
        count = 0
        for i in list(np.unique(df["time"])):
            with open(path+"/xyz_files/molecules_visualisation_" + str(count)+".xyz", 'w') as f:
                f.write(str(int(grid_size)))
                f.write('\n')
                f.write('\n')
                df.filter(pl.col("time") == i).select(["x", "y", "z", "value", "name"]).to_pandas().to_csv(f, header=False, sep='\t', index = False)
                count +=1

    elif (file == "molecules_test.csv"):
        print("Creating '.xyz' files for molecules")
        df = pl.read_csv(path + "/molecules_test.csv", sep=";", dtypes={"time":pl.Float64, "x":pl.Float64, "y":pl.Float64, "z":pl.Float64, "AMP":pl.Float64, "Complex":pl.Float64, "Defensive":pl.Float64})
        grid_size = df.filter(pl.col("time") == 0).shape[0]
        for j in range(5,8):
            count = 0
            for i in list(np.unique(df["time"])):
                with open(path+"/xyz_files/molecules_visualisation_" + str(df.columns[j]) + "_" + str(count)+".xyz", 'w') as f:
                    f.write(str(int(grid_size)))
                    f.write('\n')
                    f.write('\n')
                    df2 = df.filter(pl.col("time") == i)
                    df3 = df2.with_column(pl.lit(df2.columns[j]).alias('name'))
                    df3[:, [2, 3, 4, j, 8]].to_pandas().to_csv(f, header=False, sep='\t', index = False)
                count +=1

    elif(file == "membrane.csv"):
        df = pl.read_csv(path + "/membrane.csv", sep=";", dtypes={"x":pl.Float64, "y":pl.Float64, "z":pl.Float64})
        nb_points = df.shape[0]
        with open(path+"/xyz_files/membrane.xyz", 'w') as f:
            f.write(str(nb_points))
            f.write('\n')
            f.write('\n')
            df.select(["x", "y", "z", "TypeName"]).to_pandas().to_csv(f, header=False, sep='\t', index = False)

    elif(file == "inside.csv"):
        df = pl.read_csv(path + "/inside.csv", sep=";", dtypes={"x":pl.Float64, "y":pl.Float64, "z":pl.Float64})
        nb_points = df.shape[0]
        with open(path+"/xyz_files/inside.xyz", 'w') as f:
            f.write(str(nb_points))
            f.write('\n')
            f.write('\n')
            df.select(["x", "y", "z", "TypeName"]).to_pandas().to_csv(f, header=False, sep='\t', index = False)