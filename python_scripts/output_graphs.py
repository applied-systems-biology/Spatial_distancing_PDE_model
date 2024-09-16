#  //  Copyright by Yann Bachelot, Christoph Saffer, Paul Rudolph
#  //  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
#  //  https://www.leibniz-hki.de/en/applied-systems-biology.html
#  //  HKI-Center for Systems Biology of Infection
#  //  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
#  //  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
#  //
#  //  This code is licensed under BSD 2-Clause
#  //  See the LICENSE file provided with this code for the full license.

import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import numpy as np
import os
import sys
import scipy

path = sys.argv[1] # path to output directory

if "inside_conc.csv" in os.listdir(path):
    df = pd.read_csv(path+"/inside_conc.csv", sep=";")
    df = df.iloc[:, :-1]
    if "molecules_test.csv" in os.listdir(path):
        df3 = pd.read_csv(path+"/molecules_test.csv", sep=";")
        total_AMP = df3.groupby("time")["AMP"].sum()[0]
    else:
        total_AMP = 1
    for mol in np.unique(df["mol_type"]):
        df2 = df.where(df["mol_type"]==mol).dropna()
        df2["Uptake (norm)"] = df2["inside_concentration"]/total_AMP
        fig = px.scatter(df2, x="time", y="Uptake (norm)")
        fig.update_layout(title = " Inside concentration of " + mol + " over time", yaxis_title = "Concentration inside '"+ str(np.unique(df2["agent"]))[2:-2] +"'", xaxis_title = "time")
        with open(path+"/inside_concentration.html", "a+") as f:
            f.write(fig.to_html(full_html=False, include_plotlyjs='cdn', include_mathjax='cdn'))

    df2 = df.where(df["mol_type"]=="AMP").dropna()

    # Fit a linear function through the origin (no intercept)
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(df2['time'], df2['inside_concentration'])
    # Create a scatter plot of your data points
    fig = go.Figure()

    # Add the scatter plot of your data points
    fig.add_trace(go.Scatter(x=df2['time'], y=df2['inside_concentration'], mode='markers', name='Data Points'))

    # Add the fitted line to the plot
    fitted_line = slope * df2['time']
    fig.add_trace(go.Scatter(x=df2['time'], y=fitted_line, mode='lines', name=f'Fitted Line (Slope={slope:.2f})'))

    fig.update_layout(
        title='',
        xaxis_title='Time',
        yaxis_title='Concentration'
    )
    with open(path+"/slope_inside_conc.html", "a+") as f:
        f.write(fig.to_html(full_html=False, include_plotlyjs='cdn', include_mathjax='cdn'))

    with open(path+'/slope.txt', 'w') as file:
        file.write(f'Slope: {slope:.2f}')
        file.write(f'R2: {r_value**2:.2f}')

if "spatial_dist_center.csv" in os.listdir(path):
    df = pd.read_csv(path + "/spatial_dist_center.csv", sep=";")
    if "inside_conc.csv" in os.listdir(path):
        df_2 = (pd.read_csv(path+"/inside_conc.csv", sep=";")
                .iloc[:, :-1]
                .where(lambda y: y["mol_type"] == "AMP")
                .dropna()
                )
        merged_df = pd.merge(df, df2, on="time", how="left")
        merged_df.loc[(merged_df["dist_to_cell"] >= 3.4) | (merged_df["inside_concentration"].isnull()), "inside_concentration"] = 0
        fig = px.line(merged_df, x="dist_to_cell",
                      y= ["AMP", "Complex", "Defensive", "inside_concentration"],
                      animation_frame= "time",
                      line_shape="spline",
                      markers=True,
                      color_discrete_map={
                          "AMP": "blue",
                          "Complex": "green",
                          "Defensive": "red",
                          "inside_concentration":"lightblue"
                      })
    else:
        fig = px.line(df, x="dist_to_cell",
                  y= ["AMP", "Complex", "Defensive"],
                  animation_frame= "time",
                  line_shape="spline",
                  markers=True,
                  color_discrete_map={
                      "AMP": "blue",
                      "Complex": "green",
                      "Defensive": "red"
                  })
    fig.add_vline(x=3.5, line_width=3, line_dash="dash", line_color="orange", name = "cell membrane")
    fig.add_vrect(x0=-1.0, x1=3.5, line_width=0, line_dash = "dash", fillcolor="orange", opacity=0.2)
    fig.add_annotation(x = 1.5, y = 1/10 * np.max(df["AMP"]), text = "Pathogenic cell", showarrow=False)
    fig.update_layout(xaxis_title = "Distance to center of cell", yaxis_title = r"$Concentration (\mu m^{-3})$", legend_title = "Molecules")
    with open(path+"/spatial_dist_center.html", "a+") as f:
        f.write(fig.to_html(full_html=False, include_plotlyjs='cdn', include_mathjax='cdn', auto_play=False))

if "spatial_dist_center_updated.csv" in os.listdir(path):
    df = pd.read_csv(path + "/spatial_dist_center_updated.csv", sep=";")
    columns_to_replace = df.columns[df.columns != 'time']

    df[columns_to_replace] = df[columns_to_replace].replace(0, np.nan)

    fig = px.line(df, x="dist_to_cell",
                  y= ["AMP", "Complex", "Defensive", "AMP_uptaken"],
                  animation_frame= "time",
                  line_shape="spline",
                  markers=True,
                  color_discrete_map={
                      "AMP": "blue",
                      "Complex": "green",
                      "Defensive": "red",
                      "AMP_uptaken": "lightblue"
                  })
    fig.add_vline(x=3.5, line_width=3, line_dash="dash", line_color="orange", name = "cell membrane")
    fig.add_vrect(x0=-1.0, x1=3.5, line_width=0, line_dash = "dash", fillcolor="orange", opacity=0.1)
    fig.add_annotation(x = 1.5, y = 1/10 * np.max(df["AMP"]), text = "Pathogenic cell", showarrow=False)
    fig.update_layout(xaxis_title = "Distance to center of cell", yaxis_title = r"$Concentration (\mu m^{-3})$", legend_title = "Molecules")
    with open(path+"/spatial_dist_center_updated.html", "a+") as f:
        f.write(fig.to_html(full_html=False, include_plotlyjs='cdn', include_mathjax='cdn', auto_play=False))


if "dist_2_cells.csv" in os.listdir(path):
    df = pd.read_csv(path + "/dist_2_cells.csv", sep=";")
    radius = np.unique(df["radius"])[0]
    center = np.abs(np.unique(df["center_cell"])[0])

    condition = (((df['x'] > center - radius) & (df['x'] < center + radius))
                 | ((df['x'] > -center - radius) & (df['x'] < -center + radius)))
    df.loc[condition, ["AMP", "Complex", "Defensive"]] = np.nan

    fig = px.line(df, x="x",
                  y=["AMP", "Complex", "Defensive"],
                  animation_frame="time",
                  line_shape="spline",
                  markers=True,
                  color_discrete_map={
                      "AMP": "blue",
                      "Complex": "green",
                      "Defensive": "red"
                  })

    fig.add_vline(x=center + radius, line_width=3, line_dash="dash", line_color="orange", name = "cell membrane")
    fig.add_vline(x=center - radius, line_width=3, line_dash="dash", line_color="orange", name = "cell membrane")
    fig.add_vline(x=-center + radius, line_width=3, line_dash="dash", line_color="orange", name = "cell membrane")
    fig.add_vline(x=-center - radius, line_width=3, line_dash="dash", line_color="orange", name = "cell membrane")

    fig.add_vrect(x0=center - radius, x1=center + radius, line_width=0, line_dash = "dash", fillcolor="orange", opacity=0.2)
    fig.add_vrect(x0=-center - radius, x1=-center + radius, line_width=0, line_dash = "dash", fillcolor="orange", opacity=0.2)
    fig.add_annotation(x = center + 1.5, y = 1/10 * np.max(df["AMP"]), text = "Cell", showarrow=False)
    fig.add_annotation(x = -center - 1.5, y = 1/10 * np.max(df["AMP"]), text = "Cell", showarrow=False)
    fig.update_layout(xaxis_title = r"$x (\mu m)$", yaxis_title = r"$Concentration (\mu m^{-3})$", legend_title = "Molecules")
    with open(path+"/dist_2_cells.html", "a+") as f:
        f.write(fig.to_html(full_html=False, include_plotlyjs='cdn', include_mathjax='cdn', auto_play=False))

if "molecules.csv" in os.listdir(path):
    df = pd.read_csv(path+"/molecules.csv", sep=";")
    for mol in np.unique(df["name"]):
        df2 = df.where(df["name"]==mol)
        df_mean = df2.groupby(["time"]).mean()
        time =np.unique(df["time"].to_list())
        fig = px.scatter(x=time, y=df_mean["value"])
        fig.update_layout(title = "Concentration of " + mol + " over time", yaxis_title = "Average concentration", xaxis_title = "time")
        with open(path+"/molecules.html", "a+") as f:
            f.write(fig.to_html(full_html=False, include_plotlyjs='cdn', include_mathjax='cdn'))

if "molecules_test.csv" in os.listdir(path):
    df = pd.read_csv(path+"/molecules_test.csv", sep=";")
    df_mean = df.groupby(["time"]).mean()
    time = np.unique(df["time"].to_list())

    for i in range(4, 7):
        fig = px.scatter(x=time, y=df_mean.iloc[:, i])
        fig.update_layout(title = "Relative concentration of " + list(df.columns)[i+1] + " over time", yaxis_title = "Average concentration", xaxis_title = "time")
        with open(path+"/molecules_test.html", "a+") as f:
            f.write(fig.to_html(full_html=False, include_plotlyjs='cdn', include_mathjax='cdn'))

    if "membrane.csv" in os.listdir(path):
        df2 = pd.read_csv(path + "/membrane.csv", sep=";")
        z = np.unique(df2["z"])[int(len(np.unique(df2["z"]))/2)]
        membrane_x = df2[df2["z"] == z]["x"].to_list()
        membrane_y = df2[df2["z"] == z]["y"].to_list()
        for mol in df.columns[5:8]:
            fig = px.scatter(df[df["z"]== z], x="x",
                             y = "y", color = mol, animation_frame="time", width=800, height=800,range_color=(np.min(df[df["z"]== z][mol]),np.max(df[df["z"]== z][mol])))
            fig.add_trace(go.Scatter(x=membrane_x, y=membrane_y, marker=dict(size=2, symbol="square", color="green"), mode="markers", showlegend=False))
            fig.update_traces(marker=dict(size=20, symbol="square"),
                              selector=dict(mode='markers'))
            fig.update_layout(title = "Evolution of relative concentration of " +mol+" in 2D [z=" + str(z)+"]")
            with open(path+"/2D_visualization.html", "a+") as f:
                f.write(fig.to_html(full_html=False, include_plotlyjs='cdn', include_mathjax='cdn', auto_play=False))

    df_sum = df.groupby(["time"]).sum()
    time = np.unique(df["time"].to_list())
    df_sum["time"] = time
    fig = px.line(df_sum, x="time", y=["AMP", "Defensive", "Complex"], line_shape="linear")
    fig.update_layout(title = "Population over time", yaxis_title = "Pop", xaxis_title = "time")
    with open(path+"/population.html", "w") as f:
        f.write(fig.to_html(full_html=False, include_plotlyjs='cdn', include_mathjax='cdn'))

if "uptake.csv" in os.listdir(path):
    df = pd.read_csv(path + "/uptake.csv", sep = ";")
    df2 = df[df["molecule"] == "AMP"]
    fig = px.line(df2, x="time", y="uptake", line_shape="spline", markers=True)
    fig.update_layout(yaxis_title = "Uptake", xaxis_title = "time")
    with open(path+"/uptake.html", "a+") as f:
        f.write(fig.to_html(full_html=False, include_plotlyjs='cdn', include_mathjax='cdn'))