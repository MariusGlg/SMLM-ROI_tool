""" @author: Marius Glogger
Research Group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.
Allows user to draw a polygon (ROI) on a SMLM dataset (HDF5 files: dbscan or dbcluster). Localizations within the ROI are
either kept (if Mask is set to False in the config.ini file) or deleted (imaged is masked, Mask == True).
Modified dataset is saved as a new HDF5-file (name_Mask or name_ROI).
Coordinates defining the mask/ROI are saved in a xlsm-file."""


from mpl_point_clicker import clicker
from mpl_interactions import zoom_factory, panhandler
import h5py
import numpy as np
import matplotlib.pyplot as plt
import yaml as _yaml
from numba import njit
from configparser import ConfigParser
import pandas as pd
import os

"""
    User information:
    Left-click to place a point
    Right click to remove the nearest point
    Left click on legend to change classes
    Get positions with .get_positions()
    Scroll to zoom
    Middle-click to drag and move
    Close window to save file
"""


#  config.ini file
config = ConfigParser()
file = "config.ini"
config.read(file)
#  config file parameter
# load hdf5 files
path = config["INPUT_FILES"]["path"]
#  extract file names
file_name = config["INPUT_FILES"]["file_name"]
Mask = config["INPUT_FILES"]["Mask"]

# column header

header_dbscan = ["frame", "x", "y", "photons", "sx", "sy", "bg", "lpx", "lpy", "ellipticity", "net_gradient", "len",
                 "n", "photon_rate", "group"]  # column names of dbscan file (len = 11)
header_dbscluster = ["groups", "convex_hull", "area", "mean_frame", "com_x", "com_y", "std_frame",
                 "std_x", "std_y", "n"]  # column header of dbscluster file (len = 10)

x_coor = ()
y_coor = ()
if Mask == "True":
    appendix = "_Mask"
else:
    appendix = "_ROI"



class LoadHDF5(object):
    """ loads .hdf5 files from path.
        :return: lists containing individual dbscan_cluster information."""

    def __init__(self, path):  # path to data
        self.path = path

    def load(self):
        """load .hdf5_file"""
        with h5py.File(self.path, "r") as locs_file:
            key = list(locs_file.keys())[0]  # get key name
            locs = locs_file[str(key)][...]
        data_pd = pd.DataFrame(locs)
        return data_pd


class Save(object):
    """ Saves new .hdf5 files, corresponding .yaml file and .xlsm containing ROI coordinates"""

    def __init__(self, p_coloc, path, filename, polygon):
        self.p_coloc = p_coloc
        self.path = path
        self.filename = filename
        self.polygon = polygon

    def save_pd_to_hdf5(self):
        p1_dbscan_filtered = self.p_coloc.values.tolist()  # convert to list of lists
        p1_dbscan_filtered = [tuple(x) for x in p1_dbscan_filtered]  # convert to list of tuples
        name = self.filename + "_" + appendix +".hdf5"

        if x_coor == "com_x":
            #generates picasso-readable h5py dbcluster files
            with h5py.File(os.path.join(self.path, str(name)), "w") as locs_file:
                # define column types for locs file
                formats = ([(np.float32)] * (len(header_dbscluster)))
                formats[0] = (int)
                formats[9] = (np.uint32)
                ds_dt = np.dtype({'names': header_dbscluster, 'formats': formats})
                locs_file.create_dataset("locs", data=p1_dbscan_filtered, dtype=ds_dt)
            #define column times for dbscan file
        else:
            # generates picasso-readable h5py dbscan files
            with h5py.File(os.path.join(self.path, str(name)), "w") as locs_file:
                # define column types
                formats = ([(np.float32)] * (len(header_dbscan)))
                formats[0] = (np.uint32)
                formats[11] = (np.uint32)
                formats[12] = (np.uint32)
                formats[14] = (int)
                ds_dt = np.dtype({'names': header_dbscan, 'formats': formats})
                locs_file.create_dataset("locs", data=p1_dbscan_filtered, dtype=ds_dt)

    def save_yaml(self):
        name = self.filename + "_" + appendix + ".yaml"
        content = []
        if x_coor == "x":
            # save yaml file to reopen modified dbscan file
            with open(os.path.join(self.path, self.filename + ".yaml"), 'r') as yaml_file:
                text = _yaml.load_all(yaml_file, _yaml.FullLoader)
                with open(os.path.join(self.path, name), 'w') as outfile:
                    for doc in text:
                        content.append(doc)
                    _yaml.dump_all(content, outfile)
        else:
            # save empty yaml file to reopen dbcluster files in picasso
            with open(self.filename + ".yaml", 'w') as yaml_file:
                _yaml.dump_all(content, yaml_file)

    def save_mask(self):
        name = self.filename + appendix +".xlsx"
        df = pd.DataFrame(polygon, columns=["x", "y"])
        writer = pd.ExcelWriter(os.path.join(self.path, name), engine='xlsxwriter')
        df.to_excel(writer, sheet_name=appendix, index=False)
        writer.save()
        pass

    def main(self):
        # Check for datatype and save as hdf5 file
        self.save_pd_to_hdf5()
        self.save_yaml()
        self.save_mask()


# Load dbscan files
HDF5_file = LoadHDF5(os.path.join(path, (file_name + ".hdf5")))  # load dbscan files
HDF5_file_pd = HDF5_file.load()

# check if dbscan or dbscluster file is loaded
if "x" in HDF5_file_pd:
    x_coor = "x"
    y_coor = "y"
else:
    x_coor = "com_x"
    y_coor = "com_y"

fig, ax = plt.subplots(constrained_layout=True)
ax.scatter(x=HDF5_file_pd[x_coor], y=HDF5_file_pd[y_coor], s=0.2, color="black", label="locs")
plt.text(0, 0, "left: place point")
plt.text(0, 10, "right: remove point")
plt.text(0, 20, "scroll to zoom")
plt.text(0, 30, "middle click: drag&move")
plt.text(0, 40, "close to save")


# add zooming and middle click to pan
zoom_factory(ax)
ph = panhandler(fig, button=2)
klicker = clicker(ax, ["mask"], markers=["x"], colors=["blue"], **{"linestyle": "--"})
plt.show()

polygon_dict = klicker.get_positions()  # dict containing polygon vertices
polygon = polygon_dict['mask']  # array containing polygon vertices


@njit(nopython=True)
def check_locs_in_ROI(x, y, polygon):
    """checks if point is inside polygon. Returns list of booleans."""
    n = len(polygon)
    inside = False
    xints = 0.0
    p1x, p1y = polygon[0]
    for i in range(n + 1):
        p2x, p2y = polygon[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xints = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x, p1y = p2x, p2y
    return inside


points = np.stack((HDF5_file_pd[x_coor].to_numpy(), HDF5_file_pd[y_coor].to_numpy()), axis=1)
points_inside = [check_locs_in_ROI(point[0], point[1], polygon) for point in points]
if Mask == "True":
    points_outside = [not elem for elem in points_inside]  # reverse to get locs out of the ROI (Mask the region)
    # filter the original data using boolean mask -> remove localization within ROI
    p1_dbscan_pd_filter = HDF5_file_pd.loc[points_outside, :]
if Mask == "False":
    p1_dbscan_pd_filter = HDF5_file_pd.loc[points_inside, :]  # filter the original data using boolean mask

# Plot the data
fig, ax = plt.subplots(constrained_layout=True)
ax.scatter(x=p1_dbscan_pd_filter[x_coor], y=p1_dbscan_pd_filter[y_coor], s=0.2, color="blue", label="filter")
ax.plot(polygon[:, 0], polygon[:, 1], '--ko')
ax.plot((polygon[0, 0], polygon[-1, 0]), (polygon[0, 1], polygon[-1, 1]), "--ko")
plt.show()

# Save coloc data p1 and p2 as dbscan
p1_dbcluster_save = Save(p1_dbscan_pd_filter, path, file_name, polygon)
p1_dbcluster_save.main()
