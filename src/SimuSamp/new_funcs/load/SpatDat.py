"""
Loads data exported from HALO into Shapely objects
"""

from SimuSamp.src.SimuSamp.new_funcs.load.load_data import load_data
from SimuSamp.src.SimuSamp.new_funcs.compute.compute_hpfs import compute_hpfs


class SpatDat:
    def __init__(self, sampleid, parent_filepath, cell_name="CD8"):

        self.sampleid = sampleid
        self.parent_filepath = parent_filepath
        self.cell_name = cell_name

        self.filepath = f"{self.parent_filepath}/{self.sampleid}"

        self.object_data, self.annotation_data = load_data(self.filepath,
                                                           self.cell_name)

    def compute_fields(self, width_microns):

        hpfs = compute_hpfs(self.sampleid,
                            self.object_data,
                            self.annotation_data,
                            width_microns)

        self.hpfs = hpfs
