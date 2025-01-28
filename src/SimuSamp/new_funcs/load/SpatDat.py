"""
Loads data exported from HALO into Shapely objects
"""

from SimuSamp.new_funcs.load.load_data import load_data
from SimuSamp.new_funcs.compute.compute_hpfs import compute_hpfs
from SimuSamp.new_funcs.load.subset_cells import subset
from SimuSamp.new_funcs.load.subset_annotations import subset_anno
from SimuSamp.new_funcs.compute.poisson import poisson_cells


class SpatDat:
    def __init__(self, sampleid, parent_filepath, cell_name="CD8"):

        self.sampleid = sampleid
        self.parent_filepath = parent_filepath
        self.cell_name = cell_name

        self.filepath = f"{self.parent_filepath}/{self.sampleid}"

        self.object_data, self.annotation_data = load_data(
            self.filepath, self.cell_name
        )

        self.poisson_cells = {}

    def subset_cells(self, annotation):

        cells = subset(self.object_data, self.annotation_data, annotation)

        return cells

    def poisson_distribution(self, annotation, n_cells=None):

        random_cells = poisson_cells(
            self.object_data, self.annotation_data, annotation, n_cells
        )

        self.poisson_cells[annotation] = random_cells

        return random_cells

    def subset_annotation(self, annotation):

        anno = subset_anno(self.annotation_data, annotation)

        return anno

    def compute_fields(self, width_microns):

        hpfs = compute_hpfs(
            self.sampleid, self.object_data, self.annotation_data, width_microns
        )

        self.hpfs = hpfs
