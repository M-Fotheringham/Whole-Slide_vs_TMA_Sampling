"""
Loads data exported from HALO into Shapely objects
"""

from SimuSamp.functions.load.load_data import load_data
from SimuSamp.functions.load.load_tma_data import load_tma_data
from SimuSamp.functions.compute.compute_hpfs import compute_hpfs
from SimuSamp.functions.load.subset_cells import subset
from SimuSamp.functions.load.subset_annotations import subset_anno
from SimuSamp.functions.load.subset_cores import subset_tma_cores
from SimuSamp.functions.compute.poisson import poisson_cells


class SpatDat:
    def __init__(self, sampleid, parent_filepath, cell_name="CD8", tma=False):

        self.sampleid = sampleid
        self.parent_filepath = parent_filepath
        self.cell_name = cell_name

        self.filepath = f"{self.parent_filepath}/{self.sampleid}"

        self.object_data, self.annotation_data = load_data(
            self.filepath, self.cell_name
        )
        
        if tma:
            self.core_map, self.tma_data, self.tma_annotation = load_tma_data(
                self.sampleid, self.parent_filepath
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

    # def subset_cores(self, core_id=None, region=None):
    #     # Work in progress...
    #     tma_cells, tma_anno = subset_tma_cores(
    #         self.tma_data, self.tma_annotation, core_id, region
    #     )

        # return tma_cells, tma_anno
