"""Subset TMA cores related functions."""


def subset_tma_cores(core_map, cell_df, anno_df, core_id, region):
    """Subset cells and annotations to a specific TMA core."""

    cell_filt = cell_df.merge(core_map[["Coords", "Region"]], on="Coords")
    anno_filt = anno_df.merge(core_map[["Coords", "Region"]], on="Coords")

    if core_id is not None:

        cell_filt = cell_df[cell_df["Coords"] == core_id]
        anno_filt = anno_df[anno_df["Coords"] == core_id]

    if region is not None:

        cell_filt = cell_filt[cell_filt["Region"] == region]
        anno_filt = anno_filt[anno_filt["Region"] == region]

    return cell_filt, anno_filt
