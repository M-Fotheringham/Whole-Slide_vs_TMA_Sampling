

def subset_anno(anno_df, annotation):

    anno = anno_df[anno_df["layer"] == annotation]["geometry"].iloc[0]

    return anno
