###

def load_sim_data(block, parent_filepath):
    # Ensure proper formatting
    parent_filepath = os.path.join(parent_filepath, block).replace("\\", "/") 
    if not os.path.exists(os.path.join(parent_filepath, "Object_Data.csv").replace("\\", "/")):
        raise FileNotFoundError(os.path.join(parent_filepath, 'Object_Data.csv').replace("\\", "/") + " doesn't exist.")
    # Load object data
    object_data = pd.read_csv(os.path.join(parent_filepath, 'Object_Data.csv').replace("\\", "/"))
    object_data = object_data[['Analysis Region', 'Algorithm Name', 'Object Id', 'XMin', 'XMax', 'YMin', 'YMax', 'CD8']]
    object_data = object_data[(object_data["Analysis Region"] == "Partition Zone") & (object_data["CD8"] == 1)]
    cell_coords = []
    for idx, row in object_data.iterrows():
        cell_coord = ((row["XMin"] + row["XMax"]) / 2, (row["YMin"] + row["YMax"]) / 2)
        cell_coords.append(cell_coord)
    object_data["Cell Centre"] = cell_coords
    object_data = object_data.drop(['XMin', 'XMax', 'YMin', 'YMax'], axis=1)
    # Load annotation vertices (individual lines of file)
    if not os.path.exists(os.path.join(parent_filepath, 'Annotations.annotations').replace("\\", "/")):
        raise FileNotFoundError(os.path.join(parent_filepath, 'Annotations.annotations').replace("\\", "/") + " doesn't exist.")
    with open(os.path.join(parent_filepath, 'Annotations.annotations').replace("\\", "/"), "r") as f:
        lines = f.readlines()
    # Convert vertices to polygon
    positive_region = []
    layers = []
    polygons = []
    vertices = []
    layer_name = []
    region_status = []
    for line in lines:
        if "<Annotation Name=" in line:
            # Start a new polygon
            layer_name = line.split('Name="')[1].split('" Visible=')[0]
        elif "<Region Type=" in line:
            region_status = line.split('NegativeROA="')[1].split('">')[0] == '0'
            vertices = []
        elif "<V X=" in line:
            # Extract the x,y coords from the line
            x = int(line.split('X="')[1].split('" Y=')[0])
            y = int(line.split('Y="')[1].split('" />')[0])
            vertices.append((x, y))
        elif "</Region>" in line:
            # End current polygon and add it to the list
            if len(vertices) >= 4:
                polygons.append(shapelyPolygon(vertices))
                layers.append(layer_name)
                positive_region.append(region_status)
            else:
                print(f"Small annotation detected for {block}, excluding annotation with {len(vertices)} vertices: {vertices}.")
    # create polygon dataframe to keep track of individual layers on regions
    annotation_layers = pd.DataFrame({"Layer": layers, "Positive Region": positive_region, "Polygon": polygons, "Area": [shapelyPolygon(i).area for i in polygons]})
    annotation_layers = annotation_layers[annotation_layers["Layer"].isin(["Tumour", "Peritumoral Zone", "Partition Zone"])]
    annotation_layers = annotation_layers.sort_values(["Layer", "Positive Region", "Area"], ascending=[False, False, False])
    # Compile relevant annotation layers
    plot_annotations = annotation_layers[((annotation_layers["Layer"] == "Partition Zone") & (annotation_layers["Positive Region"] == 0)) | ((annotation_layers["Layer"] == "Tumour") & (annotation_layers["Positive Region"] == 1))]
    return object_data, plot_annotations
