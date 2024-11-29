from load_data import load_data


class load_sim_data:
    def __init__(self, block, parent_filepath):
        """
        Args:
            block: This is the unique tissue block identifier.
            parent_filepath: Where the data is stored. The object data and
                            annotations file are in folders for each tissue block, simply named block.
        """
    
        self.block = block
        self.parent_filepath = parent_filepath

        load_data(self.block, self.parent_filepath)
    
    
    # Methods derived from functions
