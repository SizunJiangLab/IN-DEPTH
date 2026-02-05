import os
import json
import logging
import numpy as np
import tifffile
from tqdm import tqdm
from shapely.geometry import shape, mapping
from shapely import affinity
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Third-party / Custom libraries found in notebook
from pyqupath.geojson import GeojsonProcessor
from pyqupath.tiff import TiffZarrReader, PyramidWriter

# Setup Logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%H:%M:%S'
)
logger = logging.getLogger(__name__)

class GeometryAwarePartitioner:
    """
    Handles the recursive partitioning of a spatial area based on object density
    to create a grid that avoids cutting through objects where possible.
    Enforces INTEGER boundaries for pixel-perfect tiling.
    """
    def __init__(self, objects_data, width, height):
        self.objects = objects_data
        self.width = int(width)   # Ensure canvas is int
        self.height = int(height) # Ensure canvas is int
        self.grid_boxes = []

    def _find_best_split(self, objects):
        """
        Finds a valid split plane (X or Y) that does not intersect any object
        and divides the set as evenly as possible.
        """
        best_split = None
        min_balance_diff = float('inf')

        n = len(objects)
        
        # --- Try X-Axis Split ---
        objects.sort(key=lambda o: o['bounds'][0]) 
        for i in range(n - 1):
            left_group = objects[:i+1]
            right_group = objects[i+1:]
            
            l_max = max(o['bounds'][2] for o in left_group)
            r_min = min(o['bounds'][0] for o in right_group)
            
            if l_max < r_min:
                diff = abs(len(left_group) - len(right_group))
                if diff < min_balance_diff:
                    min_balance_diff = diff
                    # FORCE INTEGER SPLIT
                    split_pos = int((l_max + r_min) / 2) 
                    best_split = ('x', split_pos, i+1)

        # --- Try Y-Axis Split ---
        objects.sort(key=lambda o: o['bounds'][1]) 
        for i in range(n - 1):
            bottom_group = objects[:i+1]
            top_group = objects[i+1:]
            
            b_max = max(o['bounds'][3] for o in bottom_group)
            t_min = min(o['bounds'][1] for o in top_group)
            
            if b_max < t_min:
                diff = abs(len(bottom_group) - len(top_group))
                # Prefer Y split if it balances better, or if X failed
                if diff < min_balance_diff:
                    min_balance_diff = diff
                    # FORCE INTEGER SPLIT
                    split_pos = int((b_max + t_min) / 2)
                    best_split = ('y', split_pos, i+1)

        return best_split

    def _recursive_partition(self, objects_list, x_min, y_min, x_max, y_max):
        # Ensure inputs are integers (sanity check)
        x_min, y_min = int(x_min), int(y_min)
        x_max, y_max = int(x_max), int(y_max)

        # BASE CASE: Only 1 object left
        if len(objects_list) == 1:
            self.grid_boxes.append({
                'box': (x_min, y_min, x_max, y_max),
                'object': objects_list[0]
            })
            return

        split = self._find_best_split(objects_list)
        
        if split:
            axis, pos, idx = split
            pos = int(pos) # Double check it is int
            
            if axis == 'x':
                objects_list.sort(key=lambda o: o['bounds'][0])
                left = objects_list[:idx]
                right = objects_list[idx:]
                self._recursive_partition(left, x_min, y_min, pos, y_max)
                self._recursive_partition(right, pos, y_min, x_max, y_max)
                
            elif axis == 'y':
                objects_list.sort(key=lambda o: o['bounds'][1])
                bottom = objects_list[:idx]
                top = objects_list[idx:]
                self._recursive_partition(bottom, x_min, y_min, x_max, pos)
                self._recursive_partition(top, x_min, pos, x_max, y_max)
        else:
            # FALLBACK: Force split on centroids if no clean cut found
            logger.warning(f"No clean gap found for group size {len(objects_list)}. Forcing integer centroid split.")
            objects_list.sort(key=lambda o: o['centroid'][0])
            mid = len(objects_list) // 2
            
            p1 = objects_list[mid-1]['centroid'][0]
            p2 = objects_list[mid]['centroid'][0]
            
            # FORCE INTEGER SPLIT
            pos = int((p1 + p2) / 2)
            
            self._recursive_partition(objects_list[:mid], x_min, y_min, pos, y_max)
            self._recursive_partition(objects_list[mid:], pos, y_min, x_max, y_max)

    def generate_grid(self):
        self.grid_boxes = []
        self._recursive_partition(self.objects, 0, 0, self.width, self.height)
        return self.grid_boxes

    def export_to_geojson(self, output_path):
        grid_geojson = {"type": "FeatureCollection", "features": []}

        for item in self.grid_boxes:
            # Explicitly cast to int for clean GeoJSON output
            x_min, y_min, x_max, y_max = map(int, item['box'])
            contained_object_id = item['object']['id']
            
            # Create Polygon coordinates (closing loop)
            coordinates = [[
                [x_min, y_min], [x_max, y_min], 
                [x_max, y_max], [x_min, y_max], 
                [x_min, y_min]
            ]]
            
            feature = {
                "type": "Feature",
                "properties": {
                    "contained_object_id": contained_object_id,
                    "type": "grid_partition"
                },
                "geometry": {"type": "Polygon", "coordinates": coordinates}
            }
            grid_geojson['features'].append(feature)

        with open(output_path, 'w') as f:
            json.dump(grid_geojson, f, indent=2)
        
        logger.info(f"Exported {len(grid_geojson['features'])} INTEGER grid boxes to {output_path}")

def load_and_prep_geojson(geojson_path):
    """Loads GeoJSON and calculates bounds/centroids for processing."""
    with open(geojson_path, 'r') as f:
        data = json.load(f)

    objects_data = []
    for feature in data['features']:
        poly = shape(feature['geometry'])
        objects_data.append({
            'id': feature.get('id'),
            'geometry': poly,
            'bounds': poly.bounds, 
            'centroid': (poly.centroid.x, poly.centroid.y)
        })
    return objects_data

def expand_grid_geojson(input_path, output_path, scale=1.10):
    """
    Reads a grid GeoJSON, expands every polygon by a scale factor
    from its center (overlap), and saves to a new file.
    
    NOTE: Expansion naturally creates floats (e.g. 10 * 1.1 = 11.0).
    The resulting overlapping boxes will likely have decimal coordinates,
    but the SRC crop will handle them fine. The critical part is that 
    the DST grid (unexpanded) is integer-based for perfect stitching.
    """
    logger.info(f"Expanding grid from {input_path} by factor {scale}...")
    
    with open(input_path, 'r') as f:
        geojson_data = json.load(f)
    
    expanded_features = []
    
    for i, feature in enumerate(geojson_data['features']):
        original_polygon = shape(feature['geometry'])
        
        # Scaling usually creates floats. We accept this for the expanded/overlapping grid
        # as it is only used for processing individual cores, not for stitching.
        expanded_polygon = affinity.scale(
            original_polygon, xfact=scale, yfact=scale, origin='center'
        )
        
        props = feature.get('properties', {})
        props.update({'expansion_scale': scale, 'is_overlapping': True})
        
        expanded_features.append({
            "type": "Feature",
            "id": feature.get('id', i),
            "properties": props,
            "geometry": mapping(expanded_polygon)
        })
    
    new_geojson = {"type": "FeatureCollection", "features": expanded_features}
    
    with open(output_path, 'w') as f:
        json.dump(new_geojson, f, indent=2)
    
    logger.info(f"Expanded grid saved to {output_path}")

def crop_image_by_geojson(image_path, geojson_path, output_dir):
    """
    Crops a TiffZarrReader image using the polygons defined in a GeoJSON file.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    try:
        tiff = TiffZarrReader.from_ometiff(image_path)
        geojson_processor = GeojsonProcessor.from_path(geojson_path)
        
        logger.info(f"Cropping {os.path.basename(image_path)} into {output_dir}...")
        
        # The crop generator
        # Note: PyQuPath will convert coordinates to slice indices internally.
        # If the input GeoJSON has integers, the slices will be clean.
        crop_gen = geojson_processor.crop_dict_by_polygons(tiff.zimg_dict)
        
        for name, img in tqdm(crop_gen, desc="Tiles"):
            out_file = os.path.join(output_dir, f"{name}.ome.tiff")
            out_dict = {key: item.astype(np.uint16) for key, item in img.items()}
            tiff_writer = PyramidWriter.from_dict(out_dict)
            tiff_writer.export_ometiff_pyramid(out_file)
            
    except Exception as e:
        logger.error(f"Error processing {image_path}: {e}")

def get_image_dimensions(image_path):
    """Reads OME-TIFF to get dimensions (Height, Width)."""
    tiff = TiffZarrReader.from_ometiff(image_path)
    first_key = list(tiff.zimg_dict.keys())[0]
    height, width = tiff.zimg_dict[first_key].shape
    return width, height

def list_files(directory):
    paths = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith('.ome.tiff'):
                paths.append(os.path.join(root, file))
    return paths

# =========================================================
# MAIN PIPELINE RUNNER
# =========================================================

def run_pipeline(config):
    dataset_name = config['name']
    logger.info(f"=== Starting Pipeline for {dataset_name} ===")

    # 1. Get Dimensions
    ref_img = config['ref_image']
    width, height = get_image_dimensions(ref_img)
    logger.info(f"Reference dimensions: {width}x{height}")

    # 2. Generate Integer Grid
    objects_data = load_and_prep_geojson(config['input_geojson'])
    partitioner = GeometryAwarePartitioner(objects_data, width, height)
    partitioner.generate_grid()
    partitioner.export_to_geojson(config['grid_geojson_path'])

    # 3. Process DST (Standard Integer Grid)
    logger.info("Processing DST (Standard Grid)...")
    crop_image_by_geojson(
        image_path=ref_img,
        geojson_path=config['grid_geojson_path'],
        output_dir=config['dst_output_dir']
    )

    # 4. Expand Grid (Overlapping)
    expand_grid_geojson(
        config['grid_geojson_path'], 
        config['expanded_grid_geojson_path'], 
        scale=config['expansion_factor']
    )

    # 5. Process SRC (Expanded Grid)
    logger.info("Processing SRC (Expanded Grid)...")
    src_files = list_files(config['src_input_dir'])
    
    for src_path in tqdm(src_files, desc="Processing Cycles"):
        cycle_num = os.path.basename(src_path).split('.ome.tiff')[0]
        cycle_out_dir = os.path.join(config['src_output_base_dir'], cycle_num)
        
        crop_image_by_geojson(
            image_path=src_path,
            geojson_path=config['expanded_grid_geojson_path'],
            output_dir=cycle_out_dir
        )
    
    logger.info(f"=== Finished Pipeline for {dataset_name} ===")

if __name__ == "__main__":
    
    # Base paths
    BASE_DATA_DIR = "/mnt/nfs/home/huayingqiu/for_steph/INDEPTH_NOV_26/data"
    
    # Configuration for DFCI dataset
    DFCI_CONFIG = {
        'name': 'DFCI',
        'ref_image': f"{BASE_DATA_DIR}/ometiff_lite/DFCI/dst/0.ome.tiff",
        'input_geojson': f"{BASE_DATA_DIR}/02.5_crop_cores/dearrayer_DFCI_correct.geojson",
        'grid_geojson_path': f"{BASE_DATA_DIR}/02.5_crop_cores/DFCI_grid_boxes.geojson",
        'expanded_grid_geojson_path': f"{BASE_DATA_DIR}/02.5_crop_cores/expand_DFCI_grid_boxes.geojson",
        'expansion_factor': 1.10,
        'dst_output_dir': "../data/01.5_crop_grid/DFCI/dst",
        'src_input_dir': f"{BASE_DATA_DIR}/ometiff_lite/DFCI/src",
        'src_output_base_dir': "../data/01.5_crop_grid/DFCI/src"
    }

    # Configuration for Rochester dataset
    ROCHESTER_CONFIG = {
        'name': 'Rochester',
        'ref_image': f"{BASE_DATA_DIR}/ometiff_lite/Rochester/dst/0.ome.tiff",
        'input_geojson': f"{BASE_DATA_DIR}/02.5_crop_cores/dearrayer_Rochester_correct.geojson",
        'grid_geojson_path': f"{BASE_DATA_DIR}/02.5_crop_cores/Rochester_grid_boxes.geojson",
        'expanded_grid_geojson_path': f"{BASE_DATA_DIR}/02.5_crop_cores/expand_Rochester_grid_boxes.geojson",
        'expansion_factor': 1.10,
        'dst_output_dir': "../data/01.5_crop_grid/Rochester/dst",
        'src_input_dir': f"{BASE_DATA_DIR}/ometiff_lite/Rochester/src", # Assuming similar structure
        'src_output_base_dir': "../data/01.5_crop_grid/Rochester/src" # Assuming similar structure
    }

     # Configuration for Tub97.1 dataset
    TUB971_CONFIG = {
        'name': 'Tub971',
        'ref_image': f"{BASE_DATA_DIR}/ometiff_lite/Tub97.1/dst/0.ome.tiff",
        'input_geojson': f"{BASE_DATA_DIR}/02.5_crop_cores/dearrayer_Tub97.1_correct.geojson",
        'grid_geojson_path': f"{BASE_DATA_DIR}/02.5_crop_cores/Tub97.1_grid_boxes.geojson",
        'expanded_grid_geojson_path': f"{BASE_DATA_DIR}/02.5_crop_cores/expand_Tub97.1_grid_boxes.geojson",
        'expansion_factor': 1.10,
        'dst_output_dir': "../data/01.5_crop_grid/Tub97.1/dst",
        'src_input_dir': f"{BASE_DATA_DIR}/ometiff_lite/Tub97.1/src", # Assuming similar structure
        'src_output_base_dir': "../data/01.5_crop_grid/Tub97.1/src" # Assuming similar structure
    }

     # Configuration for Tub97.2 dataset
    TUB972_CONFIG = {
        'name': 'Tub972',
        'ref_image': f"{BASE_DATA_DIR}/ometiff_lite/Tub97.2/dst/0.ome.tiff",
        'input_geojson': f"{BASE_DATA_DIR}/02.5_crop_cores/dearrayer_Tub97.2_correct.geojson",
        'grid_geojson_path': f"{BASE_DATA_DIR}/02.5_crop_cores/Tub97.2_grid_boxes.geojson",
        'expanded_grid_geojson_path': f"{BASE_DATA_DIR}/02.5_crop_cores/expand_Tub97.2_grid_boxes.geojson",
        'expansion_factor': 1.10,
        'dst_output_dir': "../data/01.5_crop_grid/Tub97.2/dst",
        'src_input_dir': f"{BASE_DATA_DIR}/ometiff_lite/Tub97.2/src", # Assuming similar structure
        'src_output_base_dir': "../data/01.5_crop_grid/Tub97.2/src" # Assuming similar structure
    }

    # Run Pipelines
    try:
        # run_pipeline(DFCI_CONFIG)
        # Uncomment to run Rochester when ready
        # run_pipeline(ROCHESTER_CONFIG)
        run_pipeline(TUB971_CONFIG)
        run_pipeline(TUB972_CONFIG)
    except KeyboardInterrupt:
        logger.warning("Pipeline interrupted by user.")
    except Exception as e:
        logger.exception("Fatal error in pipeline.")