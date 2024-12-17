import cv2 
import tifffile 
import numpy as np
from tqdm import tqdm
import pandas as pd 
import os 

cores = ["Rochester_1", "Rochester_4", "Rochester_5", "Rochester_6",
                         "Rochester_7", "Rochester_8", "Rochester_9", "Rochester_10",
                         "Rochester_11", "Rochester_12", "Rochester_13", "Rochester_14",
                         "Rochester_15", "Rochester_16", "Rochester_17", "Rochester_18",
                         "Rochester_19", "Rochester_21", "Rochester_22", "Rochester_23",
                         "Rochester_25", "Rochester_TonsilA"]


fusion_image = tifffile.imread('../data/INDEPTH_DLBCL_Left_Scan1.qptiff')[0]

print('Finished loading fusion scan.')

fusion_rotate = cv2.rotate(fusion_image, cv2.ROTATE_180)

print('Finished rotateing fusion scan')

fusion_seg = tifffile.imread('../data/Rochester_MESMER_mask.tiff')

print('Finished loading MESMER mask.')

Rochester_fusion_position = pd.read_csv('../data/Rochester_core_position.csv')

Rochester_fusion_position['x1_rotate'] = fusion_rotate.shape[1] - Rochester_fusion_position['x1'] 
Rochester_fusion_position['x2_rotate'] = fusion_rotate.shape[1] - Rochester_fusion_position['x2'] 
Rochester_fusion_position['y1_rotate'] = fusion_rotate.shape[0] - Rochester_fusion_position['y1'] 
Rochester_fusion_position['y2_rotate'] = fusion_rotate.shape[0] - Rochester_fusion_position['y2'] 

print('Finished loading Fusion core position')

geomx_image = tifffile.imread('../data/20240420_Final_DLBCL_TMA_Rochester.ome.tiff')

[min_val, max_val] = np.percentile(geomx_image, [0.1, 99.9])

geomx_syto13 = (geomx_image - min_val)/(max_val - min_val)

geomx_syto13 = geomx_syto13 * 255

geomx_syto13 = np.clip(geomx_syto13, 0, 255)

geomx_syto13 = geomx_syto13.astype('uint8')

print('Finished loading and processing GeoMx scan.')

Rochester_geomx_core_position = pd.read_csv('../data/Rochester_geomx_position_mask_making.csv')

print('Finished loading GeoMX core position')

for core in tqdm(cores):
    i = core.split('_')[1]
    x1 = Rochester_geomx_core_position[Rochester_geomx_core_position['Core'] == i]['x1'].values[0]
    x2 = Rochester_geomx_core_position[Rochester_geomx_core_position['Core'] == i]['x2'].values[0]
    y1 = Rochester_geomx_core_position[Rochester_geomx_core_position['Core'] == i]['y1'].values[0]
    y2 = Rochester_geomx_core_position[Rochester_geomx_core_position['Core'] == i]['y2'].values[0]
    geomx_core = geomx_syto13[y1:y2, x1:x2]
    x1 = Rochester_fusion_position[Rochester_fusion_position['Core'] == i]['x2_rotate'].values[0]
    x2 = Rochester_fusion_position[Rochester_fusion_position['Core'] == i]['x1_rotate'].values[0]
    y1 = Rochester_fusion_position[Rochester_fusion_position['Core'] == i]['y2_rotate'].values[0]
    y2 = Rochester_fusion_position[Rochester_fusion_position['Core'] == i]['y1_rotate'].values[0]
    fusion_core = fusion_rotate[y1:y2, x1:x2]
    fusion_core_seg = fusion_seg[y1:y2, x1:x2]

    # Initialize SIFT detector
    sift = cv2.SIFT_create(nfeatures = 10000)

    # Find keypoints and descriptors
    kp1, des1 = sift.detectAndCompute(fusion_core, None)
    kp2, des2 = sift.detectAndCompute(geomx_core, None)

    # Initialize matcher
    matcher = cv2.BFMatcher()

    # Match descriptors
    matches = matcher.knnMatch(des1, des2, k=2)

    # Apply ratio test
    good_matches = []
    for m, n in matches:
        if m.distance < 0.5 * n.distance:
            good_matches.append(m)

    print('Finished feature detection and feature matching.')

    # Extract matched keypoints
    src_pts = np.float32([kp1[m.queryIdx].pt for m in good_matches]).reshape(-1, 1, 2)
    dst_pts = np.float32([kp2[m.trainIdx].pt for m in good_matches]).reshape(-1, 1, 2)

    # Compute affine transformation
    A = cv2.estimateAffinePartial2D(src_pts, dst_pts)[0]

    output_path = f"../output/img_registration/fusion_to_geomx/Rochester_mask_making/Rochester_{i}"

    if not os.path.exists(output_path):
        os.makedirs(output_path)

    np.savetxt(f'{output_path}/Rochester_affine_matrix.txt', A)

    A_inverse = np.linalg.inv(np.vstack((A, [0, 0, 1])))

    np.savetxt(f'{output_path}/Rochester_inverse_affine_matrix.txt', A_inverse)

    print('Finished writing transformation matrices.')

    height, width = fusion_core.shape

    output_height, output_width  = geomx_core.shape
    
    # Create an empty output image
    
    output_image = np.zeros((output_height, output_width), dtype=np.uint8)
    output_seg = np.zeros((output_height, output_width), dtype=np.float64)

    # Loop over each pixel in the output image
    for y_out in tqdm(range(output_height)):
        for x_out in range(output_width):
            # Apply the homography transformation to get the corresponding coordinates in the input image
            homogeneous_coord = np.dot(A_inverse, np.array([[x_out], [y_out], [1]]))
            x_in = int(homogeneous_coord[0] / homogeneous_coord[2])
            y_in = int(homogeneous_coord[1] / homogeneous_coord[2])

            # Check if the transformed coordinate is within the bounds of the input image
            if 0 <= x_in < width and 0 <= y_in < height:
                # Copy the pixel value from the input image to the output image
                output_image[y_out, x_out] = fusion_core[y_in, x_in]
                output_seg[y_out,x_out] = fusion_core_seg[y_in, x_in]


    aligned_check = np.stack((geomx_core, output_image), axis=0)

    with tifffile.TiffWriter(f'{output_path}/Rochester_affline_aligned.ome.tiff', bigtiff=True) as tif:
        options = dict(tile = (128, 128))
        tif.write(aligned_check, subifds = 4, **options)
        tif.write(aligned_check[:, ::2, ::2], subfiletype = 1, **options)
        tif.write(aligned_check[:, ::4, ::4], subfiletype = 1, **options)
        tif.write(aligned_check[:, ::8, ::8], subfiletype = 1, **options)
        tif.write(aligned_check[:, ::16, ::16], subfiletype = 1, **options)

    print('Finished writing aligned image.')

    tifffile.imwrite(f'{output_path}/Rochester_aligned_MESMER_mask.tiff', output_seg)

    print('Finished writing aligned MESMER mask.')
        