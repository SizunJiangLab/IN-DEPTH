import cv2 
import tifffile 
import numpy as np


fusion_image = tifffile.imread('../FINAL-TONSIL_2nd_L24_Scan1.qptiff')

fusion_dapi = fusion_image[0]

fusion_rotate = cv2.rotate(fusion_dapi, cv2.ROTATE_180)

fusion_seg = tifffile.imread('../seg_result/L24_MESMER_mask.tiff')

geomx_image = tifffile.imread('../img_registration/20240312_Tonsil_L24_ReScan2.ome.tiff')

[min_val, max_val] = np.percentile(geomx_image, [0.1, 99.9])

geomx_syto13 = (geomx_image - min_val)/(max_val - min_val)

geomx_syto13 = geomx_syto13 * 255

geomx_syto13 = np.clip(geomx_syto13, 0, 255)

geomx_syto13 = geomx_syto13.astype('uint8')

# Initialize SIFT detector
sift = cv2.SIFT_create(nfeatures = 10000)

# Find keypoints and descriptors
kp1, des1 = sift.detectAndCompute(fusion_rotate, None)
kp2, des2 = sift.detectAndCompute(geomx_syto13, None)

# Initialize matcher
matcher = cv2.BFMatcher()

# Match descriptors
matches = matcher.knnMatch(des1, des2, k=2)

# Apply ratio test
good_matches = []
for m, n in matches:
    if m.distance < 0.7 * n.distance:
        good_matches.append(m)

# Extract matched keypoints
src_pts = np.float32([kp1[m.queryIdx].pt for m in good_matches]).reshape(-1, 1, 2)
dst_pts = np.float32([kp2[m.trainIdx].pt for m in good_matches]).reshape(-1, 1, 2)

# Compute affine transformation
A = cv2.estimateAffinePartial2D(src_pts, dst_pts)[0]

np.savetxt('../img_registration/L24_affine_matrix.txt', A)

A_inverse = np.linalg.inv(np.vstack((A, [0, 0, 1])))

np.savetxt('../img_registration/L24_inverse_affine_matrix.txt', A_inverse)

height, width = fusion_rotate.shape

output_height, output_width  = geomx_syto13.shape



# Create an empty output image
output_image = np.zeros((output_height, output_width), dtype=np.uint8)
output_seg = np.zeros((output_height, output_width), dtype=np.float64)

# Loop over each pixel in the output image
for y_out in range(output_height):
    for x_out in range(output_width):
        # Apply the homography transformation to get the corresponding coordinates in the input image
        homogeneous_coord = np.dot(A_inverse, np.array([[x_out], [y_out], [1]]))
        x_in = int(homogeneous_coord[0] / homogeneous_coord[2])
        y_in = int(homogeneous_coord[1] / homogeneous_coord[2])

        # Check if the transformed coordinate is within the bounds of the input image
        if 0 <= x_in < width and 0 <= y_in < height:
            # Copy the pixel value from the input image to the output image
            output_image[y_out, x_out] = fusion_rotate[y_in, x_in]
            output_seg[y_out,x_out] = fusion_seg[y_in, x_in]


aligned_check = np.stack((geomx_syto13, output_image), axis=0)

with tifffile.TiffWriter('../img_registration/L24_affline_aligned.ome.tiff', bigtiff=True) as tif:
    options = dict(tile = (128, 128))
    tif.write(aligned_check, subifds = 4, **options)
    tif.write(aligned_check[:, ::2, ::2], subfiletype = 1, **options)
    tif.write(aligned_check[:, ::4, ::4], subfiletype = 1, **options)
    tif.write(aligned_check[:, ::8, ::8], subfiletype = 1, **options)
    tif.write(aligned_check[:, ::16, ::16], subfiletype = 1, **options)

tifffile.imwrite('../img_registration/L24_aligned_MESMER_mask.tiff', output_seg)