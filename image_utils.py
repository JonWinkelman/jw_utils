import rawpy
from PIL import Image


def dng_to_tiff(dng_path, tiff_path):
    """Open DNG file create array and covert to tiff and save.


    dng_path: path to DNG file
    tiff_path: path where output tiff should be stored
    """
    #
    with rawpy.imread(dng_path) as raw:
        # Convert to an image array
        rgb_image = raw.postprocess()
    
    # Save as a TIFF file using Pillow
    image = Image.fromarray(rgb_image)
    image.save(tiff_path, format="TIFF")
    print(f"Converted {dng_path} to {tiff_path}")
