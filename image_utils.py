import rawpy
from PIL import Image
import argparse


def dng_to_tiff(dng_path, tiff_path):
    """Convert a DNG file to a TIFF file."""
    with rawpy.imread(dng_path) as raw:
        rgb_image = raw.postprocess()

    image = Image.fromarray(rgb_image)
    image.save(tiff_path, format="TIFF")
    print(f"Converted {dng_path} to {tiff_path}")


def print_metadata(dng_path):
    """Print metadata of a DNG file."""
    with rawpy.imread(dng_path) as raw:
        metadata = raw.metadata
        print("Metadata:")
        print(metadata)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Image utility functions.")
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Subparser for dng_to_tiff
    parser_dng_to_tiff = subparsers.add_parser("dng_to_tiff", help="Convert a DNG file to a TIFF file.")
    parser_dng_to_tiff.add_argument("dng_path", type=str, help="Path to the input DNG file.")
    parser_dng_to_tiff.add_argument("tiff_path", type=str, help="Path to save the output TIFF file.")

    # Subparser for print_metadata
    parser_print_metadata = subparsers.add_parser("print_metadata", help="Print metadata of a DNG file.")
    parser_print_metadata.add_argument("dng_path", type=str, help="Path to the input DNG file.")

    # Parse arguments
    args = parser.parse_args()

    # Call the appropriate function based on the subcommand
    if args.command == "dng_to_tiff":
        dng_to_tiff(args.dng_path, args.tiff_path)
    elif args.command == "print_metadata":
        print_metadata(args.dng_path)
    else:
        parser.print_help()
