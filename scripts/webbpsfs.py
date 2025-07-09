import argparse
from dataclasses import dataclass
from typing import Optional
from astropy.io import fits
from nirphot.utils import generate_filename
from nirphot.images import makewebbpsf
from nirphot.plotting.psfs import save_psf_plots


@dataclass
class Args:
    imgpath: str
    filter: str
    output_directory: str
    output_filename: str
    img_ext: int = 0
    ps: Optional[float] = None
    pa: Optional[float] = None
    ovsam: int = 3
    fov_pix: Optional[int] = None
    fov_as: Optional[float] = None


def parse_arguments() -> Args:
    parser = argparse.ArgumentParser(
        description="Generate and rotate WebbPSF for a given filter and FITS image.",
        epilog="Example: uv run webbpsfs file.fits F200W output psf_out --img_ext 0 --ps 0.03 --ovsam 3 --fov_pix 256 --fov_as 0.1 --pa 0.0",
    )

    parser.add_argument("imgpath", type=str, help="Path to the image FITS file")
    parser.add_argument("filter", type=str, help="Filter to generate the WebbPSF for")
    parser.add_argument(
        "output_directory", type=str, help="Directory to save the output FITS file"
    )
    parser.add_argument(
        "output_filename", type=str, help="Filename for the output FITS file"
    )
    parser.add_argument(
        "--img_ext",
        type=int,
        default=0,
        help="Image extension to use; default is 0 (primary HDU)",
    )
    parser.add_argument(
        "--ps",
        type=float,
        default=None,
        help="Pixel scale for the WebbPSF; default is None",
    )
    parser.add_argument(
        "--pa",
        type=float,
        default=None,
        help="Position angle for the WebbPSF; default is None",
    )
    parser.add_argument(
        "--ovsam",
        type=int,
        default=3,
        help="Oversampling factor for the WebbPSF; default is 3",
    )
    parser.add_argument(
        "--fov_pix",
        type=int,
        default=None,
        help="Field of view in pixels for the WebbPSF; default is None",
    )
    parser.add_argument(
        "--fov_as",
        type=float,
        default=2.0,
        help="Field of view in arcseconds for the WebbPSF; default is 2.0",
    )

    parsed_args = parser.parse_args()
    args_dict = vars(parsed_args)

    # Convert dictionary to Args dataclass
    args = Args(**args_dict)
    return args


def run(args: Args) -> None:
    fp = args.imgpath
    img_ext = args.img_ext
    pixscl = args.ps
    pa = args.pa
    filter = args.filter
    oversample = args.ovsam
    fov_pixels = args.fov_pix
    fov_arcsec = args.fov_as
    output_directory = args.output_directory
    output_filename = args.output_filename

    hdu = fits.open(fp)
    psf = makewebbpsf.generate_webbpsf(
        filter=filter,
        img_hdu=hdu,
        ext=img_ext,
        pixscl=pixscl,
        oversample=oversample,
        fov_pixels=fov_pixels,
        fov_arcsec=fov_arcsec,
    )

    psf_rot = makewebbpsf.rotate_webbpsf(psf, img_hdu=hdu, ext=img_ext, pa=pa)

    full_output_path = generate_filename(output_filename, "fits", output_directory)
    psf_rot.writeto(full_output_path, overwrite=True)

    # save only rotated PSF
    i = fits.PrimaryHDU(
        psf_rot["ROTATED_DET_DIST"].data, header=psf_rot["ROTATED_DET_DIST"].header
    )
    i.header["EXTNAME"] = "ROTATED_DET_DIST"

    full_output_path_rot = generate_filename(
        output_filename + "_rot", "fits", output_directory
    )
    i.writeto(full_output_path_rot, overwrite=True)

    full_output_path_pdf = generate_filename(output_filename, "pdf", output_directory)
    save_psf_plots(psf_rot, path=full_output_path_pdf, filt=filter)

    print(f"File saved as {full_output_path}")


def main():
    args = parse_arguments()
    run(args)
