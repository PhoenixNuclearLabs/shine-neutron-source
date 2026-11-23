#!/usr/bin/python3

# ******************************************************************************
# Copyright 2026, SHINE Technologies. All rights reserved.
# ******************************************************************************

"""Make the geometry plots of the Artemis models, from the plot command files mcnp/plot_*: one
PNG per plot file and code in images/.

    python plot_geometry.py            both
    python plot_geometry.py --mcnp     MCNP: images/artemis_mcnp_geom_<plot>.png
    python plot_geometry.py --openmc   OpenMC: images/artemis_openmc_geom_<plot>.png

MCNP plots need mcnp6 and Ghostscript (gswin64c on Windows, gs elsewhere) on the PATH, and Pillow:
the MCNP plotter runs in batch mode on each plot file and writes PostScript, which Ghostscript turns
into a PNG, cropped to the plot frame. OpenMC plots run openmc/artemis.py --plot (it needs OpenMC).
The OpenMC plots also outline the materials in the flux maps of plot_results.py.
"""

import argparse
import glob
import os
import platform
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MCNP_DIR = os.path.join(HERE, "mcnp")
IMAGES_DIR = os.path.join(HERE, "images")


def plot_files():
    """Names of the plot command files in mcnp/ (plot_xy_pit, ...), without the plot_ prefix."""
    files = [f for f in glob.glob(os.path.join(MCNP_DIR, "plot_*")) if os.path.isfile(f)]
    return sorted(os.path.basename(f)[len("plot_") :] for f in files)


def mcnp_commands(name):
    """The plot file's commands, without comments, on as few lines as fit in 80 columns, followed by
    scale 1 (a frame around the plot, to crop to) and end (to leave the plotter)."""
    tokens = []
    for line in open(os.path.join(MCNP_DIR, f"plot_{name}")):
        tokens += [line.split("$")[0].split()]
    lines = [""]
    for t in [x for x in tokens if x] + [["scale", "1"]]:
        text = " ".join(t)
        if lines[-1] and len(lines[-1]) + 1 + len(text) > 80:
            lines.append("")
        lines[-1] = (lines[-1] + " " + text).strip()
    return "\n".join(lines) + "\nend\n"


def crop_to_frame(image):
    """Crop a plot to the inside of its black frame (drawn by scale 1)."""
    black = np.sum(np.array(image)[:, :, :3], axis=2) == 0
    cols = np.where(black.sum(axis=0) >= 0.75 * black.sum(axis=0).max())[0]
    rows = np.where(black.sum(axis=1) >= 0.75 * black.sum(axis=1).max())[0]
    # The frame lines may be more than one pixel wide: skip to their inner edges
    while cols[1] - cols[0] == 1:
        cols = cols[1:]
    while cols[-1] - cols[-2] == 1:
        cols = cols[:-1]
    while rows[1] - rows[0] == 1:
        rows = rows[1:]
    while rows[-1] - rows[-2] == 1:
        rows = rows[:-1]
    return image.crop((cols[0] + 1, rows[0] + 1, cols[-1], rows[-1]))


def plot_mcnp(dpi):
    from PIL import Image

    mcnp = shutil.which("mcnp6")
    gs = shutil.which("gswin64c" if platform.system() == "Windows" else "gs")
    if mcnp is None or gs is None:
        raise SystemExit("MCNP plots need mcnp6 and Ghostscript on the PATH")

    # Run in mcnp/ (the deck reads its source from ../source/), with every file in a temporary
    # folder
    tmp = os.path.join(MCNP_DIR, "geometry_temp")
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp)
    try:
        for name in plot_files():
            com = os.path.join("geometry_temp", f"{name}.com")
            with open(os.path.join(MCNP_DIR, com), "w", newline="\n") as w:
                w.write(mcnp_commands(name))
            ps = os.path.join("geometry_temp", name)  # MCNP adds .ps
            subprocess.run(
                [mcnp, "ip", "i=artemis.i", f"n={ps}.", "notek", f"com={com}", f"plotm={ps}"],
                cwd=MCNP_DIR,
                check=True,
                stdout=subprocess.DEVNULL,
            )

            # Each command line draws the plot again: keep only the last page of the PostScript
            lines = open(os.path.join(MCNP_DIR, ps + ".ps")).readlines()
            ends = [i for i, x in enumerate(lines) if "stroke showpage" in x]
            start = ends[-2] + 1 if len(ends) > 1 else 0
            last = os.path.join(MCNP_DIR, ps + "_last.ps")
            with open(last, "w") as w:
                w.writelines(lines[start : ends[-1] + 1])

            png = os.path.join(tmp, f"{name}.png")
            subprocess.run([gs, "-q", "-sDEVICE=png16m", f"-r{dpi}", "-o", png, last], check=True)
            image = crop_to_frame(Image.open(png).convert("RGB").rotate(-90, expand=True))
            out = os.path.join(IMAGES_DIR, f"artemis_mcnp_geom_{name}.png")
            image.save(out, optimize=True)
            print("Wrote", os.path.relpath(out, HERE))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def plot_openmc():
    script = os.path.join(HERE, "openmc", "artemis.py")
    subprocess.run([sys.executable, script, "--plot"], check=True, stdout=subprocess.DEVNULL)
    for name in plot_files():
        print("Wrote", os.path.join("images", f"artemis_openmc_geom_{name}.png"))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--mcnp", action="store_true", help="MCNP plots only")
    parser.add_argument("--openmc", action="store_true", help="OpenMC plots only")
    parser.add_argument(
        "--dpi", type=int, default=135, help="resolution of the MCNP plots (default 135)"
    )
    args = parser.parse_args()
    both = not (args.mcnp or args.openmc)
    if args.mcnp or both:
        plot_mcnp(args.dpi)
    if args.openmc or both:
        plot_openmc()


if __name__ == "__main__":
    main()
