#!/usr/bin/python3

# ******************************************************************************
# Copyright 2026, SHINE Technologies, LLC. MIT License: see LICENSE.
# ******************************************************************************

"""Make the geometry plots of the models, from each system's plot command files
mcnp/<system>/plot_*: one PNG per plot file and code in images/<system>/.

    python plot_geometry.py                     every system, both codes
    python plot_geometry.py --system apollo     one system (artemis, apollo)
    python plot_geometry.py --mcnp              MCNP: <system>_mcnp_geom_<plot>.png
    python plot_geometry.py --openmc            OpenMC: <system>_openmc_geom_<plot>.png

The systems are the folders in mcnp/ that hold a deck of their name (mcnp/<system>/<system>.i).
MCNP plots need mcnp6 and Ghostscript (gswin64c on Windows, gs elsewhere) on the PATH, and Pillow:
the MCNP plotter runs in batch mode on each plot file and writes PostScript, which Ghostscript turns
into a PNG, cropped to the plot frame (and, for a rectangular window, shrunk back to its shape from
the square frame MCNP draws it in). OpenMC plots run openmc/<system>/<system>.py --plot (it needs
OpenMC), for the systems that have an OpenMC model. The OpenMC plots also outline the materials in
the flux maps of plot_results.py.
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
MCNP_ROOT = os.path.join(HERE, "mcnp")
IMAGES_DIR = os.path.join(HERE, "images")


def systems():
    """The systems: the folders in mcnp/ that hold a deck of their name."""
    return sorted(
        os.path.basename(d)
        for d in glob.glob(os.path.join(MCNP_ROOT, "*"))
        if os.path.isfile(os.path.join(d, f"{os.path.basename(d)}.i"))
    )


def plot_files(system):
    """Names of the plot command files in mcnp/<system>/ (plot_xy_pit, ...), without the plot_
    prefix."""
    files = glob.glob(os.path.join(MCNP_ROOT, system, "plot_*"))
    return sorted(os.path.basename(f)[len("plot_") :] for f in files if os.path.isfile(f))


def mcnp_commands(system, name):
    """The plot file's commands, without comments, on as few lines as fit in 80 columns, followed by
    scale 1 (a frame around the plot, to crop to) and end (to leave the plotter)."""
    tokens = []
    for line in open(os.path.join(MCNP_ROOT, system, f"plot_{name}")):
        tokens += [line.split("$")[0].split()]
    lines = [""]
    for t in [x for x in tokens if x] + [["scale", "1"]]:
        text = " ".join(t)
        if lines[-1] and len(lines[-1]) + 1 + len(text) > 80:
            lines.append("")
        lines[-1] = (lines[-1] + " " + text).strip()
    return "\n".join(lines) + "\nend\n"


def plot_extent(system, name):
    """The plot's half-widths along its horizontal and vertical axes (its last extent; one value
    for a square window)."""
    extent = None
    for line in open(os.path.join(MCNP_ROOT, system, f"plot_{name}")):
        f = line.split("$")[0].split()
        if f and f[0] == "extent":
            extent = [float(v) for v in f[1:3]]
    return extent * 2 if len(extent) == 1 else extent


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


def plot_mcnp(system, dpi):
    from PIL import Image

    mcnp = shutil.which("mcnp6")
    gs = shutil.which("gswin64c" if platform.system() == "Windows" else "gs")
    if mcnp is None or gs is None:
        raise SystemExit("MCNP plots need mcnp6 and Ghostscript on the PATH")

    # Run in mcnp/<system>/ (the deck reads its source from ../../source/), with every file in a
    # temporary folder
    folder = os.path.join(MCNP_ROOT, system)
    tmp = os.path.join(folder, "geometry_temp")
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp)
    try:
        for name in plot_files(system):
            com = os.path.join("geometry_temp", f"{name}.com")
            with open(os.path.join(folder, com), "w", newline="\n") as w:
                w.write(mcnp_commands(system, name))
            ps = os.path.join("geometry_temp", name)  # MCNP adds .ps
            subprocess.run(
                [mcnp, "ip", f"i={system}.i", f"n={ps}.", "notek", f"com={com}", f"plotm={ps}"],
                cwd=folder,
                check=True,
                stdout=subprocess.DEVNULL,
            )

            # Each command line draws the plot again: keep only the last page of the PostScript
            lines = open(os.path.join(folder, ps + ".ps")).readlines()
            ends = [i for i, x in enumerate(lines) if "stroke showpage" in x]
            start = ends[-2] + 1 if len(ends) > 1 else 0
            last = os.path.join(folder, ps + "_last.ps")
            with open(last, "w") as w:
                w.writelines(lines[start : ends[-1] + 1])

            png = os.path.join(tmp, f"{name}.png")
            subprocess.run([gs, "-q", "-sDEVICE=png16m", f"-r{dpi}", "-o", png, last], check=True)
            image = crop_to_frame(Image.open(png).convert("RGB").rotate(-90, expand=True))
            # MCNP draws every window in a square frame, stretching a rectangular one (extent with
            # two values): shrink the longer side back to the window's shape, averaging over area
            # so thin lines stay
            eh, ev = plot_extent(system, name)
            if eh != ev:
                w, h = image.size
                size = (w, round(h * ev / eh)) if eh > ev else (round(w * eh / ev), h)
                image = image.resize(size, Image.Resampling.BOX)
            os.makedirs(os.path.join(IMAGES_DIR, system), exist_ok=True)
            out = os.path.join(IMAGES_DIR, system, f"{system}_mcnp_geom_{name}.png")
            image.save(out, optimize=True)
            print("Wrote", os.path.relpath(out, HERE))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def plot_openmc(system):
    script = os.path.join(HERE, "openmc", system, f"{system}.py")
    if not os.path.isfile(script):
        print(f"No OpenMC model for {system} ({os.path.relpath(script, HERE)}): no OpenMC plots")
        return
    subprocess.run([sys.executable, script, "--plot"], check=True, stdout=subprocess.DEVNULL)
    for name in plot_files(system):
        print("Wrote", os.path.join("images", system, f"{system}_openmc_geom_{name}.png"))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--system", choices=systems(), help="one system only (default: all)")
    parser.add_argument("--mcnp", action="store_true", help="MCNP plots only")
    parser.add_argument("--openmc", action="store_true", help="OpenMC plots only")
    parser.add_argument(
        "--dpi", type=int, default=135, help="resolution of the MCNP plots (default 135)"
    )
    args = parser.parse_args()
    both = not (args.mcnp or args.openmc)
    for system in [args.system] if args.system else systems():
        if args.mcnp or both:
            plot_mcnp(system, args.dpi)
        if args.openmc or both:
            plot_openmc(system)


if __name__ == "__main__":
    main()
