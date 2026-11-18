"""Matplotlib's published colormap categories and gallery order.

https://matplotlib.org/stable/gallery/color/colormap_reference.html
"""

_REFERENCE_GROUPS = (
    ("Perceptually Uniform Sequential", ("viridis", "plasma", "inferno", "magma", "cividis")),
    ("Sequential", (
        "Greys", "Purples", "Blues", "Greens", "Oranges", "Reds", "YlOrBr", "YlOrRd",
        "OrRd", "PuRd", "RdPu", "BuPu", "GnBu", "PuBu", "YlGnBu", "PuBuGn", "BuGn", "YlGn",
    )),
    ("Sequential (2)", (
        "gray", "bone", "pink", "spring", "summer", "autumn", "winter", "cool",
        "Wistia", "hot", "afmhot", "gist_heat", "copper",
    )),
    ("Diverging", (
        "PiYG", "PRGn", "BrBG", "PuOr", "RdGy", "RdBu", "RdYlBu", "RdYlGn",
        "Spectral", "coolwarm", "bwr", "seismic", "berlin", "managua", "vanimo",
    )),
    ("Cyclic", ("twilight", "twilight_shifted", "hsv")),
    ("Qualitative", (
        "Pastel1", "Pastel2", "Paired", "Accent", "okabe_ito", "Dark2", "Set1",
        "Set2", "Set3", "tab10", "tab20", "tab20b", "tab20c",
    )),
    ("Miscellaneous", (
        "flag", "prism", "ocean", "gist_earth", "terrain", "gist_stern", "gnuplot",
        "gnuplot2", "CMRmap", "cubehelix", "brg", "gist_rainbow", "rainbow", "jet",
        "turbo", "nipy_spectral", "gist_ncar",
    )),
    # Legacy grayscale names listed below the reference charts, then spelling aliases.
    ("Aliases", ("gist_gray", "gist_yarg", "binary", "grey", "gist_grey", "gist_yerg", "Grays")),
)


def matplotlib_groups(registered_names):
    """Return the available groups, preserving registry order for newly added maps."""
    remaining = dict.fromkeys(name for name in registered_names
                              if not name.endswith("_r") and not name.startswith(("cmr.", "cmc.")))
    groups = []
    for category, gallery_names in _REFERENCE_GROUPS:
        names = tuple(name for name in gallery_names if name in remaining)
        if names:
            groups.append((category, names))
            for name in names:
                del remaining[name]
    if remaining:
        groups.append(("Additional registered maps", tuple(remaining)))
    return tuple(groups)
