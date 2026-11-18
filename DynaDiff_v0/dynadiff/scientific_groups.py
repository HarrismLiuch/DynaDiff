"""Display order from the CMasher and Scientific colour maps overview figures."""

# Figure 1 in https://cmasher.readthedocs.io/user/introduction.html#colormap-overview
# orders maps within each category by their lightness profile.
_CMASHER = (
    ("Sequential", (
        "lavender dusk tree cosmic ember emerald lilac nuclear sapphire toxic amber "
        "eclipse ghostlight sepia amethyst apple arctic chroma flamingo freeze gothic "
        "jungle neutral rainforest savanna sunburst swamp torch voltage fall ocean "
        "horizon bubblegum gem pepper neon tropical"
    ).split()),
    ("Diverging", (
        "iceburn redshift seaweed watermelon wildfire guppy pride prinsenvlag fusion "
        "holly viola waterlily"
    ).split()),
    ("Cyclic", "copper_s emergency_s seasons_s infinity_s infinity copper emergency seasons".split()),
)

# Order shown by cmcrameri.show_cmaps() in the project's introductory gallery:
# https://github.com/callumrollo/cmcrameri#extra-instructions
# The gallery omits the categorical S variants; keep them in the same order
# as their sequential parents in a separate category.
_CMCRAMERI_SEQUENTIAL = (
    "batlow batlowW batlowK devon lajolla bamako davos bilbao nuuk oslo grayC hawaii "
    "lapaz tokyo buda acton turku imola glasgow lipari navia"
).split()
_CMCRAMERI = (
    ("Sequential", _CMCRAMERI_SEQUENTIAL),
    ("Diverging", "broc cork vik lisbon tofino berlin roma bam vanimo managua".split()),
    ("Multi-sequential", ("oleron", "bukavu", "fes")),
    ("Cyclic", ("romaO", "bamO", "brocO", "corkO", "vikO")),
    ("Categorical", tuple(name + "S" for name in _CMCRAMERI_SEQUENTIAL)),
)


def scientific_groups(registered_names):
    """Return nonempty provider/category groups of installed forward maps."""
    available = dict.fromkeys(name for name in registered_names if not name.endswith("_r"))
    result = []
    for provider, prefix, categories in (("CMasher", "cmr.", _CMASHER), ("Crameri", "cmc.", _CMCRAMERI)):
        for category, names in categories:
            members = tuple(prefix + name for name in names if prefix + name in available)
            if members:
                result.append((f"{provider} · {category}", members))
                for name in members:
                    del available[name]
        # Keep maps added after the published gallery accessible in registry order.
        extra = tuple(name for name in available if name.startswith(prefix))
        if extra:
            result.append((f"{provider} · Additional", extra))
    return tuple(result)
