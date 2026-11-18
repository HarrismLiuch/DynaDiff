"""Checks for spot integration, detector geometry and display-only contrast."""

import unittest

import numpy as np

from dynadiff.detector import PROFILES, color_lut, colormap_groups, colormap_names, detector_image, scale_image


class DetectorTests(unittest.TestCase):
    def test_integrated_flux_and_coordinates(self):
        xy = [[0.145, -0.111], [0.145, -0.111]]
        index = np.arange(257)
        for profile in PROFILES:
            with self.subTest(profile=profile):
                image = detector_image(xy, [0.5, 0.2], 1, size=257,
                                       fwhm=0.024, profile=profile)
                flux = image.sum()
                self.assertAlmostEqual(flux, 0.7, places=8)
                self.assertAlmostEqual(index @ image.sum(axis=0) / flux, 128 * (1 + 0.145), places=6)
                self.assertAlmostEqual(index @ image.sum(axis=1) / flux, 128 * (1 + 0.111), places=6)
                np.testing.assert_allclose(image, detector_image([xy[0]], [0.7], 1, size=257,
                                                                 fwhm=0.024, profile=profile),
                                           atol=1e-15)

    def test_fwhm_and_detector_edges(self):
        for profile in PROFILES:
            with self.subTest(profile=profile):
                image = detector_image([[0, 0]], [1], 1, size=257,
                                       fwhm=0.0625, profile=profile)
                # Eight pixels FWHM. Pixel integration slightly rounds the crest.
                self.assertAlmostEqual(image[128, 132] / image[128, 128], 0.5, delta=0.012)
                edge = detector_image([[1, 0]], [1], 1, size=257,
                                      fwhm=0.0625, profile=profile)
                self.assertLess(edge.sum(), 0.65)
                self.assertGreater(edge.sum(), 0.45)
                self.assertLess(edge[128, 0], edge[128, -1] * 1e-8)

    def test_narrow_gaussian_concentrates_in_containing_pixel(self):
        image = detector_image([[0.01, 0]], [1.0], 1, size=65, fwhm=0.001)
        self.assertAlmostEqual(image[32, 32], 1.0, places=10)

    def test_clipping_and_log_preserve_input(self):
        data = np.arange(100, dtype=float).reshape(10, 10)
        data[-1, -1] = 10_000
        original = data.copy()
        linear, lower, upper = scale_image(data, clip_percent=2, logarithmic=False)
        self.assertEqual(lower, 0)
        self.assertAlmostEqual(upper, np.percentile(data, 98))
        self.assertEqual(linear[-1, -1], 1)
        levels, lower, upper = scale_image(np.array([[0., 1., 100., 10_000.]]),
                                          clip_percent=0, logarithmic=True)
        np.testing.assert_allclose(levels, [[0., 0., 0.5, 1.]])
        self.assertEqual((lower, upper), (1, 10_000))
        np.testing.assert_array_equal(data, original)

    def test_all_colormap_families_and_reversal(self):
        names = colormap_names()
        for name in ("gray", "twilight", "cmc.batlow", "cmr.eclipse"):
            with self.subTest(name=name):
                self.assertIn(name, names)
                lut = color_lut(name)
                self.assertEqual(lut.shape, (256, 4))
                np.testing.assert_array_equal(color_lut(name, reverse=True), lut[::-1])

    def test_gallery_order_and_complete_registry(self):
        from matplotlib import colormaps

        groups = colormap_groups()
        names = colormap_names()
        self.assertEqual(groups[0], ("Matplotlib · Perceptually Uniform Sequential",
                                     ("viridis", "plasma", "inferno", "magma", "cividis")))
        self.assertEqual(tuple(name for _, entries in groups for name in entries), names)
        self.assertCountEqual(names, [name for name in colormaps if not name.endswith("_r")])
        self.assertEqual(tuple(dict.fromkeys(label.split(" · ")[0] for label, _ in groups)),
                         ("Matplotlib", "CMasher", "Crameri"))
        self.assertEqual(dict(groups)["CMasher · Sequential"][:3],
                         ("cmr.lavender", "cmr.dusk", "cmr.tree"))
        self.assertEqual(dict(groups)["Crameri · Sequential"][:3],
                         ("cmc.batlow", "cmc.batlowW", "cmc.batlowK"))


if __name__ == "__main__":
    unittest.main()
