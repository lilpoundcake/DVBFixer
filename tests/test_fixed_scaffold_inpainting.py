"""Also runs via stdlib unittest in the isolated Torch runtime (no pytest needed)."""
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

MODULE_PATH = Path(__file__).resolve().parents[1] / "deploy/protpardelle-1c/fixed_scaffold_inpainting.py"
spec = importlib.util.spec_from_file_location("private_inpainting_tests", MODULE_PATH)
policy = importlib.util.module_from_spec(spec)
spec.loader.exec_module(policy)


class MaskTests(unittest.TestCase):
    def test_rigid_placement_ignores_peptide_neighbour_but_not_other_chain(self):
        from scipy.spatial import cKDTree
        points = np.array([[1., 0., 0.], [2., 0., 0.], [9., 0., 0.]])
        tree = cKDTree([[1.1, 0., 0.], [0., 0., 0.]])
        shifts = policy.residue_shifts(points, np.array([1, 1, 3]), tree, [0, 2], 'AABB')
        np.testing.assert_allclose(shifts, [[.35, 0., 0.], [.35, 0., 0.], [0., 0., 0.]])

    def test_overlapping_or_empty_ownership_is_refused(self):
        known = np.zeros((1, 4, 37), dtype=bool)
        known[0, 0, 0] = True
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError):
                policy.FixedScaffold(known, known, Path(tmp) / "overlap")
            with self.assertRaises(ValueError):
                policy.FixedScaffold(known, np.zeros_like(known), Path(tmp) / "empty")


@unittest.skipUnless(importlib.util.find_spec("torch") is not None, "run isolated sampling-runtime unittest for Torch controls")
class TorchTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import torch
        cls.torch = torch

    def fixture(self, output, device="cpu"):
        torch = self.torch
        known = np.zeros((1, 4, 37), dtype=bool)
        known[0, [0, 3], 0] = True
        known[0, [0, 3], 1] = True
        known[0, [0, 3], 36] = True  # input OXT is known, never generated
        generated = np.zeros_like(known)
        generated[0, [1, 2], 0] = True
        generated[0, [1, 2], 1] = True
        controller = policy.FixedScaffold(known, generated, output, expected_steps=3)
        state = torch.arange(4*37*3, dtype=torch.float32, device=device).reshape(1, 4, 37, 3)/10
        coords = torch.ones((1, 2, 37, 3), device=device)*7
        mask = torch.as_tensor(known[:, [0, 3]], device=device)
        return controller, state, coords, mask

    def test_noise_consistency_generated_bits_and_zero_noise_endpoint(self):
        torch = self.torch
        devices = ["cpu"]+(["mps"] if torch.backends.mps.is_available() else [])
        for device in devices:
            with self.subTest(device=device), tempfile.TemporaryDirectory() as tmp:
                c, state, coords, mask = self.fixture(Path(tmp) / "inpainting", device)
                rng = torch.get_rng_state().clone()
                c.initialize(state, [[0, 3]], coords, mask, 100.)
                self.assertTrue(torch.equal(rng, torch.get_rng_state()))
                state = c.noisy(state, 100., -1, stage="initial")
                for i, (sigma, next_sigma) in enumerate(((100., 50.), (50., 2.), (2., 0.))):
                    denoiser_input = c.noisy(state+3, sigma, i)
                    target = c.clean+sigma*c.epsilon
                    self.assertTrue(torch.equal(denoiser_input[c.known], target[c.known]))
                    raw_prediction = denoiser_input+5
                    prediction = c.prediction(raw_prediction, i)
                    self.assertTrue(torch.equal(prediction[c.known], c.clean[c.known]))
                    self.assertTrue(torch.equal(prediction[c.generated], raw_prediction[c.generated]))
                    sc = c.self_condition(raw_prediction, 2., i)
                    self.assertTrue(torch.equal(sc[c.known], (c.clean/2)[c.known]))
                    before = state+17
                    state = c.after_step(before, next_sigma, i, 3)
                    self.assertTrue(torch.equal(before[c.generated], state[c.generated]))
                    self.assertTrue(torch.equal(state[c.known], (c.clean+next_sigma*c.epsilon)[c.known]))
                self.assertTrue(torch.equal(state[c.known], c.clean[c.known]))
                native_mask = torch.as_tensor(c.generated_numpy, device=device)
                out_mask = c.output_mask(native_mask)
                self.assertTrue(bool(out_mask[0, 0, 36]))
                self.assertFalse(bool(out_mask[0, 1, 36]))
                self.assertEqual(c.summary()["accepted_steps"], 3)
                self.assertEqual(c.maximum_post_error, 0.)
                from scripts.audit_cc94_inpainting import audit_trajectory
                (Path(tmp) / "inpainting-summary.json").write_text(json.dumps(c.summary()))
                proof = audit_trajectory(Path(tmp), c.known_numpy, c.generated_numpy,
                                          state[0].detach().cpu().numpy(), steps=3)
                self.assertTrue(proof["native_endpoint_known_bits_exact"])
                self.assertEqual(proof["projection_payloads_checked"], 10)
                path = c.output / "projection-accepted-step-0001.npz"
                with np.load(path) as data:
                    changed = {k: data[k].copy() for k in data.files}
                changed["known_after"][0, 0] += .1
                np.savez_compressed(path, **changed)
                with self.assertRaisesRegex(ValueError, "noise path"):
                    audit_trajectory(Path(tmp), c.known_numpy, c.generated_numpy,
                                     state[0].detach().cpu().numpy(), steps=3)

    def test_incomplete_native_motif_and_bad_step_order_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            c, state, coords, mask = self.fixture(Path(tmp) / "mask")
            mask[0, 0, 36] = False
            with self.assertRaises(ValueError):
                c.initialize(state, [[0, 3]], coords, mask, 100.)
        with tempfile.TemporaryDirectory() as tmp:
            c, state, coords, mask = self.fixture(Path(tmp) / "order")
            c.initialize(state, [[0, 3]], coords, mask, 100.)
            with self.assertRaises(ValueError):
                c.after_step(state, 0., 2, 3)
            with self.assertRaises(ValueError):
                c.output_mask(c.known)

    def test_unsupported_noise_or_sequence_modes_are_not_silently_changed(self):
        with tempfile.TemporaryDirectory() as tmp:
            c, _, _, _ = self.fixture(Path(tmp) / "mode")
            options = dict(num_steps=3, step_scale=1., s_churn=0., jump_steps=False, uniform_steps=True,
                gt_aatype=object(), cc=SimpleNamespace(enabled=True, num_recurrence_steps=1,
                    reconstruction_guidance=SimpleNamespace(enabled=False), replacement_guidance=SimpleNamespace(enabled=False)),
                pd=SimpleNamespace(enabled=False), sidechain_mode=False)
            c.validate_mode(**options)
            options["s_churn"] = 1.
            with self.assertRaises(ValueError):
                c.validate_mode(**options)

    def test_actual_pinned_sample_patch_is_scoped_and_source_files_unchanged(self):
        import inspect

        from protpardelle.core.models import Protpardelle

        torch = self.torch
        model = Protpardelle.__new__(Protpardelle)
        torch.nn.Module.__init__(model)
        original = model.sample.__func__
        source_path = Path(inspect.getsourcefile(original))
        original_bytes = source_path.read_bytes()
        with tempfile.TemporaryDirectory() as tmp:
            controller, _, _, _ = self.fixture(Path(tmp) / "patch")
            with self.assertRaisesRegex(RuntimeError, "synthetic failure"):
                with policy.installed(model, controller):
                    self.assertIsNot(model.sample.__func__, original)
                    self.assertIn("_dvbfixer_inpaint", model.sample.__func__.__globals__)
                    self.assertEqual(source_path.read_bytes(), original_bytes)
                    raise RuntimeError("synthetic failure")
        self.assertIs(model.sample.__func__, original)
        self.assertNotIn("sample", model.__dict__)
        self.assertEqual(source_path.read_bytes(), original_bytes)


if __name__ == "__main__":
    unittest.main()
