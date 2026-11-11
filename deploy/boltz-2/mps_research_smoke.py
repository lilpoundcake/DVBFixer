"""Separate fork-specific post-update reinjection pilot; no public trace claims."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from dataclasses import asdict
from pathlib import Path

import numpy as np
import torch
from atom_mapping import model_atom_axis, partition_requested_atoms, validate_feature_axis
from mps_research_reinjection import MPSFixedAtomCallback

from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
from dvbfixer.model.diffusion.pdb_materialize import materialize_candidate_pdb


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def run(args):
    from boltz.data.module.inferencev2 import (
        Boltz2InferenceDataModule,
        PredictionDataset,
        load_input,
    )
    from boltz.data.tokenize.boltz2 import Boltz2Tokenizer
    from boltz.data.types import Manifest
    from boltz.main import (
        Boltz2DiffusionParams,
        BoltzSteeringParams,
        MSAModuleArgs,
        PairformerArgsV2,
        process_inputs,
    )
    from boltz.model.models.boltz2 import Boltz2
    from pytorch_lightning import Trainer, seed_everything
    from rdkit import Chem

    if not torch.backends.mps.is_available():
        raise RuntimeError('native MPS required')
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=args.engine, text=True).strip()
    if revision != '1431d573c3a02a2cb6b5c3b775e7128c7a98d539':
        raise ValueError('private hook requires its reviewed fork revision')
    if digest(args.cache / 'boltz2_conf.ckpt') != '090e82ac8c92f5e943fa1b39e7410a44027bea7243c0bbb3caa67a77fc1428e1':
        raise ValueError('checkpoint digest mismatch')
    request = DiffusionRequest.from_json(args.request.read_text())
    if request.seeds != (7,) or request.candidate_count != 1:
        raise ValueError('this prospective fixation pilot is one seed-7 sample')
    source_path = args.request.absolute().parent / request.normalized_pdb.path
    if digest(source_path) != request.normalized_pdb.sha256:
        raise ValueError('source digest mismatch')
    args.output.mkdir(parents=True, exist_ok=False)
    seed_everything(7, workers=True)
    Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)
    processed_root = args.output / 'boltz'
    process_inputs(data=[args.input.absolute()], out_dir=processed_root,
        ccd_path=args.cache / 'ccd.pkl', mol_dir=args.cache / 'mols',
        msa_server_url='https://api.colabfold.com', msa_pairing_strategy='greedy',
        boltz2=True, preprocessing_threads=1, use_msa_server=False)
    processed = processed_root / 'processed'
    manifest = Manifest.load(processed / 'manifest.json')
    if len(manifest.records) != 1:
        raise ValueError('expected one processed target')
    record = manifest.records[0]
    paths = [processed / p for p in ('structures', 'msa', 'constraints', 'templates', 'mols')]
    tokenized = Boltz2Tokenizer().tokenize(load_input(record, *paths))
    axis, names, token_indices = model_atom_axis(tokenized, request,
        chain_name_by_asym_id={c.chain_id: c.chain_name for c in record.chains})
    represented, passthrough = partition_requested_atoms(request, axis)
    if passthrough:
        raise ValueError('this private pilot requires every observation on the model axis')
    dataset = PredictionDataset(manifest, paths[0], paths[1], args.cache / 'mols', *paths[2:])
    features = dataset[0]
    numbers = validate_feature_axis(axis, token_indices, features)
    fixed = {}
    for line in source_path.read_text().splitlines():
        if not line.startswith(('ATOM  ', 'HETATM')):
            continue
        atom = AtomIdentity(line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip())
        if atom in request.fixed_atoms:
            if atom in fixed or line[16].strip() or float(line[54:60]) != 1:
                raise ValueError('ambiguous fixed observation')
            fixed[atom] = np.array([float(line[i:i+8]) for i in (30, 38, 46)])
    if fixed.keys() != set(represented):
        raise ValueError('fixed coordinates do not match validated real axis')
    callback = MPSFixedAtomCallback(axis, len(features['atom_pad_mask']), fixed)
    receipt = {'source_revision': revision, 'source_sha256': request.normalized_pdb.sha256,
        'real_atoms': len(axis), 'padded_atoms': callback.padded_count, 'represented_fixed': len(fixed),
        'generated': len(request.generated_atoms), 'atomic_axis': [asdict(a) for a in axis],
        'validation': 'token traversal, feature names/elements, real unpadded prefix, atom-to-token',
        'patch_sha256': digest(Path(__file__).with_name('mps-fork-per-step-callback.patch')),
        'model_source_sha256': digest(args.engine / 'src/boltz/model/models/boltz2.py'),
        'sampler_source_sha256': digest(args.engine / 'src/boltz/model/modules/diffusionv2.py')}
    (args.output / 'axis-preflight.json').write_text(json.dumps(receipt, indent=2)+'\n')
    params = Boltz2DiffusionParams()
    params.step_scale = 1.5
    steering = BoltzSteeringParams()
    steering.fk_steering = steering.physical_guidance_update = steering.contact_guidance_update = False
    # Only the checksum-verified pinned upstream checkpoint above is trusted.
    # Its Lightning metadata contains OmegaConf objects, not tensor-only state.
    model = Boltz2.load_from_checkpoint(args.cache / 'boltz2_conf.ckpt', strict=True, map_location='cpu', weights_only=False,
        predict_args={'recycling_steps': 1, 'sampling_steps': 200, 'diffusion_samples': 1,
                      'max_parallel_samples': 1, 'write_confidence_summary': True,
                      'write_full_pae': False, 'write_full_pde': False},
        diffusion_process_args=asdict(params), ema=False, use_kernels=False,
        pairformer_args=asdict(PairformerArgsV2()), msa_args=asdict(MSAModuleArgs(use_paired_feature=True)),
        steering_args=asdict(steering))
    model.dvbfixer_step_callback = callback
    model.eval()
    def final_axis_check(_module, positional):
        batch = positional[0]
        required = ('atom_pad_mask', 'atom_to_token', 'ref_atom_name_chars', 'ref_element')
        if batch['atom_pad_mask'].shape[0] != 1:
            raise ValueError('private pilot requires one final-consumer target')
        actual = {key: batch[key][0] for key in required}
        if validate_feature_axis(axis, token_indices, actual) != numbers:
            raise ValueError('final model feature elements differ from preflight')
        (args.output / 'final-consumer-axis.json').write_text(json.dumps({
            'feature_names_elements_atom_to_token_and_padding_verified': True,
            'device': str(batch['atom_pad_mask'].device), 'real_atoms': len(axis)})+'\n')
    model.register_forward_pre_hook(final_axis_check)
    data = Boltz2InferenceDataModule(manifest=manifest, target_dir=paths[0], msa_dir=paths[1],
        mol_dir=args.cache / 'mols', num_workers=0, constraints_dir=paths[2], template_dir=paths[3], extra_mols_dir=paths[4])
    seed_everything(7, workers=True)
    trainer = Trainer(default_root_dir=args.output, accelerator='mps', devices=1, precision='32-true',
                      logger=False, enable_checkpointing=False)
    try:
        predictions = trainer.predict(model, datamodule=data, return_predictions=True)
    finally:
        (args.output / 'reinjection-steps.json').write_text(json.dumps(callback.steps, indent=2)+'\n')
        if callback.last_pre_overwrite is not None:
            np.savez_compressed(args.output / 'last-pre-overwrite.npz', coordinates=callback.last_pre_overwrite)
    if len(predictions) != 1 or predictions[0].get('exception'):
        raise RuntimeError('model prediction failed')
    xyz = predictions[0]['coords']
    np.savez_compressed(args.output / 'returned-prediction.npz', coordinates=xyz.detach().cpu().numpy())
    (args.output / 'returned-prediction.json').write_text(json.dumps({
        'device': str(xyz.device), 'shape': list(xyz.shape),
        'note': 'Lightning may move returned predictions to CPU after MPS execution'})+'\n')
    if len(callback.steps) != 200 or any(s['post_overwrite_max_angstrom'] for s in callback.steps):
        raise RuntimeError('post-update fixed-coordinate conformance failed')
    if xyz.shape != (1, callback.padded_count, 3):
        raise ValueError('unexpected returned prediction shape')
    final_consumer = json.loads((args.output / 'final-consumer-axis.json').read_text())
    if not final_consumer['device'].startswith('mps'):
        raise ValueError('model input did not execute on MPS')
    target = torch.as_tensor(callback.target, dtype=xyz.dtype, device=xyz.device)
    if not torch.equal(xyz[:, callback.indices, :], target[None]):
        raise ValueError('returned prediction changed final fixed coordinates')
    native = xyz[0, :len(axis)].detach().cpu().numpy().astype(np.float64)
    np.savez_compressed(args.output / 'native-arrays.npz', coordinates=native,
                        atom_name=np.array([a.atom_name for a in axis]))
    text = materialize_candidate_pdb(source_path.read_bytes().decode(), axis, native, names, request,
        elements=tuple(Chem.GetPeriodicTable().GetElementSymbol(n) for n in numbers))
    (args.output / 'candidate.pdb').write_text(text)
    (args.output / 'execution-evidence.json').write_text(json.dumps({
        'arm': 'private-boltz-mps-post-update-reinjection-v1', 'callbacks': 200,
        'device': str(xyz.device), 'cpu_only': 'explicit 3x3 covariance SVD',
        'fixed_after_updates_not_during_noise_augmentation': True,
        'scientific_qualification': False}, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('request', 'input', 'cache', 'engine', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    run(parser.parse_args())
