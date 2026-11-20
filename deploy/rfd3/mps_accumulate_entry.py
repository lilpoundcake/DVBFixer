"""Private 1NON technical replay: consume reference tiles in the owned motif buffer."""
import importlib.util
import inspect
import json
import os
import runpy
import textwrap
from pathlib import Path

import torch
from rfd3.model.layers.encoders import TokenInitializer

ROOT = Path(os.environ['DVBFIXER_RFD3_TRANSPORT_OUTPUT']).parent
SPEC = importlib.util.spec_from_file_location('private_tile_accumulation', Path(__file__).with_name('mps_pairwise_extended.py'))
HELPER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(HELPER)


def accumulate(module, positions, mask, owned):
    if positions.device.type != 'mps' or positions.shape != (5592, 3):
        raise ValueError('expected complete 1NON native atom axis')
    result = HELPER.reference_rows(module, positions, mask, rows=64, accumulate_into=owned)
    torch.mps.synchronize()
    torch.mps.empty_cache()
    (ROOT/'combined.json').write_text(json.dumps({'shape':list(result.shape), 'same_buffer':result is owned,
        'driver_bytes_after_trim':torch.mps.driver_allocated_memory()})+'\n')
    return result


if __name__=='__main__':
    native = TokenInitializer.forward
    source = textwrap.dedent(inspect.getsource(native))
    old = 'P_LL = P_LL + self.ref_pos_embedder(f["ref_pos"], valid_mask)'
    if source.count(old) != 1:
        raise ValueError('unexpected initializer implementation')
    namespace = dict(native.__globals__, _combine=accumulate)
    code = source.replace(old, 'P_LL = _combine(self.ref_pos_embedder, f["ref_pos"], valid_mask, P_LL)')
    if os.environ.get('DVBFIXER_PRIVATE_PAIR_FINISH') == '1':
        from mps_pair_finish import finish, install
        namespace['_finish'] = finish
        code = install(code)
        if os.environ.get('DVBFIXER_PRIVATE_GATHER_BIAS') == '1':
            from mps_pair_finish import install_local_bias
            install_local_bias()
    if os.environ.get('DVBFIXER_PRIVATE_PAIR_DIAGNOSTIC') == '1':
        from mps_pair_diagnostic import instrument, measure
        namespace['_measure'] = measure
        code = instrument(code)
    (ROOT/'executed-initializer.py').write_text(code)
    exec(compile(code,'<private-pair-accumulator>','exec'),namespace)
    TokenInitializer.forward = namespace['forward']
    runpy.run_path(str(ROOT/'memory_entry.py'),run_name='__main__')
