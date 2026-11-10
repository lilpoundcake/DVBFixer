"""Logging-only native component-checkpoint amendment, all original starts.

Enable native verbose derivative diagnostics, check ALL exposed numerical option
getters against the quiet checker, capture fd output with native flush controls.
No Python override of C++ virtual functions or numerical profile change.
"""
from __future__ import annotations

import argparse
import ctypes
import json
from pathlib import Path
from unittest.mock import patch

from scripts import research_chi_checkpoints as parent
from scripts.research_chi_termination import capture_native
from scripts.research_chi_termination_flush import flush_checked
from scripts.research_dunbrack_prior import sha, write

NUMERICAL_GETTERS=('min_type','minimize_tolerance','use_nblist','nblist_auto_update','max_iter',
    'gmax_cutoff_for_convergence','armijo_max_step_limit','linmin_deriv_cutoff',
    'ax_init','bx_init','xx_init','brent_abs_tolerance','cmaes_lambda','cmaes_rgsigma','ga_mutation_probability','silent')
SETTABLE=('nblist_auto_update','max_iter','gmax_cutoff_for_convergence','armijo_max_step_limit',
          'linmin_deriv_cutoff','cmaes_lambda','cmaes_rgsigma','ga_mutation_probability','silent')


def check_equal(before,after):
    if set(before)!=set(NUMERICAL_GETTERS) or before!=after:
        raise ValueError('native verbose diagnostics changed numerical options')


def plan(parent_path,output):
    raw=parent_path.read_bytes()
    record=json.loads(raw)
    for p,d in {**record['inputs'],**record['source_hashes']}.items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('bound quiet-checkpoint state changed')
    for p in (parent_path,Path('tmp/pdbfixer-research/chi-checkpoints-run-01/checkpoint-evidence.json')):
        record['inputs'][str(p)]=sha(p.read_bytes())
    for filename in ('scripts/research_chi_checkpoint_verbose.py','scripts/research_chi_termination.py',
                     'scripts/research_chi_termination_flush.py'):
        record['source_hashes'][filename]=sha(Path(filename).read_bytes())
    record['init_flags']+=' -unmute core.optimization'
    record['verbose_checkpoint_amendment']=dict(native_verbose=True,native_stdout=True,
        outer_checkpoint_trace_flags='quiet outer adapter flags are recorded by original runner; actual native override flags/options are bound in verbose-adapter.json',
        numerical_getters=list(NUMERICAL_GETTERS),
        enrollment='ALL42 original starts; no narrowing to the two unresolved endpoints',
        interpretation='native derivative component checkpoints, not full function/iteration/return-branch trace; all original exact replay comparisons remain required',
        controls='checked native stdio flushes and synthetic real-tracer BEGIN/END controls per native call')
    write(output,record)


def run(plan_path,output):
    record=json.loads(plan_path.read_bytes())
    if 'verbose_checkpoint_amendment' not in record:
        raise ValueError('separate native verbose amendment plan required')
    import pyrosetta.rosetta.core.optimization as optimization
    from pyrosetta.rosetta.basic import Tracer

    original=optimization.AtomTreeMinimizer
    option_type=optimization.MinimizerOptions
    inventory=[(case,cell) for case in record['cases'] for cell in case['cells']]
    created=[]
    captures=[]
    libc=ctypes.CDLL(None)
    flush=libc.fflush
    flush.argtypes=[ctypes.c_void_p]
    flush.restype=ctypes.c_int

    class VerboseNative:
        def __init__(self):
            self.native=original()
            self.case,self.cell=inventory[len(created)]
            created.append(self)

        def run(self,pose,movemap,sf,options):
            before={name:getattr(options,name)() for name in NUMERICAL_GETTERS}
            verbose=option_type(options.min_type(),options.minimize_tolerance(),options.use_nblist(),True,True)
            for name in SETTABLE:
                getattr(verbose,name)(getattr(options,name)())
            verbose.deriv_check_to_stdout(True)
            after={name:getattr(verbose,name)() for name in NUMERICAL_GETTERS}
            check_equal(before,after)
            if not verbose.deriv_check() or not verbose.deriv_check_verbose() or not verbose.deriv_check_to_stdout():
                raise ValueError('native verbose diagnostic flags differ')
            path=output/self.case['id']/self.cell['declaration']['id']/'native-components.log'
            tracer=Tracer('core.optimization')
            flush_checked(flush)
            with capture_native(path):
                text='DVBFIXER_SYNTHETIC_COMPONENT_CAPTURE_BEGIN\n'
                tracer(200).write(text,len(text))
                tracer.flush_all_channels()
                try:
                    result=self.native.run(pose,movemap,sf,verbose)
                finally:
                    text='DVBFIXER_SYNTHETIC_COMPONENT_CAPTURE_END\n'
                    tracer(200).write(text,len(text))
                    tracer.flush_all_channels()
                    flush_checked(flush)
            text=path.read_text()
            if text.count('DVBFIXER_SYNTHETIC_COMPONENT_CAPTURE_BEGIN')!=1 or text.count('DVBFIXER_SYNTHETIC_COMPONENT_CAPTURE_END')!=1 or 'deriv_dev:' not in text:
                raise ValueError('native component capture controls/data missing')
            captures.append(dict(id=self.case['id'],index=self.cell['declaration']['index'],log_path=str(path),
                log_sha256=sha(path.read_bytes()),numerical_options_before=before,numerical_options_after=after,
                actual_native_flags=dict(deriv_check=True,deriv_check_verbose=True,deriv_check_to_stdout=True)))
            return result

        def deriv_check_result(self):
            return self.native.deriv_check_result()

    with patch.object(optimization,'AtomTreeMinimizer',VerboseNative):
        parent.run(plan_path,output)
    if optimization.AtomTreeMinimizer is not original or len(captures)!=len(inventory):
        raise ValueError('verbose native adapter not restored/inventory incomplete')
    write(output/'verbose-adapter.json',dict(captures=captures,hook_restored=True,numerical_settings_changed=False,
        native_coordinates_read=False,python_virtual_callback=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('plan')
    p.add_argument('--parent',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p=sub.add_parser('run')
    p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.command=='plan':
        plan(args.parent,args.output)
    else:
        run(args.plan,args.output)
