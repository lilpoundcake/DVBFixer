"""Explicit native-stdio flush amendment; preserve failed buffered-log attempt.

Keep logging/numerical replay unchanged; flush C stdio before fd boundaries and
check synthetic tracer markers inside every capture. No synthetic optimizer
warnings, changed numerical settings or attribution of old unscoped messages.
"""
from __future__ import annotations

import argparse
import ctypes
import json
from contextlib import contextmanager
from pathlib import Path
from unittest.mock import patch

from scripts import research_chi_termination as parent
from scripts.research_dunbrack_prior import sha, write


def flush_checked(flush):
    if flush(None)!=0:
        raise OSError('native stdio fflush failed')


def plan(parent_path,output):
    raw=parent_path.read_bytes()
    record=json.loads(raw)
    for p,d in {**record['inputs'],**record['source_hashes']}.items():
        if sha(Path(p).read_bytes())!=d:
            raise ValueError('original logging state changed')
    for p in (parent_path,Path('tmp/pdbfixer-research/chi-termination-run-01.log'),
              Path('tmp/pdbfixer-research/chi-termination-run-01/termination-evidence.json'),
              Path('tmp/pdbfixer-research/chi-termination-flush-diagnostic-01.log'),
              Path('tmp/pdbfixer-research/chi-termination-flush-diagnostic-02.log'),
              Path('tmp/pdbfixer-research/chi-termination-flush-native-02.log')):
        record['inputs'][str(p)]=sha(p.read_bytes())
    record['source_hashes']['scripts/research_chi_termination_flush.py']=sha(Path(__file__).read_bytes())
    record['native_flush_amendment']='explicit pre/post C-stdio fflush around each fd capture plus real tracer begin/end control markers; failed unscoped original logs retained; no numerical setting changes'
    write(output,record)


def run(plan_path,output):
    record=json.loads(plan_path.read_bytes())
    if 'native_flush_amendment' not in record:
        raise ValueError('separate native-flush amendment plan required')
    libc=ctypes.CDLL(None)
    flush=libc.fflush
    flush.argtypes=[ctypes.c_void_p]
    flush.restype=ctypes.c_int
    original=parent.capture_native
    captures=[]

    @contextmanager
    def flushed(path):
        from pyrosetta.rosetta.basic import Tracer

        tracer=Tracer('core.optimization.Minimizer')
        flush_checked(flush)
        with original(path):
            text='DVBFIXER_SYNTHETIC_CAPTURE_CONTROL_BEGIN\n'
            tracer(200).write(text,len(text))
            tracer.flush_all_channels()
            try:
                yield
            finally:
                text='DVBFIXER_SYNTHETIC_CAPTURE_CONTROL_END\n'
                tracer(200).write(text,len(text))
                tracer.flush_all_channels()
                flush_checked(flush)
        text=path.read_text()
        if text.count('DVBFIXER_SYNTHETIC_CAPTURE_CONTROL_BEGIN')!=1 or text.count('DVBFIXER_SYNTHETIC_CAPTURE_CONTROL_END')!=1:
            raise ValueError('native tracer capture controls missing/duplicated')
        captures.append(dict(path=str(path),sha256=sha(path.read_bytes()),synthetic_markers_verified=True))

    with patch.object(parent,'capture_native',flushed):
        parent.run(plan_path,output)
    if parent.capture_native is not original:
        raise ValueError('native flush adapter not restored')
    if len(captures)!=sum(len(c['cells']) for c in record['cases']):
        raise ValueError('native capture inventory incomplete')
    write(output/'native-flush-adapter.json',dict(hook_restored=True,captures=captures,
        numerical_settings_changed=False,synthetic_optimizer_warnings=False))


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
