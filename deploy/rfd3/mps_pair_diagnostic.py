"""Instrument the two native pair additions after reference accumulation."""
import ast
import json
import os
from pathlib import Path

import torch


def measure(label, operation):
    path = Path(os.environ['DVBFIXER_RFD3_TRANSPORT_OUTPUT']).parent/'pair-operations.jsonl'
    def record(stage, result=None):
        torch.mps.synchronize()
        row = {'operation':label, 'stage':stage, 'driver':torch.mps.driver_allocated_memory(),
               'live':torch.mps.current_allocated_memory()}
        if result is not None:
            row['shape'] = list(result.shape)
        with path.open('a') as stream:
            stream.write(json.dumps(row)+'\n')
        if row['driver'] > 6*1024**3:
            raise RuntimeError('diagnostic boundary exceeded frozen 6 GiB driver threshold')
    record('enter')
    value = operation()
    record('exit', value)
    return value


def instrument(code):
    tree = ast.parse(code)
    found = []
    def wrap(label, node):
        call = ast.parse(f'_measure({label!r}, lambda: None)',mode='eval').body
        call.args[1].body = node
        return call
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign) or ast.unparse(node.targets[0]) != 'P_LL':
            continue
        expr = node.value
        if not isinstance(expr, ast.BinOp) or ast.unparse(expr.left) != 'P_LL':
            continue
        rhs = ast.unparse(expr.right)
        if rhs.startswith('self.process_single_l('):
            label = 'single-atom-broadcast'
        elif rhs.startswith('self.process_z('):
            label = 'token-to-atom-gather'
        else:
            continue
        found.append(label)
        expr.right = wrap(label, expr.right)
        node.value = wrap(label+'-add-to-pair',expr)
    if found != ['single-atom-broadcast','token-to-atom-gather']:
        raise ValueError(f'unexpected native operations: {found}')
    return ast.unparse(ast.fix_missing_locations(tree))+'\n'
