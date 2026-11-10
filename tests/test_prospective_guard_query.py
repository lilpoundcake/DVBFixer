"""Dependency-isolated query entry must not eagerly import scientific adapters."""
import ast
from pathlib import Path


def test_query_entry_has_no_eager_standard_science_import():
    tree=ast.parse(Path('scripts/research_prospective_guard_query.py').read_text())
    modules={node.module for node in tree.body if isinstance(node,ast.ImportFrom)}
    assert 'scripts.prospective_chi_guard_cohort' not in modules
    assert 'scripts.research_prospective_guard_views' not in modules
    assert 'scripts.research_dunbrack_prior' in modules
