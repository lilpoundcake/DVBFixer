"""Target-bound PDB headers must not retain excluded sequence or foreign chains."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from dvbfixer.model.diffusion.contract import (
    AtomIdentity,
    ResidueIdentity,
    SequencePlacement,
    TargetSequence,
)
from dvbfixer.model.diffusion.pdb_materialize import materialize_candidate_pdb
from dvbfixer.model.diffusion.pdb_metadata import reconcile_pdb_metadata
from tests.test_diffusion_pdb_materialize import _atom_line, _request


def record(kind: str, *fields: tuple[int, str]) -> str:
    chars = list(kind.ljust(80))
    for start, text in fields:
        chars[start:start + len(text)] = text
    return ''.join(chars) + '\n'


def coordinates() -> str:
    return ''.join(_atom_line(i, i, float(i), occupancy=1, bfactor=0) for i in range(10, 15))


def helix(chain='A', first=10, last=14):
    return record('HELIX', (7, '  1'), (11, ' H1'), (15, 'ALA'), (19, chain), (21, f'{first:4}'),
        (27, 'ALA'), (31, chain), (33, f'{last:4}'), (38, ' 1'), (71, f'{last-first+1:5}'))


def sheet(identifier=' S1', chain='A', first=10, last=14, count=1, strand=1, registered=False):
    fields = [(7, f'{strand:3}'), (11, identifier), (14, f'{count:2}'), (17, 'ALA'), (21, chain), (22, f'{first:4}'),
        (28, 'ALA'), (32, chain), (33, f'{last:4}'), (38, ' 1' if registered else ' 0')]
    if registered:
        fields += [(45, 'ALA'), (49, chain), (50, f'{first:4}'), (60, 'ALA'), (64, 'A'), (65, '  11')]
    return record('SHEET', *fields)


def test_authoritative_target_replaces_stale_seqres_without_changing_coordinates():
    source = 'HEADER    KEEP ORIGINAL EXPERIMENT\nSEQRES   1 A   20  ' + ' '.join(['ALA'] * 20) + '\n' + coordinates() + 'END\n'
    normalized = reconcile_pdb_metadata(source, _request())
    assert 'SEQRES   1 A    5  ALA ALA ALA ALA ALA\n' in normalized.text
    assert 'A   20' not in normalized.text
    assert coordinates() in normalized.text
    assert source.startswith('HEADER    KEEP ORIGINAL EXPERIMENT\n')
    assert 'METADATA_SYNC' in normalized.text
    assert not reconcile_pdb_metadata(normalized.text, _request()).changes


def test_missing_optional_seqres_is_added_by_writer_but_not_required_by_auditor():
    source = coordinates() + 'END\n'
    assert 'SEQRES' in reconcile_pdb_metadata(source, _request()).text
    assert not reconcile_pdb_metadata(source, _request(), add_missing_seqres=False).changes


def test_header_endpoints_use_case_sensitive_chain_and_insertion_codes():
    source = helix() + helix('a') + helix(first=9) + coordinates() + 'END\n'
    normalized = reconcile_pdb_metadata(source, _request()).text
    assert helix() in normalized
    assert helix('a') not in normalized
    assert helix(first=9) not in normalized
    # Correct number, wrong insertion code: must not resolve to plain A10.
    inserted = helix()[:25] + 'A' + helix()[26:]
    assert inserted not in reconcile_pdb_metadata(inserted + coordinates(), _request()).text


def test_keeps_complete_sheets_but_not_partial_registration_groups():
    valid = sheet()
    partial = sheet(' S2', count=2, last=11) + sheet(' S2', 'B', first=13, count=2, strand=2, registered=True)
    normalized = reconcile_pdb_metadata(valid + partial + coordinates(), _request()).text
    assert valid in normalized
    assert ' S2' not in normalized


def test_cispep_and_sites_referencing_absent_residues_are_not_passed_through():
    good = record('CISPEP', (11, 'ALA'), (15, 'A'), (17, '  11'), (25, 'ALA'), (29, 'A'), (31, '  12'))
    bad = good[:29] + 'B' + good[30:]
    site = record('SITE', (7, '  1'), (11, ' S1'), (15, ' 2'), (18, 'ALA'), (22, 'A'), (23, '  10'), (29, 'ALA'), (33, 'B'), (34, '   1'))
    normalized = reconcile_pdb_metadata(good + bad + site + coordinates(), _request()).text
    assert good in normalized and bad not in normalized and 'SITE' not in normalized


def test_does_not_guess_database_mapping_for_an_excluded_range():
    good = record('DBREF', (7, 'TEST'), (12, 'A'), (14, '  10'), (20, '  14'))
    bad = record('DBREF', (7, 'TEST'), (12, 'A'), (14, '   1'), (20, '  20'))
    normalized = reconcile_pdb_metadata(good + bad + coordinates(), _request()).text
    assert good in normalized and bad not in normalized


def test_missing_atom_annotations_keep_only_genuinely_missing_atoms():
    missing_residue = record('REMARK', (7, '465'), (15, 'ALA'), (19, 'A'), (21, '   11'))
    excluded_residue = record('REMARK', (7, '465'), (15, 'ALA'), (19, 'A'), (21, '    1'))
    missing_atom = record('REMARK', (7, '470'), (15, 'ALA'), (19, 'A'), (20, '  10'), (25, '   CA   CB'))
    normalized = reconcile_pdb_metadata(missing_residue + excluded_residue + missing_atom + coordinates(), _request()).text
    assert 'REMARK 465' not in normalized
    lines = [line for line in normalized.splitlines() if line.startswith('REMARK 470')]
    assert len(lines) == 1 and lines[0][25:].split() == ['CB']


def test_dbref_extended_records_require_a_complete_valid_pair():
    first = record('DBREF1', (7, 'TEST'), (12, 'A'), (14, '  10'), (20, '  14'))
    second = record('DBREF2', (7, 'TEST'), (12, 'A'))
    source = first + second + coordinates()
    normalized = reconcile_pdb_metadata(source, _request()).text
    assert first in normalized and second in normalized
    assert 'DBREF' not in reconcile_pdb_metadata(first + coordinates(), _request()).text
    assert 'DBREF' not in reconcile_pdb_metadata(second + coordinates(), _request()).text


def test_link_requires_existing_atoms_and_retains_valid_record_verbatim():
    link = record('LINK', (12, ' CA '), (17, 'ALA'), (21, 'A'), (22, '  10'),
                  (42, ' CA '), (47, 'ALA'), (51, 'A'), (52, '  14'))
    absent_atom = link[:42] + ' CB ' + link[46:]
    output = reconcile_pdb_metadata(link + absent_atom + coordinates(), _request()).text
    assert link in output and absent_atom not in output


def test_assemblies_are_kept_whole_without_collapsing_chain_case():
    def assembly(identifier, chains):
        return (f'REMARK 350 BIOMOLECULE: {identifier}\n'
                f'REMARK 350 APPLY THE FOLLOWING TO CHAINS: {chains}\n'
                'REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000  0.00000\n'
                'REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000  0.00000\n'
                'REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000  0.00000\n')
    good, bad = assembly('1', 'A'), assembly('2', 'A, a')
    output = reconcile_pdb_metadata(good + bad + coordinates(), _request()).text
    assert good in output and 'BIOMOLECULE: 2' not in output
    assert not reconcile_pdb_metadata(output, _request()).changes


def test_terminal_annotations_follow_generated_terminus_without_renumbering_serials():
    source = coordinates() + 'TER      99      ALA A  10\nEND\n'
    normalized = reconcile_pdb_metadata(source, _request()).text
    terminal = next(line for line in normalized.splitlines() if line.startswith('TER'))
    assert terminal[6:11].strip() == '99' and terminal[22:26].strip() == '14'
    assert coordinates() in normalized
    assert 'TER\n' in reconcile_pdb_metadata(coordinates() + 'TER\nEND\n', _request()).text


def test_seqres_normalization_is_case_sensitive_and_preserves_consistent_records():
    upper = 'SEQRES   1 A    5  ALA ALA ALA ALA ALA\n'
    lower = 'SEQRES   1 a   20  ' + ' '.join(['ALA'] * 20) + '\n'
    atom = _atom_line(100, 1, 0., occupancy=1, bfactor=0)
    atom = atom[:21] + 'a' + atom[22:]
    request = replace(_request(), target_sequences=(*_request().target_sequences, TargetSequence('a', 'A')),
        sequence_placements=(*_request().sequence_placements, SequencePlacement('a', 1, (0,), (ResidueIdentity('a', '1'),))),
        fixed_atoms=(*_request().fixed_atoms, AtomIdentity('a', '1', '', 'CA')))
    normalized = reconcile_pdb_metadata(upper + lower + coordinates() + atom, request).text
    assert upper in normalized
    assert 'SEQRES   1 a    1  ALA\n' in normalized


def test_9cqd_fasta_and_no_terminal_output_metadata_match_the_requested_construct(tmp_path):
    from dvbfixer.diagnose.structural import check_seqres_vs_atom
    from dvbfixer.model.diffusion.scope import _CANONICAL_ONE_TO_THREE
    from dvbfixer.model.diffusion_cli import _target_sequences, build_cli_diffusion_request

    source = Path(__file__).parent / 'fixtures/protein_context/9cqd_split.pdb'
    fasta = source.with_suffix('.fasta')
    original = source.read_bytes()
    request = build_cli_diffusion_request(source, _target_sequences(source, str(fasta)), seeds=(7,),
        profile='protpardelle-1c-mps-local-v1', no_terminal=True)
    assert {t.chain: len(t.sequence) for t in request.target_sequences} == {'M': 27, 'W': 214, 'X': 221}
    assert len(request.generated_atoms) == 30
    names: dict[tuple[str, str, str], str] = {}
    for gap in request.gaps:
        sequence = next(target.sequence for target in request.target_sequences if target.chain == gap.chain)
        for offset, residue in enumerate(gap.generated_residues, gap.target_interval.start):
            names[(residue.chain, residue.residue_number, residue.insertion_code)] = _CANONICAL_ONE_TO_THREE[sequence[offset]]
    rendered = materialize_candidate_pdb(original.decode(), request.generated_atoms,
        [(float(i), 0., 0.) for i in range(len(request.generated_atoms))],
        tuple(names[(atom.chain, atom.residue_number, atom.insertion_code)] for atom in request.generated_atoms), request)
    output = tmp_path / 'output.pdb'
    output.write_text(rendered)
    assert not check_seqres_vs_atom(output)
    assert [line for line in rendered.splitlines() if line.startswith('SSBOND')] == [line for line in original.decode().splitlines() if line.startswith('SSBOND')]
    assert all(line in rendered for line in original.decode().splitlines(keepends=True) if line.startswith(('ATOM', 'ANISOU')))
    assert source.read_bytes() == original
    assert not reconcile_pdb_metadata(rendered, request, add_missing_seqres=False).changes
