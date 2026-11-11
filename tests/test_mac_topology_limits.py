import numpy as np
from scripts.run_mac_topology_limits_controls import (
    disulfide_control,
    linked_control,
    trajectory_control,
)


def test_legitimate_linked_geometry_is_not_a_fatal_topological_defect():
    result = linked_control()
    assert not result["segments"]["alarms"]
    assert len(result["surface"]["piercings"]) == 1
    assert not result["surface"]["publication_gate"]
    assert not result["surface"]["piercings"][0]["topological_defect_established"]
    moved = linked_control(np.eye(3), np.array([10., -20., 3.]))
    assert len(moved["surface"]["piercings"]) == 1


def test_clean_endpoints_do_not_certify_unsaved_trajectory():
    result = trajectory_control()
    assert not result["snapshots"]["start"]["piercings"]
    assert not result["snapshots"]["end"]["piercings"]
    assert result["snapshots"]["intermediate"]["piercings"]
    assert result["endpoint_cleanliness_proves_no_passage"] is False


def test_shared_explicit_disulfide_parser_preserves_case_and_icode():
    result = disulfide_control()
    assert result["SSBOND_preserved"] and result["chain_case_and_insertion_code_preserved"]
    assert not result["disulfide_surface_or_knot_qualification"]
