import ast
from pathlib import Path
from types import SimpleNamespace

from scripts.run_1non_diffusion_full import SEED, native_arguments
from scripts.run_mac_reconstruction_guarded import CHILD


def test_native_input_is_complete_joint_tetramer_not_local_or_modeller():
    plan = {"parent_request": "parent.json", "source": "protein-only.pdb", "checkpoint": "cc94.pth",
            "refinement_python": "cpu-python"}
    args = native_arguments(plan, Path("native"), prepare_only=True)
    assert "--complete-observed-context" in args and "--whole-structure" in args
    assert args[args.index("--gap-indices")+1:args.index("--whole-structure")] == ["0", "1", "2", "3"]
    assert args[args.index("--source")+1] == args[args.index("--baseline")+1] == "protein-only.pdb"
    assert args[args.index("--seed")+1] == str(SEED)
    assert "--legacy-local-union" not in args and "--prepare-only" in args


def test_existing_device_guard_accounts_for_protpardelle_tensors():
    # Execute the actual embedded hook functions without importing Torch or
    # starting its monitoring thread/model; test the newly used engine prefix.
    nodes = [n for n in ast.parse(CHILD).body if isinstance(n, ast.FunctionDef) and n.name in {"tensors", "hook"}]
    class Tensor:
        device = SimpleNamespace(type="mps")
    report = {"model_mps_calls": 0, "module_classes": {}}
    namespace = {"torch": SimpleNamespace(Tensor=Tensor), "report": report}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), "guard-hook", "exec"), namespace)
    module = type("Block", (), {"__module__": "protpardelle.core.models"})()
    namespace["hook"](module, (), {"x": Tensor()})
    assert report["model_mps_calls"] == 1
    assert report["module_classes"] == {"protpardelle.core.models.Block": 1}


def test_inpainting_is_explicit_private_flag_not_a_change_to_default_arguments():
    plan = {"parent_request": "parent.json", "source": "protein-only.pdb", "checkpoint": "cc94.pth",
            "refinement_python": "cpu-python"}
    assert "--fixed-scaffold-inpainting" not in native_arguments(plan, Path("native"))
    plan["fixed_scaffold_inpainting"] = True
    assert "--fixed-scaffold-inpainting" in native_arguments(plan, Path("native"))
