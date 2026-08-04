from dvbfixer.doctor import collect_capabilities, main


def test_capability_report_has_stable_sections():
    report = collect_capabilities()
    assert set(report) == {"python_packages", "executables", "openmm_platforms"}
    assert "OpenMM" in report["python_packages"]
    assert "tleap" in report["executables"]


def test_json_output_is_machine_readable(capsys):
    main(["--format", "json"])
    assert '"python_packages"' in capsys.readouterr().out
