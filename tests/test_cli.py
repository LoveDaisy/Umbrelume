import json

import numpy as np
import pytest

from eclipse.cli import main


def test_cli_writes_physical_spectra_and_reproducible_settings(tmp_path):
    main(["simulate", "--out", str(tmp_path), "--no-plots", "--radial-step-km", "500",
          "--wavelength-step-nm", "20", "--height-step-m", "500"])
    result = np.load(tmp_path/"spectra.npz")
    metadata = json.loads((tmp_path/"metadata.json").read_text())
    assert result["irradiance_W_m2_nm"].shape == (19, 21)
    assert result["radius_m"][-1] == 9e6
    assert np.all(result["irradiance_W_m2_nm"] >= 0)
    assert np.all(result["irradiance_W_m2_nm"] <= result["irradiance_no_ozone_W_m2_nm"])
    assert metadata["atmosphere"]["ozone_du"] == 300
    assert metadata["aerosol_preset"] == "opac"
    layers = metadata["aerosol"]["layers"]
    assert [layer["aod550"] for layer in layers] == [.013, .005]
    assert layers[0]["source"]["candidate"] == "free_troposphere_reference"
    assert layers[1]["source"]["candidate"] == "stratosphere_background"
    assert (layers[1]["source"]["opac1998_sha256"]
            == "68a20ba9488db1e0f65a9976d4ee46e6c879f89054fcff583ba67bb5d68d425c")
    assert "stratospheric-background" in metadata["approximations"][3]
    assert metadata["numerics"]["height_step_m"] == 500
    assert (tmp_path/"rays.npz").is_file()
    assert (tmp_path/"profile.csv").is_file()


def test_cli_default_matches_explicit_opac_and_npz_layer_shape(tmp_path):
    default_out, explicit_out = tmp_path/"default", tmp_path/"explicit"
    common = ["--no-plots", "--radial-step-km", "1000", "--wavelength-step-nm", "50",
              "--height-step-m", "500"]
    main(["simulate", "--out", str(default_out), *common])
    main(["simulate", "--out", str(explicit_out), "--aerosol", "opac", *common])
    default = np.load(default_out/"rays.npz")
    explicit = np.load(explicit_out/"rays.npz")
    assert set(default.files) == set(explicit.files)
    assert "aerosol_layer_taus" in default.files
    assert default["aerosol_layer_taus"].shape == (2, default["tangent_altitude_m"].size)
    for key in default.files:
        np.testing.assert_array_equal(default[key], explicit[key])


def test_cli_opac_ft_keeps_single_layer_semantics(tmp_path):
    main(["simulate", "--out", str(tmp_path), "--aerosol", "opac-ft", "--no-plots",
          "--radial-step-km", "1000", "--wavelength-step-nm", "50", "--height-step-m", "500"])
    metadata = json.loads((tmp_path/"metadata.json").read_text())
    rays = np.load(tmp_path/"rays.npz")
    assert metadata["aerosol_preset"] == "opac-ft"
    assert len(metadata["aerosol"]["layers"]) == 1
    assert metadata["aerosol"]["layers"][0]["aod550"] == .013
    assert rays["aerosol_layer_taus"].shape == (1, rays["tangent_altitude_m"].size)


def test_cli_aod550_override_touches_only_free_troposphere_layer(tmp_path):
    main(["simulate", "--out", str(tmp_path), "--aerosol", "opac", "--aod550", "0.5",
          "--no-plots", "--radial-step-km", "1000", "--wavelength-step-nm", "50",
          "--height-step-m", "500"])
    metadata = json.loads((tmp_path/"metadata.json").read_text())
    assert [layer["aod550"] for layer in metadata["aerosol"]["layers"]] == [.5, .005]
    assert metadata["aerosol_aod_source"] == "User override; not an independent prior"


def test_cli_aod550_zero_is_the_molecular_counterfactual(tmp_path):
    zero_out, none_out, ft_zero_out = tmp_path/"zero", tmp_path/"none", tmp_path/"ft-zero"
    common = ["--no-plots", "--radial-step-km", "1000", "--wavelength-step-nm", "50",
              "--height-step-m", "500"]
    main(["simulate", "--out", str(zero_out), "--aerosol", "opac", "--aod550", "0", *common])
    main(["simulate", "--out", str(none_out), "--aerosol", "none", *common])
    main(["simulate", "--out", str(ft_zero_out), "--aerosol", "opac-ft",
          "--aod550", "0", *common])
    zero = np.load(zero_out/"rays.npz")
    none = np.load(none_out/"rays.npz")
    ft_zero = np.load(ft_zero_out/"rays.npz")
    for key in zero.files:
        np.testing.assert_array_equal(zero[key], none[key])
        np.testing.assert_array_equal(zero[key], ft_zero[key])
    zero_metadata = json.loads((zero_out/"metadata.json").read_text())
    none_metadata = json.loads((none_out/"metadata.json").read_text())
    ft_zero_metadata = json.loads((ft_zero_out/"metadata.json").read_text())
    assert zero_metadata["aerosol"] is None and none_metadata["aerosol"] is None
    # The preset field keeps the requested choice while the medium is
    # molecular-only; the disabled source field is the disambiguator, so the
    # combination's semantics stay pinned against silent drift — for both
    # aerosol presets that accept --aod550.
    assert zero_metadata["aerosol_preset"] == "opac"
    assert zero_metadata["aerosol_aod_source"] == "Disabled"
    assert "No aerosol extinction" in zero_metadata["approximations"]
    assert ft_zero_metadata["aerosol"] is None
    assert ft_zero_metadata["aerosol_preset"] == "opac-ft"
    assert ft_zero_metadata["aerosol_aod_source"] == "Disabled"


def test_cli_molecular_counterfactual_has_consistent_metadata(tmp_path):
    main(["simulate","--out",str(tmp_path),"--aerosol","none","--no-plots",
          "--radial-step-km","1000","--wavelength-step-nm","50","--height-step-m","500"])
    metadata=json.loads((tmp_path/"metadata.json").read_text())
    assert metadata["aerosol"] is None
    assert metadata["aerosol_aod_source"]=="Disabled"
    assert "No aerosol extinction" in metadata["approximations"]
    np.testing.assert_array_equal(np.load(tmp_path/"rays.npz")["aerosol_tau550"],0)


def test_cli_calibrated_preset_metadata_layers_and_marker(tmp_path):
    import eclipse.cli as cli_module
    main(["simulate", "--out", str(tmp_path), "--aerosol", "calibrated:june2011",
          "--no-plots", "--radial-step-km", "1000", "--wavelength-step-nm", "50",
          "--height-step-m", "500"])
    metadata = json.loads((tmp_path/"metadata.json").read_text())
    assert metadata["aerosol_preset"] == "calibrated:june2011"
    layers = metadata["aerosol"]["layers"]
    assert [layer["aod550"] for layer in layers] == [.013, .005, 0.025865877276886398]
    assert layers[2]["bottom_m"] == 0.0 and layers[2]["top_m"] == 10000.0
    source = layers[2]["source"]
    assert source["source_class"] == "calibrated_registry"
    assert source["preset"] == "june2011"
    assert source["event"] == "june_2011"
    assert source["kind"] == "calibrated"
    assert source["calibrated_not_physical_prior"] is True
    assert "Calibrated" in metadata["aerosol_aod_source"]
    assert "not an independent prior" in metadata["aerosol_aod_source"]
    marker = metadata["approximations"][3]
    assert "not a physical prior" in marker and "never a default" in marker
    rays = np.load(tmp_path/"rays.npz")
    assert rays["aerosol_layer_taus"].shape == (3, rays["tangent_altitude_m"].size)
    # The registry must be the only calibration input and the event mapping
    # must come from its 'event' field (review round 2 item 4).
    registry = json.loads(cli_module.CALIBRATED_PRESETS_JSON.read_text())
    assert registry["presets"]["june2011"]["event"] == "june_2011"


def test_cli_calibrated_negative_paths_exit_loud(tmp_path):
    for argv in (
        ["--aerosol", "calibrated:june2011", "--aod550", "0.1"],
        ["--aerosol", "calibrated:does_not_exist"],
        ["--aerosol", "calibrated:"],
        ["--aerosol", "bogus-preset"],
    ):
        with pytest.raises(SystemExit):
            main(["simulate", "--out", str(tmp_path), "--no-plots", *argv,
                  "--radial-step-km", "1000", "--wavelength-step-nm", "50",
                  "--height-step-m", "500"])


def test_cli_default_and_opac_ft_paths_never_touch_registry(tmp_path, monkeypatch):
    import eclipse.cli as cli_module

    def boom(*args, **kwargs):
        raise AssertionError("the default paths must not load the preset registry")
    monkeypatch.setattr(cli_module, "load_calibrated_preset", boom)
    common = ["--no-plots", "--radial-step-km", "1000", "--wavelength-step-nm", "50",
              "--height-step-m", "500"]
    main(["simulate", "--out", str(tmp_path/"default"), *common])
    main(["simulate", "--out", str(tmp_path/"ft"), "--aerosol", "opac-ft", *common])
    main(["simulate", "--out", str(tmp_path/"none"), "--aerosol", "none", *common])
    metadata = json.loads((tmp_path/"default"/"metadata.json").read_text())
    assert metadata["aerosol_preset"] == "opac"
    assert [layer["aod550"] for layer in metadata["aerosol"]["layers"]] == [.013, .005]


def test_cli_calibrated_provenance_negative_paths_fail_loud(tmp_path, monkeypatch):
    """Major review round 1 fix: the loader's provenance gate is field-level
    and its rejection paths are test-covered (marker missing/false, sha256
    missing/empty, malformed event, structurally broken entry)."""
    import eclipse.cli as cli_module
    registry = json.loads(cli_module.CALIBRATED_PRESETS_JSON.read_text())
    entry = json.loads(json.dumps(registry["presets"]["june2011"]))
    common = ["--no-plots", "--radial-step-km", "1000", "--wavelength-step-nm",
              "50", "--height-step-m", "500"]

    def use(preset_entry):
        path = tmp_path/"registry.json"
        path.write_text(json.dumps({"presets": {"june2011": preset_entry}}))
        monkeypatch.setattr(cli_module, "CALIBRATED_PRESETS_JSON", path)

    def broken(mutate):
        variant = json.loads(json.dumps(entry))
        mutate(variant)
        use(variant)
        with pytest.raises(SystemExit):
            main(["simulate", "--out", str(tmp_path/"out"), "--aerosol",
                  "calibrated:june2011", *common])

    broken(lambda e: e["provenance"].pop("calibrated_not_physical_prior"))
    broken(lambda e: e["provenance"].__setitem__("calibrated_not_physical_prior", False))
    broken(lambda e: e["provenance"].pop("sha256"))
    broken(lambda e: e["provenance"].__setitem__("sha256", ""))
    broken(lambda e: e["provenance"].pop("inversion_product"))
    broken(lambda e: e.pop("event"))
    broken(lambda e: e.__setitem__("event", ""))
    broken(lambda e: e.pop("increment_layer"))
    broken(lambda e: e.__setitem__("increment_layer", {"aod550": 0.1, "unknown_key": 1}))
    # A structurally broken registry file also exits through parser.error.
    path = tmp_path/"registry.json"
    path.write_text('{"presets": null}')
    monkeypatch.setattr(cli_module, "CALIBRATED_PRESETS_JSON", path)
    with pytest.raises(SystemExit):
        main(["simulate", "--out", str(tmp_path/"out"), "--aerosol",
              "calibrated:june2011", *common])
