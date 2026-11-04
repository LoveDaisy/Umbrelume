import json

import numpy as np

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
    assert metadata["aerosol"]["aod550"] == .013
    assert metadata["aerosol_preset"] == "opac-ft"
    assert metadata["numerics"]["height_step_m"] == 500
    assert (tmp_path/"rays.npz").is_file()
    assert (tmp_path/"profile.csv").is_file()


def test_cli_molecular_counterfactual_has_consistent_metadata(tmp_path):
    main(["simulate","--out",str(tmp_path),"--aerosol","none","--no-plots",
          "--radial-step-km","1000","--wavelength-step-nm","50","--height-step-m","500"])
    metadata=json.loads((tmp_path/"metadata.json").read_text())
    assert metadata["aerosol"] is None
    assert metadata["aerosol_aod_source"]=="Disabled"
    assert "No aerosol extinction" in metadata["approximations"]
    np.testing.assert_array_equal(np.load(tmp_path/"rays.npz")["aerosol_tau550"],0)
