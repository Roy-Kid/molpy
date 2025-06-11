import molpy as mp

class TestAmberACReader:

    def test_amber_ac(self, test_data_path):
        fpath = test_data_path / "data/ac/H.ac"
        frame = mp.io.read_amber_ac(fpath)
        assert isinstance(frame, mp.Frame)
        assert "atoms" in frame
        assert len(frame["atoms"].data_vars) == 8  # check n_fields
        assert frame["atoms"].sizes["index"] == 18
        assert frame["bonds"].sizes["index"] == 18
