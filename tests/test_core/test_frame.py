import pytest
import molpy as mp
import numpy.testing as npt
import xarray as xr

class TestFrame:

    @pytest.fixture()
    def frame(self):
        atoms = xr.Dataset(
            {
                'id': ('index', [1, 2, 3, 4]),
                'type': ('index', [1, 2, 3, 4]),
                'x': ('index', [0, 1, 2, 3]),
                'y': ('index', [0, 1, 2, 3]),
                'z': ('index', [0, 1, 2, 3]),
            }
        )
        return mp.Frame({'atoms': atoms}, style="atomic")
    
    def test_slice(self, frame):
        assert isinstance(frame['atoms'], xr.Dataset)
        assert set(frame['atoms'].data_vars.keys()) == {'id','type','x','y','z'}
        assert frame['atoms']['id'].to_numpy().tolist() == [1, 2, 3, 4]
        
    def test_init_with_dataframe(self):
        data = {
            'atoms': xr.Dataset(
                {
                    'id': ('index', [1, 2]),
                    'type': ('index', [1, 2]),
                    'x': ('index', [0.0, 1.0]),
                    'y': ('index', [0.0, 1.0]),
                    'z': ('index', [0.0, 1.0]),
                }
            )
        }
        frame = mp.Frame(data)
        assert 'atoms' in frame
        assert isinstance(frame['atoms'], xr.Dataset)

    def test_init_with_dict(self):
        data = {
            'atoms': xr.Dataset(
                {
                    'id': ('index', [1, 2]),
                    'type': ('index', [1, 2]),
                    'x': ('index', [0.0, 1.0]),
                    'y': ('index', [0.0, 1.0]),
                    'z': ('index', [0.0, 1.0]),
                }
            )
        }
        frame = mp.Frame(data)
        assert 'atoms' in frame
        assert isinstance(frame['atoms'], xr.Dataset)

    def test_concat(self, frame):
        frame2 = mp.Frame({
            'atoms': xr.Dataset(
                {
                    'id': ('index', [5, 6]),
                    'type': ('index', [5, 6]),
                    'x': ('index', [4, 5]),
                    'y': ('index', [4, 5]),
                    'z': ('index', [4, 5]),
                }
            )
        })
        concatenated = mp.Frame.from_frames([frame, frame2])
        npt.assert_equal(concatenated['atoms']['id'],  [1, 2, 3, 4, 5, 6])

    def test_split(self, frame):
        split_frames = frame.split([1, 1, 2, 2])
        assert len(split_frames) == 2
        npt.assert_equal(split_frames[0]['atoms']['id'],  [1, 2])
        npt.assert_equal(split_frames[1]['atoms']['id'], [3, 4])

    def test_box_property(self):
        box = mp.Box()
        frame = mp.Frame()
        frame.box = box
        assert frame.box == box
        frame.box = None
        assert frame.box is None

    def test_to_struct(self, frame):
        frame['bonds'] = xr.Dataset({'i': ('index', [1]), 'j': ('index', [2])})
        struct = frame.to_struct()
        assert 'atoms' in struct
        assert 'bonds' in struct
        assert len(struct['atoms']) == 4
        assert len(struct['bonds']) == 1

    def test_copy(self, frame):
        frame_copy = frame.copy()
        assert frame_copy is not frame
        assert frame_copy['atoms'] == frame['atoms']

    def test_add_operator(self, frame):
        frame2 = mp.Frame({
            'atoms': xr.Dataset(
                {
                    'id': ('index', [5, 6]),
                    'type': ('index', [5, 6]),
                    'x': ('index', [4, 5]),
                    'y': ('index', [4, 5]),
                    'z': ('index', [4, 5]),
                }
            )
        })
        combined = frame + frame2
        assert combined['atoms']['id'].tolist() == [1, 2, 3, 4, 5, 6]

    def test_mul_operator(self, frame):
        multiplied = frame * 2
        assert len(multiplied) == 5

    def test_init_all_atom_frame(self):
        frame = mp.Frame(style="atomic")
        assert isinstance(frame, mp.AllAtomFrame)
        frame["a"] = 0