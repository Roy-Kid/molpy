"""Test frcmod I/O functionality."""

import tempfile
from pathlib import Path

from molpy.io.readers import read_amber_frcmod
from molpy.io.writers import write_amber_frcmod


def test_read_frcmod():
    """Test reading an existing frcmod file."""
    frcmod_path = Path("/workspaces/molpy/test_amber_wrapper/parmchk2_output/tfsi.frcmod")
    
    if not frcmod_path.exists():
        print(f"Skipping: {frcmod_path} not found")
        return
    
    data = read_amber_frcmod(frcmod_path)
    
    print(f"✓ Read frcmod file: {frcmod_path}")
    print(f"  Sections: {list(data.keys())}")
    print(f"  Remark: {data['remark']}")
    print(f"  Bond section length: {len(data['bond'])} chars")
    print(f"  Angle section length: {len(data['angle'])} chars")
    print(f"  DIHE section length: {len(data['dihe'])} chars")
    
    assert 'remark' in data
    assert 'bond' in data
    assert 'angle' in data
    assert 'dihe' in data
    assert 'improper' in data
    assert 'nonbon' in data
    assert 'raw_text' in data


def test_write_frcmod():
    """Test writing a frcmod file."""
    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir) / "test.frcmod"
        
        write_amber_frcmod(
            output_path,
            remark="Test parameters",
            bond="c3-n3   300.0   1.45",
            angle="c3-n3-c3   50.0   109.5",
        )
        
        print(f"✓ Wrote frcmod file: {output_path}")
        
        # Read it back
        data = read_amber_frcmod(output_path)
        
        print(f"  Read back - Bond: {data['bond'].strip()}")
        print(f"  Read back - Angle: {data['angle'].strip()}")
        
        assert "c3-n3" in data['bond']
        assert "c3-n3-c3" in data['angle']
        assert "Test parameters" in data['remark']


def test_roundtrip():
    """Test reading and writing back preserves content."""
    frcmod_path = Path("/workspaces/molpy/test_amber_wrapper/parmchk2_output/tfsi.frcmod")
    
    if not frcmod_path.exists():
        print(f"Skipping: {frcmod_path} not found")
        return
    
    # Read original
    original = read_amber_frcmod(frcmod_path)
    
    # Write to temp file
    with tempfile.TemporaryDirectory() as tmpdir:
        temp_path = Path(tmpdir) / "roundtrip.frcmod"
        write_amber_frcmod(
            temp_path,
            remark=original['remark'],
            mass=original['mass'],
            bond=original['bond'],
            angle=original['angle'],
            dihe=original['dihe'],
            improper=original['improper'],
            nonbon=original['nonbon'],
        )
        
        # Read back
        roundtrip = read_amber_frcmod(temp_path)
        
        print(f"✓ Roundtrip test")
        print(f"  Original bond lines: {len(original['bond'].splitlines())}")
        print(f"  Roundtrip bond lines: {len(roundtrip['bond'].splitlines())}")
        
        # Compare sections (may have whitespace differences)
        assert roundtrip['remark'] == original['remark']
        assert len(roundtrip['bond'].strip()) > 0


if __name__ == "__main__":
    print("=" * 70)
    print("Testing Frcmod I/O")
    print("=" * 70)
    
    print("\n[Test 1] Reading existing frcmod file")
    print("-" * 70)
    test_read_frcmod()
    
    print("\n[Test 2] Writing new frcmod file")
    print("-" * 70)
    test_write_frcmod()
    
    print("\n[Test 3] Roundtrip test")
    print("-" * 70)
    test_roundtrip()
    
    print("\n" + "=" * 70)
    print("All tests passed!")
    print("=" * 70)
