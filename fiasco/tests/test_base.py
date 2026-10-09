"""
Test base ion functionality
"""
import h5py
import logging
import numpy as np
import os
import pytest

from packaging.version import Version

import fiasco

from fiasco.util import build_hdf5_dbase
from fiasco.util.exceptions import MissingIonError


@pytest.mark.parametrize('ion_name', [
    "fe 21",
    "Fe 21",
    "Iron 21",
    "iron XXI",
    "Fe xxi",
    "Fe 20+",
    (26, 21),
    ("Fe", "XXI"),
])
def test_create_ion_input_formats(hdf5_dbase_root, ion_name):
    ion = fiasco.base.IonBase(ion_name, hdf5_dbase_root=hdf5_dbase_root)
    assert ion.element_name == 'iron'
    assert ion.atomic_symbol == 'Fe'
    assert ion.atomic_number == 26
    assert ion.ionization_stage == 21
    assert ion.charge_state == 20
    assert ion.isoelectronic_sequence == 'C'
    assert ion._ion_name == 'fe_21'
    assert ion.ion_name_roman == 'Fe XXI'
    assert ion.ionization_stage_roman == 'XXI'
    assert ion.period == 4


def test_create_invalid_ion_raises_missing_ion_error(hdf5_dbase_root):
    with pytest.raises(MissingIonError):
        fiasco.base.IonBase('hydrogen 3', hdf5_dbase_root=hdf5_dbase_root)


@pytest.mark.parametrize(('ion_name', 'sequence_name'), [
    ("Fe 21", "C"),
    ("Fe 26", "H"),
    ("Fe 25", "He"),
    ("C 7", None),
])
def test_isoelectronic_sequence(ion_name, sequence_name, hdf5_dbase_root):
    ion = fiasco.base.IonBase(ion_name, hdf5_dbase_root=hdf5_dbase_root)
    assert ion.isoelectronic_sequence == sequence_name


def test_old_fiasco_version_warning(tmpdir_factory, caplog):
    path = tmpdir_factory.mktemp('fiasco').join('chianti_dbase_empty_old_version.h5')
    build_hdf5_dbase(None, path, files=[], check_hash=False, overwrite=True)
    test_ion = 'H 1'
    with h5py.File(path, mode='a') as hf:
        current_version = Version(hf.attrs['fiasco_version'])
        old_version = Version(
            f'{current_version.major}.{current_version.minor-1}.{current_version.micro}',
        )
        hf.attrs['fiasco_version'] = str(old_version)
        del hf['ion_index']
        ds = hf.create_dataset('ion_index', data=np.array([test_ion,]).astype(np.bytes_))
        ds.attrs['unit'] = 'SKIP'
    with caplog.at_level(logging.WARN):
        fiasco.base.IonBase(test_ion, hdf5_dbase_root=path)
    if current_version.is_devrelease:
        assert caplog.text == ''
    else:
        assert 'was produced with an earlier version of fiasco' in caplog.text


@pytest.mark.parametrize(('second_ion', 'second_version'), [
    ('C 1', '0.9.0'),
    ('H 1', '0.8.0'),
])
def test_database_cache_with_relative_paths(tmp_path, monkeypatch, caplog, second_ion, second_version):
    monkeypatch.setattr(fiasco, '__version__', '0.9.0')
    paths = []
    for directory, ion_name, version in [('first', 'H 1', '0.9.0'), ('second', second_ion, second_version)]:
        path = tmp_path / directory / 'dbase.h5'
        path.parent.mkdir()
        with h5py.File(path, 'w') as hf:
            hf.attrs['fiasco_version'] = version
            dataset = hf.create_dataset('ion_index', data=np.array([ion_name], dtype='S3'))
            dataset.attrs['unit'] = 'SKIP'
        os.utime(path, ns=(1_000_000_000, 1_000_000_000))
        paths.append(path)
    assert paths[0].stat().st_size == paths[1].stat().st_size
    monkeypatch.chdir(paths[0].parent)
    fiasco.base.IonBase('H 1', hdf5_dbase_root='dbase.h5')
    monkeypatch.chdir(paths[1].parent)
    caplog.clear()
    with caplog.at_level(logging.WARN):
        ion = fiasco.base.IonBase(second_ion, hdf5_dbase_root='dbase.h5')
    assert ion.ion_name == second_ion
    assert ('was produced with an earlier version' in caplog.text) == (second_version == '0.8.0')
