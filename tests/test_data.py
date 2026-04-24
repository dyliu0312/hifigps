import os
import numpy as np
import pytest
from hifigps.data import (
    is_exist,
    get_file_dir,
    get_filename,
    delete_files,
    del_data,
    read_h5,
    save_h5,
    save_attrs,
    load_attrs,
    split_data_generator,
)

def test_data_manipulation(tmp_path):
    filepath = str(tmp_path / "test.h5")
    
    # test is_exist
    assert not is_exist(filepath)
    
    # generate test data
    np.random.seed(42)
    test_mask = np.random.choice([True, False], (10, 10, 10), p=[0.2, 0.8])
    
    # test save_h5
    save_h5(filepath, ["mask"], [test_mask])
    assert is_exist(filepath)
    
    # test read_h5
    data_read = read_h5(filepath, "mask")
    assert np.array_equal(data_read, test_mask)
    
    # test save_attrs and load_attrs
    file_attributes = {
        "author": "dyliu",
        "version": "0.1",
    }
    dataset_attrs = {
        "mask": {"units": "None", "description": "The boolen index for masking"}
    }
    save_attrs(filepath, file_attributes, dataset_attrs)
    
    attrs = load_attrs(filepath)
    assert attrs['file_attrs']['author'] == 'dyliu'
    assert attrs['dataset_attrs']['/mask']['units'] == 'None'
    
    # test split_data_generator
    gen = split_data_generator(5, test_mask, test_mask)
    for id1, id2 in gen:
        assert id1.shape[0] == 5
        assert id2.shape[0] == 5
    
    # test del_data
    del_data(filepath, "mask")
    # After deleting dataset, read_h5 should fail or return None/Empty 
    # (Checking if it still exists in file)
    with pytest.raises(Exception):
        read_h5(filepath, "mask")
        
    # test path functions
    assert get_filename(filepath) == "test"
    # tmp_path might have different format on different OS, but get_file_dir should work
    assert get_file_dir(filepath) == str(tmp_path)
    
    # test delete_files
    delete_files([filepath])
    assert not is_exist(filepath)
