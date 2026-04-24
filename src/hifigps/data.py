"""
tools for HDF5 data manipulation.

"""

import os
from typing import Any, Callable, Dict, Iterable, List, Optional, Union  # pyright: ignore[reportDeprecated]

import h5py as h5
import numpy as np


# --- os ---
def get_file_dir(path: str) -> str:
    """
    get the file parrent dir
    """
    return os.path.dirname(path)


def get_filename(path: str) -> str:
    """
    get the file name
    """
    file_name = os.path.splitext(os.path.basename(path))[0]
    return file_name


def delete_files(file_paths: List[str]) -> None:
    """
    Deletes files from the given list of file paths.

    Args:
        file_paths (List[str]): A list of file paths to be deleted.

    Prints:
        Messages indicating whether each file was deleted, not found, or caused an error.
    """
    for file_path in file_paths:
        try:
            os.remove(file_path)
            print(f"Deleted: {file_path}")
        except FileNotFoundError:
            print(f"File not found: {file_path}")
        except PermissionError:
            print(f"Permission error: {file_path}")
        except Exception as e:
            print(f"Error deleting {file_path}: {e}")


def is_exist(file_path: str) -> bool:
    """
    check file exists.
    """
    return os.path.exists(file_path)


# --- hdf5 file ---
def info_items(name: str, item: Any):
    """
    Prints detailed information (type, shape, dtype, and all attributes) for a
    Group or Dataset object encountered during h5py.File.visititems().

    This function is designed to be passed as the `func` argument to visititems().

    Args:
        name (str): The full path name of the item (Group or Dataset) in the file.
        item (h5py.Group or h5py.Dataset): The h5py object itself.
    """
    print("=" * 40)
    print(f"Item: {name}")

    # 1. Print item type
    item_type = (
        "Group"
        if isinstance(item, h5.Group)
        else "Dataset"
        if isinstance(item, h5.Dataset)
        else "Other"
    )
    print(f"Type: {item_type}")

    # 2. Print shape and dtype for Datasets
    if isinstance(item, h5.Dataset):
        print(f"Shape: {item.shape}")
        print(f"Dtype: {item.dtype}")

    # 3. Print all attributes
    print(f"Attributes ({len(item.attrs)}):")
    if item.attrs:
        for attr_name, value in item.attrs.items():
            # Post-process attribute value for clean display
            display_value = value

            # Decode byte strings (often used for HDF5 strings)
            if isinstance(value, bytes):
                display_value = value.decode("utf-8")
            # Handle numpy arrays of byte strings
            elif isinstance(value, np.ndarray) and value.dtype.kind == "S":
                display_value = value.astype(str)
            # Limit array display length
            elif isinstance(value, np.ndarray):
                # For very long arrays, only show the start and type
                if value.size > 5:
                    display_value = (
                        f"np.array(shape={value.shape}, dtype={value.dtype})"
                    )
                # Otherwise, show the array

            # Print the attribute
            print(f"- {attr_name}: {display_value}")
    else:
        print(" None")


def vitems(
    file_path: str,
    func: Callable[[str, Union[h5.Group, h5.Dataset]], None] = info_items,
):
    """
    Recursively visit the group and dataset information in an HDF5 file using the visititems() method.

    Args:
        file_path (str): Path to the HDF5 file.
        func (callable): Function to be applied to each item (group or dataset) in the HDF5 file.
                        Defaults to `print` which will print the name and item type.
    """
    try:
        # Open the HDF5 file in read-only mode
        with h5.File(file_path, "r") as f:
            # Visit all items in the file (groups and datasets) and apply `func` to each
            f.visititems(func)
    except Exception as e:
        print(f"Error accessing HDF5 file: {e}")


def del_data(file_path: str, key: str):
    """
    delete dataset and it's attributes from h5 file.
    """
    with h5.File(file_path, "r+") as f:
        # del attributes
        for k in f[key].attrs.keys():
            del f[key].attrs[k]
        # del data
        del f[key]


def print_dset_attrs(dataset: h5.Dataset):
    """
    Prints the attributes of a given dataset if any exist.
    """
    attrs = dataset.attrs
    if attrs:
        print(f"Attributes for '{dataset.name}':")
        for attr_name, attr_value in attrs.items():
            print(f"  {attr_name}: {attr_value}")


def read_h5(
    file_path: str, key: Union[str, Iterable[str]], print_attrs: bool = False
) -> Union[np.ndarray, List[np.ndarray]]:
    """
    Read dataset from an HDF5 file.

    Args:
        file_path (str): Path to the HDF5 file.
        key (str or list): Key(s) to read from the HDF5 file.
            It can be a single key (str) or a list of keys.
        print_attrs (bool): If True, prints the attributes of the datasets
            when reading them.

    Returns:
        out (ndarray or list): The dataset(s) read from the HDF5 file.
    """
    # Ensure the key is either a string or a list of strings
    if not isinstance(key, (str, list)):
        raise TypeError("Key must be a string or a list of strings.")

    # Open the HDF5 file in read-only mode
    with h5.File(file_path, "r") as f:
        available_keys = list(f.keys())

        # Helper function to handle reading and printing attributes for each key
        def get_data_for_key(k: str) -> np.ndarray:
            if k not in f:
                raise ValueError(
                    f"Key '{k}' not found. Available keys are: {available_keys}"
                )
            obj = f[k]
            if not isinstance(obj, h5.Dataset):
                raise KeyError(f"Key '{k}' is a group, not a dataset.")
            else:
                if k in f and not obj.shape == ():
                    dataset = obj[:]
                else:
                    dataset = obj[()]
                if print_attrs:
                    print_dset_attrs(obj)
                return dataset

        # If a single key is provided, return the corresponding data
        if isinstance(key, str):
            return get_data_for_key(key)

        # If a list of keys is provided, return a list of corresponding data
        elif isinstance(key, list):
            data_list = [get_data_for_key(k) for k in key]
            return data_list


def save_h5(
    output: str,
    keys: List[str],
    data: List[np.ndarray],
    group_name: Union[str, None] = None,
    **kwargs,
) -> None:
    """
    Save data to an HDF5 file and create group if group_name is provided.

    **For single dataset, keys and datas should be provided list with one element.
    e.g. keys = ['data',], data = [data,]**

    Args:
        output (str): The path to the output HDF5 file.
        keys (list): A list of keys for the datasets.
        data (list): A list of data arrays to be saved.
        group_name (str, optional): A group paths to be created. Defaults to None.
        kwargs (dict, optional): Additional arguments for `f.create_dataset`.

    """
    with h5.File(output, "a") as f:
        if group_name is not None:
            f.create_group(group_name)

        for k, d in zip(keys, data):
            if group_name is not None:
                dataset_path = f"{group_name}/{k}"
            else:
                dataset_path = k
            f.create_dataset(dataset_path, data=d, **kwargs)


def save_attrs(filepath, file_attrs=None, dataset_attrs=None):
    """
    Save specified attributes in an HDF5 file, with an option to add dataset-level attributes.

    Parameters:
    filepath (str): Path to the HDF5 file. It will create the file if it doesn't exist.
    attributes (dict, optional): A dictionary containing file-level key-value pairs to store as attributes.
    dataset_attrs (dict, optional): A dictionary containing keys as dataset names and values as dictionaries of dataset-level attributes.
                                     For example: {'dataset1': {'units': 'm', 'description': 'Height data'}}.
    """

    if file_attrs is None and dataset_attrs is None:
        raise ValueError("Please provided a kind of attributes at least.")

    with h5.File(filepath, "a") as f:
        if file_attrs is not None:
            # Store file-level attributes in the root of the HDF5 file
            for key, value in file_attrs.items():  # type: ignore
                try:
                    f.attrs[key] = value
                except Exception as e:
                    print(f"Error storing attributes: {e}")
            print(f"Successfully stored file-level attributes in {filepath}")

        if dataset_attrs is not None:
            # Store dataset-specific attributes (if provided)
            for dataset_name, attrs in dataset_attrs.items():
                if dataset_name not in f:
                    raise ValueError(
                        f"Dataset '{dataset_name}' not found in {filepath}"
                    )
                dataset = f[dataset_name]
                for attr_name, attr_value in attrs.items():
                    dataset.attrs[attr_name] = attr_value
                print(
                    f"Successfully stored dataset-level attributes for '{dataset_name}'"
                )


def _decode_attribute_value(value: Any) -> Any:
    """Helper function to decode HDF5 byte strings."""
    if isinstance(value, bytes):
        return value.decode("utf-8")
    elif isinstance(value, np.ndarray) and value.dtype.kind == "S":
        # Decode numpy arrays of byte strings
        return value.astype(str)
    return value


def load_attrs(file_path: str) -> Optional[Dict[str, Any]]:
    """
    Load attributes from an HDF5 file, structured to mirror the save_h5_attributes function.

    It returns file-level attributes and a nested dictionary of dataset-level attributes.

    Args:
        file_path (str): The path to the HDF5 file.

    Returns:
        Optional[Dict[str, Any]]: A dictionary with 'file_attributes' and 'dataset_attrs' keys,
                                  or None if the file cannot be opened.
                                  Example structure:
                                  {
                                      'file_attributes': {'attr1': val1, 'attr2': val2},
                                      'dataset_attrs': {
                                          '/dataset1': {'d_attr1': d_val1},
                                          '/group/dataset2': {'d_attr2': d_val2}
                                      }
                                  }
    """
    try:
        with h5.File(file_path, "r") as f:
            # 1. Read File-Level Attributes
            # These are attributes stored directly on the root object '/'.
            file_attrs = {}
            for name, value in f.attrs.items():
                file_attrs[name] = _decode_attribute_value(value)

            # 2. Read Dataset-Level Attributes
            # We use f.visititems to find all datasets and read their attributes.
            dataset_attrs = {}

            def visitor_func(name, item):
                """Visitor function called by f.visititems()."""
                # Check if the item is a Dataset (and not a Group)
                if isinstance(item, h5.Dataset):
                    attrs = {}
                    if item.attrs:
                        for attr_name, attr_value in item.attrs.items():
                            attrs[attr_name] = _decode_attribute_value(attr_value)

                        # Store the attributes using the full path name as the key
                        dataset_attrs[f"/{name}"] = attrs

            # Traverse the file structure
            f.visititems(visitor_func)

            # Return the combined, structured result
            return {"file_attrs": file_attrs, "dataset_attrs": dataset_attrs}

    except FileNotFoundError:
        print(f"Error: File not found at {file_path}")
        return None
    except Exception as e:
        print(f"An error occurred while loading HDF5 attributes: {e}")
        return None


# --- generator ---
def split_data_generator(splitsize: int, *data: np.ndarray):
    """
    Splits multiple datasets into chunks of size `splitsize`, ensuring that all datasets have the same
    first dimension size. The generator yields corresponding slices for each dataset.

    :param data: Variable number of datasets (arrays).
    :param splitsize: The size of each split (chunk).
    :return: A generator yielding tuples of slices from each dataset.
    """
    # Check that all input datasets have the same first dimension (number of samples)
    num_samples = data[0].shape[0]
    if not all(d.shape[0] == num_samples for d in data):
        raise ValueError(
            "All input datasets must have the same length along the first dimension."
        )

    num_splits = num_samples // splitsize
    remainder = num_samples % splitsize

    # Generate splits for the integer multiple parts
    for i in range(num_splits):
        start_idx = i * splitsize
        end_idx = (i + 1) * splitsize
        yield tuple(d[start_idx:end_idx] for d in data)

    # Generate split for the remainder part, if any
    if remainder > 0:
        yield tuple(d[num_splits * splitsize :] for d in data)


# --- stack result---
def get_stacked_result(
    file_path: str, s_key: str = "Signal", m_key: str = "Mask"
) -> np.ndarray:
    """
    Add up and average to get the stacked result from the HDF5 file.
    Those datasets are seperately saved in groups with keys "Signal" and "Mask", e.g. "/0_500/Signal", "/0_500/Mask".

    Args:
        file_path (str): Path to the HDF5 file.
        s_key (str): Key for the signal dataset. Defaults to "Signal".
        m_key (str): Key for the mask dataset. Defaults to "Mask".

    Returns:
        result_array (np.ndarray): The stacked result (3d array).
    """
    result = []

    def add_up(name, obj):
        if isinstance(obj, h5.Group):
            si = obj[s_key][:]  # pyright: ignore[reportIndexIssue]
            mi = obj[m_key][:]  # pyright: ignore[reportIndexIssue]
            s = np.ma.array(si, mask=mi)
            result.append(s)

    vitems(file_path, add_up)
    result_array = np.ma.array(result)
    return result_array.mean(axis=0)
