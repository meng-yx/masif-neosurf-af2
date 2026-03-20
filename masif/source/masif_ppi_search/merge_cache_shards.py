import argparse
import os

import numpy as np

from cache_shard_common import FINAL_CACHE_FILES, ensure_dir, load_params, read_json


def parse_args():
    parser = argparse.ArgumentParser(description='Merge subset cache shards into final cache.')
    parser.add_argument('custom_params_module', nargs='?', default=None)
    parser.add_argument('--subset-root', default=None)
    parser.add_argument('--num-subsets', type=int, required=True)
    parser.add_argument('--cache-dir', default=None)
    return parser.parse_args()


def concat_or_empty(chunks, dtype=None):
    if len(chunks) == 0:
        return np.asarray([], dtype=dtype if dtype is not None else float)
    return np.concatenate(chunks, axis=0)


def load_legacy_subset_arrays(subset_dir):
    payload = {name: np.load(os.path.join(subset_dir, name), allow_pickle=True) for name in FINAL_CACHE_FILES}
    return payload


def load_parts_arrays(subset_dir, manifest):
    parts_dir = os.path.join(subset_dir, manifest.get('parts_dir', 'parts'))
    parts = manifest.get('parts', [])
    payload_lists = {name: [] for name in FINAL_CACHE_FILES}
    for part in parts:
        prefix = part['part_prefix']
        for name in FINAL_CACHE_FILES:
            fn = os.path.join(parts_dir, '{}_{}'.format(prefix, name))
            payload_lists[name].append(np.load(fn, allow_pickle=True))
    payload = {}
    for name, arrays in payload_lists.items():
        if 'names' in name:
            payload[name] = concat_or_empty(arrays, dtype=object)
        elif '_idx' in name:
            payload[name] = concat_or_empty(arrays, dtype=int).astype(int)
        else:
            payload[name] = concat_or_empty(arrays)
    return payload


def main():
    args = parse_args()
    params = load_params(args.custom_params_module)
    subset_root = args.subset_root or os.path.join(params['cache_dir'], 'subsets')
    cache_dir = args.cache_dir or params['cache_dir']
    ensure_dir(cache_dir)

    combined = {name: [] for name in FINAL_CACHE_FILES}
    pos_offset = 0
    neg_offset = 0
    catalog_fingerprint = None

    for subset_id in range(args.num_subsets):
        subset_dir = os.path.join(subset_root, 'subset_{}'.format(subset_id))
        success = os.path.join(subset_dir, '_SUCCESS')
        manifest_path = os.path.join(subset_dir, 'manifest.json')
        if not os.path.exists(success):
            raise RuntimeError('Subset not ready: missing {}'.format(success))
        if not os.path.exists(manifest_path):
            raise RuntimeError('Subset not ready: missing {}'.format(manifest_path))

        manifest = read_json(manifest_path)
        if catalog_fingerprint is None:
            catalog_fingerprint = manifest.get('catalog_fingerprint')
        elif manifest.get('catalog_fingerprint') != catalog_fingerprint:
            raise RuntimeError('Catalog fingerprint mismatch across subset manifests.')

        if int(manifest.get('format_version', 1)) >= 2:
            payload = load_parts_arrays(subset_dir, manifest)
        else:
            payload = load_legacy_subset_arrays(subset_dir)

        subset_pos_count = len(payload['pos_names.npy'])
        subset_neg_count = len(payload['neg_names.npy'])

        payload['pos_training_idx.npy'] = payload['pos_training_idx.npy'].astype(int) + pos_offset
        payload['pos_val_idx.npy'] = payload['pos_val_idx.npy'].astype(int) + pos_offset
        payload['pos_test_idx.npy'] = payload['pos_test_idx.npy'].astype(int) + pos_offset
        payload['neg_training_idx.npy'] = payload['neg_training_idx.npy'].astype(int) + neg_offset
        payload['neg_val_idx.npy'] = payload['neg_val_idx.npy'].astype(int) + neg_offset
        payload['neg_test_idx.npy'] = payload['neg_test_idx.npy'].astype(int) + neg_offset

        for name in FINAL_CACHE_FILES:
            combined[name].append(payload[name])

        pos_offset += subset_pos_count
        neg_offset += subset_neg_count

    merged = {}
    for name, arrays in combined.items():
        if 'names' in name:
            merged[name] = concat_or_empty(arrays, dtype=object)
        elif '_idx' in name:
            merged[name] = concat_or_empty(arrays, dtype=int).astype(int)
        else:
            merged[name] = concat_or_empty(arrays)

    if len(merged['pos_training_idx.npy']) == 0 or len(merged['neg_training_idx.npy']) == 0:
        raise RuntimeError('Insufficient train samples after merge.')
    if len(merged['pos_names.npy']) != len(merged['binder_rho_wrt_center.npy']):
        raise RuntimeError('Positive binder length mismatch after merge.')
    if len(merged['pos_names.npy']) != len(merged['pos_rho_wrt_center.npy']):
        raise RuntimeError('Positive target length mismatch after merge.')
    if len(merged['neg_names.npy']) != len(merged['neg_rho_wrt_center.npy']):
        raise RuntimeError('Negative length mismatch after merge.')

    assert len(set(merged['pos_training_idx.npy']) & set(merged['pos_val_idx.npy'])) == 0
    assert len(set(merged['pos_training_idx.npy']) & set(merged['pos_test_idx.npy'])) == 0
    assert len(set(merged['pos_val_idx.npy']) & set(merged['pos_test_idx.npy'])) == 0
    assert len(set(merged['neg_training_idx.npy']) & set(merged['neg_val_idx.npy'])) == 0
    assert len(set(merged['neg_training_idx.npy']) & set(merged['neg_test_idx.npy'])) == 0
    assert len(set(merged['neg_val_idx.npy']) & set(merged['neg_test_idx.npy'])) == 0

    for name in FINAL_CACHE_FILES:
        np.save(os.path.join(cache_dir, name), merged[name])

    print('Merged subsets into {} with pos={} neg={}'.format(cache_dir, len(merged['pos_names.npy']), len(merged['neg_names.npy'])))


if __name__ == '__main__':
    main()
