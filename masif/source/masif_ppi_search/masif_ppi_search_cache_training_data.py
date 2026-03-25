import argparse
import os
from collections import OrderedDict

import numpy as np

from cache_shard_common import (
    FINAL_CACHE_FILES,
    concat_or_empty,
    ensure_dir,
    load_params,
    normalize_mix,
    read_json,
    remap_indices,
    save_json,
    split_counts,
    write_success_marker,
)


FEATURE_TENSOR_FILES = {
    'p1_rho': 'p1_rho_wrt_center.npy',
    'p1_theta': 'p1_theta_wrt_center.npy',
    'p1_input': 'p1_input_feat.npy',
    'p1_mask': 'p1_mask.npy',
    'p2_rho': 'p2_rho_wrt_center.npy',
    'p2_theta': 'p2_theta_wrt_center.npy',
    'p2_input': 'p2_input_feat.npy',
    'p2_mask': 'p2_mask.npy',
}


class PairArrayCache:
    def __init__(self, max_records=16, mmap_mode='r'):
        self.max_records = max(1, int(max_records))
        self.mmap_mode = mmap_mode
        self._store = OrderedDict()

    def get(self, record):
        rec_key = record['ppi_pair_id']
        if rec_key in self._store:
            self._store.move_to_end(rec_key)
            return self._store[rec_key]

        in_dir = record['in_dir']
        arrays = {
            name: np.load(os.path.join(in_dir, filename), mmap_mode=self.mmap_mode)
            for name, filename in FEATURE_TENSOR_FILES.items()
        }
        self._store[rec_key] = arrays
        self._store.move_to_end(rec_key)
        while len(self._store) > self.max_records:
            self._store.popitem(last=False)
        return arrays


class SubsetPartWriter:
    def __init__(self, subset_dir, flush_every_pos):
        self.subset_dir = subset_dir
        self.parts_dir = os.path.join(subset_dir, 'parts')
        ensure_dir(self.parts_dir)
        self.flush_every_pos = max(1, int(flush_every_pos))

        self.part_idx = 0
        self.global_pos_count = 0
        self.global_neg_count = 0
        self.total_mix_counts = {'cross': 0, 'within': 0, 'hard': 0, 'fallback': 0}
        self.parts = []
        self._reset_buffer()

    def _reset_buffer(self):
        self.binder_rho_wrt_center = []
        self.binder_theta_wrt_center = []
        self.binder_input_feat = []
        self.binder_mask = []
        self.pos_rho_wrt_center = []
        self.pos_theta_wrt_center = []
        self.pos_input_feat = []
        self.pos_mask = []
        self.neg_rho_wrt_center = []
        self.neg_theta_wrt_center = []
        self.neg_input_feat = []
        self.neg_mask = []
        self.pos_names = []
        self.neg_names = []

        self.pos_training_idx = []
        self.pos_val_idx = []
        self.pos_test_idx = []
        self.neg_training_idx = []
        self.neg_val_idx = []
        self.neg_test_idx = []

        self.local_pos_count = 0
        self.local_neg_count = 0

    def add_pos_block(self, split, binder_rows, pos_rows, pos_names):
        n_pos = len(pos_names)
        self.binder_rho_wrt_center.append(binder_rows[0])
        self.binder_theta_wrt_center.append(binder_rows[1])
        self.binder_input_feat.append(binder_rows[2])
        self.binder_mask.append(binder_rows[3])
        self.pos_rho_wrt_center.append(pos_rows[0])
        self.pos_theta_wrt_center.append(pos_rows[1])
        self.pos_input_feat.append(pos_rows[2])
        self.pos_mask.append(pos_rows[3])
        self.pos_names.extend(pos_names)

        pos_indices = np.arange(self.local_pos_count, self.local_pos_count + n_pos)
        if split == 'train':
            self.pos_training_idx.extend(pos_indices.tolist())
        elif split == 'val':
            self.pos_val_idx.extend(pos_indices.tolist())
        else:
            self.pos_test_idx.extend(pos_indices.tolist())
        self.local_pos_count += n_pos

    def add_neg_block(self, split, neg_rows, neg_names):
        n_neg = len(neg_names)
        if n_neg == 0:
            return
        self.neg_rho_wrt_center.append(neg_rows[0])
        self.neg_theta_wrt_center.append(neg_rows[1])
        self.neg_input_feat.append(neg_rows[2])
        self.neg_mask.append(neg_rows[3])
        self.neg_names.extend(neg_names)

        neg_indices = np.arange(self.local_neg_count, self.local_neg_count + n_neg)
        if split == 'train':
            self.neg_training_idx.extend(neg_indices.tolist())
        elif split == 'val':
            self.neg_val_idx.extend(neg_indices.tolist())
        else:
            self.neg_test_idx.extend(neg_indices.tolist())
        self.local_neg_count += n_neg

    def accumulate_mix_counts(self, mix_counts):
        for key in self.total_mix_counts:
            self.total_mix_counts[key] += int(mix_counts.get(key, 0))

    def should_flush(self):
        return self.local_pos_count >= self.flush_every_pos

    def _save_part_array(self, part_prefix, filename, array):
        np.save(os.path.join(self.parts_dir, '{}_{}'.format(part_prefix, filename)), array)

    def flush(self):
        if self.local_pos_count == 0 and self.local_neg_count == 0:
            return

        binder_rho_wrt_center = concat_or_empty(self.binder_rho_wrt_center, 2)
        binder_theta_wrt_center = concat_or_empty(self.binder_theta_wrt_center, 2)
        binder_input_feat = concat_or_empty(self.binder_input_feat, 3)
        binder_mask = concat_or_empty(self.binder_mask, 2)
        pos_rho_wrt_center = concat_or_empty(self.pos_rho_wrt_center, 2)
        pos_theta_wrt_center = concat_or_empty(self.pos_theta_wrt_center, 2)
        pos_input_feat = concat_or_empty(self.pos_input_feat, 3)
        pos_mask = concat_or_empty(self.pos_mask, 2)
        neg_rho_wrt_center = concat_or_empty(self.neg_rho_wrt_center, 2)
        neg_theta_wrt_center = concat_or_empty(self.neg_theta_wrt_center, 2)
        neg_input_feat = concat_or_empty(self.neg_input_feat, 3)
        neg_mask = concat_or_empty(self.neg_mask, 2)

        pos_not_nan = ~np.isnan(binder_input_feat).any(axis=(1, 2))
        pos_not_nan = pos_not_nan & (~np.isnan(pos_input_feat).any(axis=(1, 2)))
        neg_not_nan = ~np.isnan(neg_input_feat).any(axis=(1, 2))

        pos_names = np.asarray([x for i, x in enumerate(self.pos_names) if pos_not_nan[i]], dtype=object)
        neg_names = np.asarray([x for i, x in enumerate(self.neg_names) if neg_not_nan[i]], dtype=object)

        binder_rho_wrt_center = binder_rho_wrt_center[pos_not_nan, :]
        binder_theta_wrt_center = binder_theta_wrt_center[pos_not_nan, :]
        binder_input_feat = binder_input_feat[pos_not_nan, ...]
        binder_mask = binder_mask[pos_not_nan, :]
        pos_rho_wrt_center = pos_rho_wrt_center[pos_not_nan, :]
        pos_theta_wrt_center = pos_theta_wrt_center[pos_not_nan, :]
        pos_input_feat = pos_input_feat[pos_not_nan, ...]
        pos_mask = pos_mask[pos_not_nan, :]
        neg_rho_wrt_center = neg_rho_wrt_center[neg_not_nan, :]
        neg_theta_wrt_center = neg_theta_wrt_center[neg_not_nan, :]
        neg_input_feat = neg_input_feat[neg_not_nan, ...]
        neg_mask = neg_mask[neg_not_nan, :]

        pos_training_idx = remap_indices(self.pos_training_idx, pos_not_nan)
        pos_val_idx = remap_indices(self.pos_val_idx, pos_not_nan)
        pos_test_idx = remap_indices(self.pos_test_idx, pos_not_nan)
        neg_training_idx = remap_indices(self.neg_training_idx, neg_not_nan)
        neg_val_idx = remap_indices(self.neg_val_idx, neg_not_nan)
        neg_test_idx = remap_indices(self.neg_test_idx, neg_not_nan)

        assert len(set(pos_training_idx) & set(pos_val_idx)) == 0
        assert len(set(pos_training_idx) & set(pos_test_idx)) == 0
        assert len(set(pos_val_idx) & set(pos_test_idx)) == 0
        assert len(set(neg_training_idx) & set(neg_val_idx)) == 0
        assert len(set(neg_training_idx) & set(neg_test_idx)) == 0
        assert len(set(neg_val_idx) & set(neg_test_idx)) == 0

        pos_base = self.global_pos_count
        neg_base = self.global_neg_count

        pos_training_idx = np.asarray(pos_training_idx, dtype=int) + pos_base
        pos_val_idx = np.asarray(pos_val_idx, dtype=int) + pos_base
        pos_test_idx = np.asarray(pos_test_idx, dtype=int) + pos_base
        neg_training_idx = np.asarray(neg_training_idx, dtype=int) + neg_base
        neg_val_idx = np.asarray(neg_val_idx, dtype=int) + neg_base
        neg_test_idx = np.asarray(neg_test_idx, dtype=int) + neg_base

        part_prefix = 'part_{:05d}'.format(self.part_idx)
        payload = {
            'pos_names.npy': pos_names,
            'neg_names.npy': neg_names,
            'binder_rho_wrt_center.npy': binder_rho_wrt_center,
            'binder_theta_wrt_center.npy': binder_theta_wrt_center,
            'binder_input_feat.npy': binder_input_feat,
            'binder_mask.npy': binder_mask,
            'pos_training_idx.npy': pos_training_idx,
            'pos_val_idx.npy': pos_val_idx,
            'pos_test_idx.npy': pos_test_idx,
            'pos_rho_wrt_center.npy': pos_rho_wrt_center,
            'pos_theta_wrt_center.npy': pos_theta_wrt_center,
            'pos_input_feat.npy': pos_input_feat,
            'pos_mask.npy': pos_mask,
            'neg_training_idx.npy': neg_training_idx,
            'neg_val_idx.npy': neg_val_idx,
            'neg_test_idx.npy': neg_test_idx,
            'neg_rho_wrt_center.npy': neg_rho_wrt_center,
            'neg_theta_wrt_center.npy': neg_theta_wrt_center,
            'neg_input_feat.npy': neg_input_feat,
            'neg_mask.npy': neg_mask,
        }
        for filename in FINAL_CACHE_FILES:
            self._save_part_array(part_prefix, filename, payload[filename])

        new_pos = int(len(pos_names))
        new_neg = int(len(neg_names))
        self.global_pos_count += new_pos
        self.global_neg_count += new_neg
        self.parts.append(
            {
                'part_idx': self.part_idx,
                'part_prefix': part_prefix,
                'pos_count': new_pos,
                'neg_count': new_neg,
                'pos_base': int(pos_base),
                'neg_base': int(neg_base),
            }
        )
        self.part_idx += 1
        self._reset_buffer()

    def finalize(self, subset_manifest):
        self.flush()
        manifest = dict(subset_manifest)
        manifest.update(
            {
                'format_version': 2,
                'parts_dir': 'parts',
                'parts_count': len(self.parts),
                'parts': self.parts,
                'pos_count': int(self.global_pos_count),
                'neg_count': int(self.global_neg_count),
                'mix_counts': self.total_mix_counts,
            }
        )
        save_json(os.path.join(self.subset_dir, 'manifest.json'), manifest)
        write_success_marker(self.subset_dir)
        return manifest


def parse_args():
    parser = argparse.ArgumentParser(description='Generate one cache subset from a global catalog.')
    parser.add_argument('custom_params_module', nargs='?', default=None)
    parser.add_argument('--subset-id', type=int, required=True)
    parser.add_argument('--num-subsets', type=int, required=True)
    parser.add_argument('--progress-every-records', type=int, default=1)
    parser.add_argument('--pair-cache-max-records', type=int, default=16)
    parser.add_argument('--flush-every-pos', type=int, default=256)
    return parser.parse_args()


def load_catalog(catalog_dir):
    success = os.path.join(catalog_dir, '_SUCCESS')
    if not os.path.exists(success):
        raise RuntimeError('Catalog not ready: missing {}'.format(success))

    manifest = read_json(os.path.join(catalog_dir, 'manifest.json'))
    records = np.load(os.path.join(catalog_dir, 'records.npy'), allow_pickle=True)

    pools = {}
    for key in ['trainval', 'test']:
        pools[key] = {
            'source_record_idx': np.load(
                os.path.join(catalog_dir, '{}_source_record_idx.npy'.format(key)),
                mmap_mode='r',
            ),
            'local_within_idx': np.load(
                os.path.join(catalog_dir, '{}_local_within_idx.npy'.format(key)),
                mmap_mode='r',
            ),
        }
    return manifest, records, pools


def records_to_lookup_arrays(records):
    ppi = []
    pdb = []
    for rec in records:
        item = rec.item() if hasattr(rec, 'item') else rec
        ppi.append(item['ppi_pair_id'])
        pdb.append(item['pdb_id'])
    return np.asarray(ppi, dtype=object), np.asarray(pdb, dtype=object)


def get_record(records, idx):
    rec = records[idx]
    return rec.item() if hasattr(rec, 'item') else rec


def get_cross_rows_and_names(pool, picks, records, pair_cache):
    rho_chunks = []
    theta_chunks = []
    input_chunks = []
    mask_chunks = []
    names = []

    for pick in picks:
        src_rec_idx = int(pool['source_record_idx'][pick])
        local_within_idx = int(pool['local_within_idx'][pick])
        src_rec = get_record(records, src_rec_idx)
        src_arrays = pair_cache.get(src_rec)
        src_k_neg2 = np.asarray(src_rec['k_neg2'], dtype=int)
        absolute_idx = int(src_k_neg2[local_within_idx])

        rho_chunks.append(src_arrays['p2_rho'][absolute_idx : absolute_idx + 1])
        theta_chunks.append(src_arrays['p2_theta'][absolute_idx : absolute_idx + 1])
        input_chunks.append(src_arrays['p2_input'][absolute_idx : absolute_idx + 1])
        mask_chunks.append(src_arrays['p2_mask'][absolute_idx : absolute_idx + 1])
        names.append('{}_p2_{}'.format(src_rec['ppi_pair_id'], absolute_idx))

    return rho_chunks, theta_chunks, input_chunks, mask_chunks, names


def draw_without_replacement_cycle(state, n_samples):
    if n_samples <= 0:
        return np.asarray([], dtype=int)

    base_idx = state['base_idx']
    if len(base_idx) == 0:
        return np.asarray([], dtype=int)

    picks = []
    while len(picks) < n_samples:
        perm = state['perm']
        cursor = state['cursor']
        remaining = len(perm) - cursor
        if remaining <= 0:
            state['perm'] = np.random.permutation(base_idx)
            state['cursor'] = 0
            continue

        take = min(n_samples - len(picks), remaining)
        picked = perm[cursor : cursor + take]
        if take == 1:
            picks.append(int(picked[0]))
        else:
            picks.extend(picked.astype(int).tolist())
        state['cursor'] += take

    return np.asarray(picks, dtype=int)


def main():
    args = parse_args()
    if args.subset_id < 0 or args.subset_id >= args.num_subsets:
        raise ValueError('--subset-id must be in [0, num_subsets).')

    params = load_params(args.custom_params_module)
    # Keep all cache artifacts under params["cache_dir"] (custom_params.py).
    catalog_dir = os.path.join(params['cache_dir'], 'catalog')
    subset_root = os.path.join(params['cache_dir'], 'subsets')
    subset_dir = os.path.join(subset_root, 'subset_{}'.format(args.subset_id))
    ensure_dir(subset_dir)

    catalog_manifest, records, pools = load_catalog(catalog_dir)
    record_ppi, record_pdb = records_to_lookup_arrays(records)
    assigned_record_indices = [i for i in range(len(records)) if i % args.num_subsets == args.subset_id]

    neg_ratio = max(1, int(params.get('neg_ratio', 1)))
    mix = normalize_mix(
        params.get('neg_mix_cross_complex', 0.0),
        params.get('neg_mix_within_complex', 1.0),
        params.get('neg_mix_hard', 0.0),
    )
    enforce_diff_pdb_for_cross = bool(params.get('enforce_diff_pdb_for_cross', True))
    hard_negative_topk = int(params.get('hard_negative_topk', 200))

    print('Subset job config: subset_id={} num_subsets={}'.format(args.subset_id, args.num_subsets))
    print('Assigned records: {}'.format(len(assigned_record_indices)))
    print('pair_cache_max_records={}'.format(args.pair_cache_max_records))
    print('flush_every_pos={}'.format(args.flush_every_pos))

    pair_cache = PairArrayCache(max_records=args.pair_cache_max_records, mmap_mode='r')
    writer = SubsetPartWriter(subset_dir=subset_dir, flush_every_pos=args.flush_every_pos)

    for local_count, rec_idx in enumerate(assigned_record_indices):
        if local_count % max(1, args.progress_every_records) == 0:
            print('Record {}/{} (global idx {})'.format(local_count, len(assigned_record_indices), rec_idx))

        rec = get_record(records, rec_idx)
        rec_arrays = pair_cache.get(rec)

        k1 = np.asarray(rec['k1'], dtype=int)
        k2 = np.asarray(rec['k2'], dtype=int)
        k_neg2 = np.asarray(rec['k_neg2'], dtype=int)
        within_dneg = np.asarray(rec['within_dneg'], dtype=float)

        n_pos = len(k1)
        if n_pos == 0:
            continue

        writer.add_pos_block(
            split=rec['split'],
            binder_rows=(
                rec_arrays['p1_rho'][k1],
                rec_arrays['p1_theta'][k1],
                rec_arrays['p1_input'][k1],
                rec_arrays['p1_mask'][k1],
            ),
            pos_rows=(
                rec_arrays['p2_rho'][k2],
                rec_arrays['p2_theta'][k2],
                rec_arrays['p2_input'][k2],
                rec_arrays['p2_mask'][k2],
            ),
            pos_names=['{}_p1_{}'.format(rec['ppi_pair_id'], int(ii)) for ii in k1],
        )

        pool_key = 'test' if rec['split'] == 'test' else 'trainval'
        pool = pools[pool_key]
        cross_all_idx = np.arange(len(pool['source_record_idx']))

        if len(cross_all_idx) == 0:
            cross_valid_idx = np.asarray([], dtype=int)
        else:
            pool_src = np.asarray(pool['source_record_idx'])
            mask = record_ppi[pool_src] != rec['ppi_pair_id']
            if enforce_diff_pdb_for_cross:
                mask = mask & (record_pdb[pool_src] != rec['pdb_id'])
            cross_valid_idx = cross_all_idx[mask]

        within_all_idx = np.arange(len(k_neg2))
        hard_order = np.argsort(within_dneg)
        hard_idx = hard_order[: min(hard_negative_topk, len(hard_order))]
        cross_valid_state = {
            'base_idx': cross_valid_idx,
            'perm': np.random.permutation(cross_valid_idx),
            'cursor': 0,
        }
        cross_all_state = {
            'base_idx': cross_all_idx,
            'perm': np.random.permutation(cross_all_idx),
            'cursor': 0,
        }
        within_state = {
            'base_idx': within_all_idx,
            'perm': np.random.permutation(within_all_idx),
            'cursor': 0,
        }
        hard_state = {
            'base_idx': hard_idx,
            'perm': np.random.permutation(hard_idx),
            'cursor': 0,
        }

        rec_mix_counts = {'cross': 0, 'within': 0, 'hard': 0, 'fallback': 0}

        for _ in range(n_pos):
            bucket_counts = split_counts(neg_ratio, mix)
            neg_rho_chunks = []
            neg_theta_chunks = []
            neg_input_chunks = []
            neg_mask_chunks = []
            neg_names = []

            n_cross = bucket_counts[0]
            cross_pick = draw_without_replacement_cycle(cross_valid_state, n_cross)
            if len(cross_pick) > 0:
                rows = get_cross_rows_and_names(pool, cross_pick, records, pair_cache)
                neg_rho_chunks.extend(rows[0])
                neg_theta_chunks.extend(rows[1])
                neg_input_chunks.extend(rows[2])
                neg_mask_chunks.extend(rows[3])
                neg_names.extend(rows[4])
                rec_mix_counts['cross'] += len(cross_pick)

            n_within = bucket_counts[1]
            within_pick = draw_without_replacement_cycle(within_state, n_within)
            if len(within_pick) > 0:
                abs_idx = k_neg2[within_pick]
                neg_rho_chunks.append(rec_arrays['p2_rho'][abs_idx])
                neg_theta_chunks.append(rec_arrays['p2_theta'][abs_idx])
                neg_input_chunks.append(rec_arrays['p2_input'][abs_idx])
                neg_mask_chunks.append(rec_arrays['p2_mask'][abs_idx])
                neg_names.extend(['{}_p2_{}'.format(rec['ppi_pair_id'], int(ii)) for ii in abs_idx])
                rec_mix_counts['within'] += len(within_pick)

            n_hard = bucket_counts[2]
            hard_pick = draw_without_replacement_cycle(hard_state, n_hard)
            if len(hard_pick) > 0:
                abs_idx = k_neg2[hard_pick]
                neg_rho_chunks.append(rec_arrays['p2_rho'][abs_idx])
                neg_theta_chunks.append(rec_arrays['p2_theta'][abs_idx])
                neg_input_chunks.append(rec_arrays['p2_input'][abs_idx])
                neg_mask_chunks.append(rec_arrays['p2_mask'][abs_idx])
                neg_names.extend(['{}_p2_{}'.format(rec['ppi_pair_id'], int(ii)) for ii in abs_idx])
                rec_mix_counts['hard'] += len(hard_pick)

            if len(neg_names) < neg_ratio:
                n_missing = neg_ratio - len(neg_names)
                fallback_source = 'cross_valid'
                fallback_pick = draw_without_replacement_cycle(cross_valid_state, n_missing)
                if len(fallback_pick) == 0:
                    fallback_source = 'within'
                    fallback_pick = draw_without_replacement_cycle(within_state, n_missing)
                if len(fallback_pick) == 0:
                    fallback_source = 'cross_all'
                    fallback_pick = draw_without_replacement_cycle(cross_all_state, n_missing)
                if len(fallback_pick) == 0:
                    raise RuntimeError('Could not sample fallback negatives for {}'.format(rec['ppi_pair_id']))

                if fallback_source in ['cross_valid', 'cross_all']:
                    rows = get_cross_rows_and_names(pool, fallback_pick, records, pair_cache)
                    neg_rho_chunks.extend(rows[0])
                    neg_theta_chunks.extend(rows[1])
                    neg_input_chunks.extend(rows[2])
                    neg_mask_chunks.extend(rows[3])
                    neg_names.extend(rows[4])
                else:
                    abs_idx = k_neg2[fallback_pick]
                    neg_rho_chunks.append(rec_arrays['p2_rho'][abs_idx])
                    neg_theta_chunks.append(rec_arrays['p2_theta'][abs_idx])
                    neg_input_chunks.append(rec_arrays['p2_input'][abs_idx])
                    neg_mask_chunks.append(rec_arrays['p2_mask'][abs_idx])
                    neg_names.extend(['{}_p2_{}'.format(rec['ppi_pair_id'], int(ii)) for ii in abs_idx])
                rec_mix_counts['fallback'] += len(fallback_pick)

            writer.add_neg_block(
                split=rec['split'],
                neg_rows=(
                    concat_or_empty(neg_rho_chunks, 2),
                    concat_or_empty(neg_theta_chunks, 2),
                    concat_or_empty(neg_input_chunks, 3),
                    concat_or_empty(neg_mask_chunks, 2),
                ),
                neg_names=neg_names,
            )

        writer.accumulate_mix_counts(rec_mix_counts)
        if writer.should_flush():
            writer.flush()

    subset_manifest = writer.finalize(
        {
            'subset_id': args.subset_id,
            'num_subsets': args.num_subsets,
            'catalog_fingerprint': catalog_manifest.get('fingerprint'),
            'assigned_records': len(assigned_record_indices),
        }
    )
    print('Subset complete: pos={} neg={} parts={}'.format(subset_manifest['pos_count'], subset_manifest['neg_count'], subset_manifest['parts_count']))


if __name__ == '__main__':
    main()
