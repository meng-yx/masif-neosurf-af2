custom_params = {}
custom_params['cache_dir'] = 'nn_models/sc05/cache/'
custom_params['model_dir'] = 'nn_models/sc05/all_feat/model_data/'
custom_params['desc_dir'] = 'descriptors/sc05/all_feat/'
custom_params['feat_mask'] = [1.0, 1.0, 1.0, 1.0, 1.0]
# Negative sampling defaults (baseline-compatible).
custom_params['neg_ratio'] = 1
custom_params['neg_mix_cross_complex'] = 0.0
custom_params['neg_mix_within_complex'] = 1.0
custom_params['neg_mix_hard'] = 0.0
custom_params['enforce_diff_pdb_for_cross'] = True
custom_params['hard_negative_topk'] = 200
custom_params['neg_loss_weight'] = 1.0
