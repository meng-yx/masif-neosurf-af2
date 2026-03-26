custom_params = {}
custom_params["cache_dir"] = "nn_models/sc05/cache/model_data_cross_3to1/"
custom_params["model_dir"] = "nn_models/sc05/all_feat/model_data_cross_3to1/"
custom_params["desc_dir"] = "descriptors/sc05/all_feat/"
custom_params["feat_mask"] = [1.0, 1.0, 1.0, 1.0, 1.0]

# Baseline behavior (legacy-compatible).
custom_params["neg_ratio"] = 3
custom_params["neg_mix_cross_complex"] = 1.0        # pair of patches coming from different complexes
custom_params["neg_mix_within_complex"] = 0.0       # pair of patches coming from the same complex
custom_params["neg_mix_hard"] = 0.0                 # pair of patches from within the same complex but are closest to the interface cutoff (just outside the interface region) - maybe not a good idea
custom_params["enforce_diff_pdb_for_cross"] = True
custom_params["hard_negative_topk"] = 200
custom_params["neg_loss_weight"] = 1.0
