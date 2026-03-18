import time
import math
from sklearn import metrics
import numpy as np
import sys
import os
from IPython.core.debugger import set_trace
from sklearn.metrics import accuracy_score, roc_auc_score

# Features and theta are flipped for the binder in construct_batch (except for hydrophobicity).
def construct_batch(
    binder_rho_wrt_center,
    binder_theta_wrt_center,
    binder_input_feat,
    binder_mask,
    c_pos_training_idx,
    pos_rho_wrt_center,
    pos_theta_wrt_center,
    pos_input_feat,
    pos_mask,
    c_neg_training_idx,
    neg_rho_wrt_center,
    neg_theta_wrt_center,
    neg_input_feat,
    neg_mask,
    neg_ratio=1,
    c_neg_training_idx_2=None,
):
    neg_ratio = max(1, int(neg_ratio))
    c_pos_training_idx = np.asarray(c_pos_training_idx, dtype=int)
    c_neg_training_idx = np.asarray(c_neg_training_idx, dtype=int)
    # Keep the 4 equal-sized chunks required by the loss.
    # Repeat positive/binder anchors to match the number of negatives.
    if len(c_neg_training_idx) == 0:
        raise ValueError("c_neg_training_idx is empty.")
    if len(c_pos_training_idx) == 0:
        raise ValueError("c_pos_training_idx is empty.")
    if len(c_neg_training_idx) != len(c_pos_training_idx):
        rep = int(np.ceil(float(len(c_neg_training_idx)) / float(len(c_pos_training_idx))))
        c_pos_training_idx = np.repeat(c_pos_training_idx, rep)[: len(c_neg_training_idx)]

    batch_rho_coords_binder = np.expand_dims(
        binder_rho_wrt_center[c_pos_training_idx], 2
    )
    batch_theta_coords_binder = np.expand_dims(
        binder_theta_wrt_center[c_pos_training_idx], 2
    )
    batch_input_feat_binder = binder_input_feat[c_pos_training_idx]
    batch_mask_binder = binder_mask[c_pos_training_idx]

    batch_rho_coords_pos = np.expand_dims(pos_rho_wrt_center[c_pos_training_idx], 2)
    batch_theta_coords_pos = np.expand_dims(pos_theta_wrt_center[c_pos_training_idx], 2)
    batch_input_feat_pos = pos_input_feat[c_pos_training_idx]
    batch_mask_pos = pos_mask[c_pos_training_idx]

    # Negate the input_features of the binder, except the last column.
    batch_input_feat_binder = -batch_input_feat_binder
    # TODO: This should not be like this ... it is a hack.
    if batch_input_feat_binder.shape[2] == 5 or batch_input_feat_binder.shape[2] == 3:
        batch_input_feat_binder[:, :, -1] = -batch_input_feat_binder[
            :, :, -1
        ]  # Do not negate hydrophobicity.
    # Also negate the theta coords for the binder.
    batch_theta_coords_binder = 2 * np.pi - batch_theta_coords_binder

    batch_rho_coords_neg = np.expand_dims(neg_rho_wrt_center[c_neg_training_idx], 2)
    batch_theta_coords_neg = np.expand_dims(neg_theta_wrt_center[c_neg_training_idx], 2)
    batch_input_feat_neg = neg_input_feat[c_neg_training_idx]
    batch_mask_neg = neg_mask[c_neg_training_idx]

    batch_rho_coords_neg_2 = batch_rho_coords_binder.copy()
    batch_theta_coords_neg_2 = batch_theta_coords_binder.copy()
    batch_input_feat_neg_2 = batch_input_feat_binder.copy()
    batch_mask_neg_2 = batch_mask_binder.copy()

    batch_rho_coords = np.concatenate(
        [
            batch_rho_coords_pos,
            batch_rho_coords_binder,
            batch_rho_coords_neg,
            batch_rho_coords_neg_2,
        ],
        axis=0,
    )
    batch_theta_coords = np.concatenate(
        [
            batch_theta_coords_pos,
            batch_theta_coords_binder,
            batch_theta_coords_neg,
            batch_theta_coords_neg_2,
        ],
        axis=0,
    )
    batch_input_feat = np.concatenate(
        [
            batch_input_feat_pos,
            batch_input_feat_binder,
            batch_input_feat_neg,
            batch_input_feat_neg_2,
        ],
        axis=0,
    )
    batch_mask = np.concatenate(
        [batch_mask_pos, batch_mask_binder, batch_mask_neg, batch_mask_neg_2], axis=0
    )
    # expand the last dimension of the mask (batch_size, max_points_patch, 1)
    batch_mask = np.expand_dims(batch_mask, 2)

    return batch_rho_coords, batch_theta_coords, batch_input_feat, batch_mask


def compute_dists(descs1, descs2):
    dists = np.sqrt(np.sum(np.square(descs1 - descs2), axis=1))
    return dists


def construct_batch_val_test(
    c_idx, rho_wrt_center, theta_wrt_center, input_feat, mask, flip=False
):
    batch_rho_coords = np.expand_dims(rho_wrt_center[c_idx], 2)
    batch_theta_coords = np.expand_dims(theta_wrt_center[c_idx], 2)
    batch_input_feat = input_feat[c_idx]
    batch_mask = mask[c_idx]
    batch_mask = np.expand_dims(batch_mask,2)
    # Flip features and theta (except hydrophobicity)
    if flip:
        batch_input_feat = -batch_input_feat
        batch_theta_coords = 2 * np.pi - batch_theta_coords
        assert len(batch_input_feat.shape) == 3
        # Hydrophobicity is not flipped. -- FIx this.
        if batch_input_feat.shape[2] == 5 or batch_input_feat.shape[2] == 3:
            batch_input_feat[:, :, -1] = -batch_input_feat[:, :, -1]

    return batch_rho_coords, batch_theta_coords, batch_input_feat, batch_mask


def compute_val_test_desc(
    learning_obj,
    idx,
    rho_wrt_center,
    theta_wrt_center,
    input_feat,
    mask,
    batch_size=100,
    flip=False,
):
    all_descs = []
    num_batches = int(np.ceil(float(len(idx)) / float(batch_size)))
    # Compute all desc for positive shapes.
    for kk in range(num_batches):
        c_idx = idx[np.arange(kk * batch_size, min((kk + 1) * batch_size, len(idx)))]

        batch_rho_coords, batch_theta_coords, batch_input_feat, batch_mask = construct_batch_val_test(
            c_idx, rho_wrt_center, theta_wrt_center, input_feat, mask, flip=flip
        )

        feed_dict = {
            learning_obj.rho_coords: batch_rho_coords,
            learning_obj.theta_coords: batch_theta_coords,
            learning_obj.input_feat: batch_input_feat,
            learning_obj.mask: batch_mask,
            learning_obj.keep_prob: 1.0,
        }

        desc = learning_obj.session.run([learning_obj.global_desc], feed_dict=feed_dict)
        desc = np.squeeze(desc)
        if len(desc.shape) == 1:
            desc = np.expand_dims(desc, 0)
        all_descs.append(desc)

    if len(all_descs) > 1:
        all_descs = np.concatenate(all_descs, axis=0)
    else:
        all_descs = all_descs[0]
    return all_descs


def compute_roc_auc(pos, neg):
    labels = np.concatenate([np.ones((len(pos))), np.zeros((len(neg)))])
    dist_pairs = np.concatenate([pos, neg])
    return metrics.roc_auc_score(labels, dist_pairs)


def sample_equal_length(descs1, descs2):
    """Return two descriptor arrays with identical length via random subsampling."""
    n1 = len(descs1)
    n2 = len(descs2)
    if n1 == 0 or n2 == 0:
        return descs1[:0], descs2[:0]
    n = min(n1, n2)
    idx1 = np.random.choice(n1, size=n, replace=False)
    idx2 = np.random.choice(n2, size=n, replace=False)
    return descs1[idx1], descs2[idx2]


class Logger:
    def __init__(self, logfile):
        self.logfile = logfile

    def write(self, text):
        with open(self.logfile, 'a') as f:
            f.write(text)
        print(text, flush=True)


# Randomly pick
def train_ppi_search(
    learning_obj,
    params,
    binder_rho_wrt_center,
    binder_theta_wrt_center,
    binder_input_feat,
    binder_mask,
    pos_training_idx,
    pos_val_idx,
    pos_test_idx,
    pos_rho_wrt_center,
    pos_theta_wrt_center,
    pos_input_feat,
    pos_mask,
    neg_training_idx,
    neg_val_idx,
    neg_test_idx,
    neg_rho_wrt_center,
    neg_theta_wrt_center,
    neg_input_feat,
    neg_mask,
    num_iterations=1000000,
    num_iter_test=1000,
    batch_size=32,
    batch_size_val_test=1000,
):

    out_dir = params["model_dir"]
    logfile = Logger(out_dir + "log.txt")
    neg_ratio = max(1, int(params.get("neg_ratio", 1)))
    logfile.write("Configured neg_ratio: {}\n".format(neg_ratio))
    logfile.write(
        "Number of training positive shapes: {}\n".format(len(pos_training_idx))
    )
    logfile.write("Number of validation positive shapes: {}\n".format(len(pos_val_idx)))
    logfile.write("Number of testing positive shapes: {}\n".format(len(pos_test_idx)))

    logfile.write(
        "Number of training negative shapes: {}\n".format(len(neg_training_idx))
    )
    logfile.write("Number of validation negative shapes: {}\n".format(len(neg_val_idx)))
    logfile.write("Number of testing negative shapes: {}\n".format(len(neg_test_idx)))

    list_training_loss = []
    list_training_norm_grad = []
    iter_time = []
    best_val_auc = 0

    pos_training_idx_copy = np.copy(pos_training_idx)
    neg_training_idx_copy = np.copy(neg_training_idx)

    logfile.write("Number of iterations: {}\n".format(num_iterations))

    iter_training_loss = []
    iter_pos_score = []
    iter_neg_score = []

    test_iterations_log = [
        0,
        1,
        2,
        4,
        8,
        16,
        32,
        64,
        128,
        256,
        512,
        1024,
        2048,
        4096,
        8192,
        16384,
        32178,
    ]

    train_start_time = time.time()
    for num_iter in range(num_iterations):
        # Read dataset for training.
        tic = time.time()

        np.random.shuffle(pos_training_idx_copy)
        np.random.shuffle(neg_training_idx_copy)

        # Keep the effective feed size close to batch_size while allowing neg_ratio > 1.
        n_anchor = max(1, batch_size // (4 * neg_ratio))
        n_pairs = max(1, n_anchor * neg_ratio)
        c_pos_training_idx = pos_training_idx_copy[:n_anchor]
        c_neg_training_idx = neg_training_idx_copy[:n_pairs]
        if len(c_pos_training_idx) == 0:
            c_pos_training_idx = np.random.choice(pos_training_idx_copy, size=1, replace=True)
        if len(c_neg_training_idx) == 0:
            c_neg_training_idx = np.random.choice(neg_training_idx_copy, size=1, replace=True)

        c_neg_training_idx_2 = None

        # Features and theta are flipped for the binder in construct_batch (except for hydrophobicity).
        batch_rho_coords, batch_theta_coords, batch_input_feat, batch_mask = construct_batch(
            binder_rho_wrt_center,
            binder_theta_wrt_center,
            binder_input_feat,
            binder_mask,
            c_pos_training_idx,
            pos_rho_wrt_center,
            pos_theta_wrt_center,
            pos_input_feat,
            pos_mask,
            c_neg_training_idx,
            neg_rho_wrt_center,
            neg_theta_wrt_center,
            neg_input_feat,
            neg_mask,
            neg_ratio=neg_ratio,
            c_neg_training_idx_2=c_neg_training_idx_2,
        )

        feed_dict = {
            learning_obj.rho_coords: batch_rho_coords,
            learning_obj.theta_coords: batch_theta_coords,
            learning_obj.input_feat: batch_input_feat,
            learning_obj.mask: batch_mask,
            learning_obj.keep_prob: 0.5,
        }

        # Do not train during the first iteration
        if num_iter == 0:
            [score] = learning_obj.session.run(
                [learning_obj.score], feed_dict=feed_dict
            )
            training_loss = 0
        else:
            _, training_loss, norm_grad, score = learning_obj.session.run(
                [
                    learning_obj.optimizer,
                    learning_obj.data_loss,
                    learning_obj.norm_grad,
                    learning_obj.score,
                ],
                feed_dict=feed_dict,
            )

        n = len(score) // 2
        try:
            pos_score = score[:n]
            neg_score = score[n:]
        except:
            print(score)
            sys.exit(1)
        iter_pos_score = np.concatenate([pos_score, iter_pos_score], axis=0)
        iter_neg_score = np.concatenate([neg_score, iter_neg_score], axis=0)
        list_training_loss.append(training_loss)

        if num_iter % num_iter_test == 0:
            logfile.write("{} training iterations finished after {:.2f} seconds.\n".format(num_iter, time.time() - train_start_time))

            logfile.write("Validating and testing.\n ")
            roc_auc = 1 - compute_roc_auc(iter_pos_score, iter_neg_score)
            logfile.write("training_loss: {}\n".format(np.mean(list_training_loss)))
            # print("Approx Training ROC AUC: {}\n ".format(roc_auc))
            logfile.write(
                "Mean training positive score: {} \n".format(
                    np.mean(1.0 / iter_pos_score)
                )
            )
            logfile.write(
                "Mean training negative score: {} \n".format(
                    np.mean(1.0 / iter_neg_score)
                )
            )
            logfile.write("Training ROC AUC: {}\n ".format(roc_auc))

            iter_pos_score = []
            iter_neg_score = []
            list_training_loss = []

            training_time_entry = time.time() - tic
            logfile.write(
                "Training {} entries took {:.2f}s \n".format(
                    len(batch_mask) * num_iter_test, training_time_entry
                )
            )

            tic = time.time()
            # Compute validation descriptors.
            pos_desc = compute_val_test_desc(
                learning_obj,
                pos_val_idx,
                pos_rho_wrt_center,
                pos_theta_wrt_center,
                pos_input_feat,
                pos_mask,
                batch_size=batch_size_val_test,
            )

            binder_desc = compute_val_test_desc(
                learning_obj,
                pos_val_idx,
                binder_rho_wrt_center,
                binder_theta_wrt_center,
                binder_input_feat,
                binder_mask,
                batch_size=batch_size_val_test,
                flip=True,
            )

            neg_desc = compute_val_test_desc(
                learning_obj,
                neg_val_idx,
                neg_rho_wrt_center,
                neg_theta_wrt_center,
                neg_input_feat,
                neg_mask,
                batch_size=batch_size_val_test,
            )

            neg_desc_2 = binder_desc.copy()
            neg_desc, neg_desc_2 = sample_equal_length(neg_desc, neg_desc_2)

            # Compute val ROC AUC.
            pos_dists = compute_dists(pos_desc, binder_desc)
            neg_dists = compute_dists(neg_desc, neg_desc_2)
            try:
                val_auc = 1 - compute_roc_auc(pos_dists, neg_dists)
            except:
                print(np.min(pos_dists))
                print(np.min(neg_dists))
                sys.exit(1)
            val_time = time.time() - tic

            logfile.write(
                "Iteration {} validation roc auc: {}\n".format(num_iter, val_auc)
            )

            logfile.write(
                "Mean validation positive score: {} \n".format(np.mean(pos_dists))
            )
            logfile.write(
                "Mean validation negative score: {} \n".format(np.mean(neg_dists))
            )
            # Compute TEST ROC AUC.
            tic = time.time()
            pos_desc = compute_val_test_desc(
                learning_obj,
                pos_test_idx,
                pos_rho_wrt_center,
                pos_theta_wrt_center,
                pos_input_feat,
                pos_mask,
                batch_size=batch_size_val_test,
            )

            binder_desc = compute_val_test_desc(
                learning_obj,
                pos_test_idx,
                binder_rho_wrt_center,
                binder_theta_wrt_center,
                binder_input_feat,
                binder_mask,
                batch_size=batch_size_val_test,
                flip=True,
            )

            neg_desc = compute_val_test_desc(
                learning_obj,
                neg_test_idx,
                neg_rho_wrt_center,
                neg_theta_wrt_center,
                neg_input_feat,
                neg_mask,
                batch_size=batch_size_val_test,
            )

            neg_desc_2 = binder_desc.copy()

            # Compute test ROC AUC.
            pos_dists = compute_dists(pos_desc, binder_desc)
            neg_desc, neg_desc_2 = sample_equal_length(neg_desc, neg_desc_2)
            neg_dists = compute_dists(neg_desc, neg_desc_2)
            test_auc = 1 - compute_roc_auc(pos_dists, neg_dists)

            test_time = time.time() - tic
            logfile.write("Iteration {} test roc auc: {}\n".format(num_iter, test_auc))
            logfile.write(
                "Iteration time: {} validation time: {} test time: {}\n".format(
                    np.mean(iter_time), val_time, test_time
                )
            )
            logfile.write("Mean test positive score: {} \n".format(np.mean(pos_dists)))
            logfile.write("Mean test negative score: {} \n".format(np.mean(neg_dists)))

            tic = time.time()

            if val_auc > best_val_auc:
                logfile.write(">>> Saving model.\n")
                best_val_auc = val_auc
                output_model = out_dir + "model"
                learning_obj.saver.save(learning_obj.session, output_model)
                np.save(out_dir + "pos_dists.npy", pos_dists)
                np.save(out_dir + "pos_test_idx.npy", neg_test_idx)
                np.save(out_dir + "neg_dists.npy", neg_dists)
                np.save(out_dir + "pos_desc.npy", pos_desc)
                np.save(out_dir + "binder_desc.npy", binder_desc)
                np.save(out_dir + "neg_desc.npy", neg_desc)
                np.save(out_dir + "neg_test_idx.npy", neg_test_idx)
                np.save(out_dir + "neg_desc_2.npy", neg_desc_2)

