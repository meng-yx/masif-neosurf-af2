import math

import tensorflow as tf
from tensorflow.keras import layers

from protein_interaction.models.point_encoder import PointEncoder


def masked_mean(x, mask):
    mask = tf.cast(mask, tf.float32)
    denom = tf.reduce_sum(mask, axis=1, keepdims=True) + 1e-6
    weighted = x * tf.expand_dims(mask, axis=-1)
    return tf.reduce_sum(weighted, axis=1) / denom


def masked_max(x, mask):
    mask = tf.cast(mask, tf.float32)
    minus_inf = tf.fill(tf.shape(x), tf.cast(-1e9, x.dtype))
    cond = tf.tile(tf.expand_dims(mask > 0.0, axis=-1), [1, 1, tf.shape(x)[-1]])
    x_masked = tf.where(cond, x, minus_inf)
    return tf.reduce_max(x_masked, axis=1)


class CrossSetInteractionV1(tf.keras.Model):
    """Permutation-invariant protein-pair classifier with masked cross-attention."""

    def __init__(self, hidden_dim=128, classifier_hidden_dim=128, dropout=0.2, name="cross_set_v1"):
        super().__init__(name=name)
        self.encoder = PointEncoder(hidden_dim=hidden_dim, dropout=dropout)
        self.classifier = tf.keras.Sequential(
            [
                layers.Dense(classifier_hidden_dim, activation="relu"),
                layers.Dropout(dropout),
                layers.Dense(classifier_hidden_dim // 2, activation="relu"),
                layers.Dropout(dropout),
                layers.Dense(1, activation="sigmoid"),
            ]
        )

    def _cross_attention(self, query_h, query_mask, key_h, key_mask):
        dk = tf.cast(tf.shape(query_h)[-1], tf.float32)
        logits = tf.matmul(query_h, key_h, transpose_b=True) / tf.sqrt(dk + 1e-6)

        attn_mask = tf.expand_dims(query_mask, 2) * tf.expand_dims(key_mask, 1)
        minus_inf = tf.fill(tf.shape(logits), tf.cast(-1e9, logits.dtype))
        logits = tf.where(attn_mask > 0.0, logits, minus_inf)

        weights = tf.nn.softmax(logits, axis=-1)
        context = tf.matmul(weights, key_h)
        return context

    def call(self, inputs, training=False):
        q_desc = inputs["query_desc"]
        q_xyz = inputs["query_xyz"]
        q_mask = inputs["query_mask"]
        m_desc = inputs["matched_desc"]
        m_xyz = inputs["matched_xyz"]
        m_mask = inputs["matched_mask"]

        q_in = tf.concat([q_desc, q_xyz], axis=-1)
        m_in = tf.concat([m_desc, m_xyz], axis=-1)

        q_h = self.encoder(q_in, training=training)
        m_h = self.encoder(m_in, training=training)

        q_ctx = self._cross_attention(q_h, q_mask, m_h, m_mask)
        m_ctx = self._cross_attention(m_h, m_mask, q_h, q_mask)

        q_joint = tf.concat([q_h, q_ctx, tf.abs(q_h - q_ctx)], axis=-1)
        m_joint = tf.concat([m_h, m_ctx, tf.abs(m_h - m_ctx)], axis=-1)

        q_mean = masked_mean(q_joint, q_mask)
        q_max = masked_max(q_joint, q_mask)
        m_mean = masked_mean(m_joint, m_mask)
        m_max = masked_max(m_joint, m_mask)

        global_repr = tf.concat([q_mean, q_max, m_mean, m_max], axis=-1)
        out = self.classifier(global_repr, training=training)
        return out

