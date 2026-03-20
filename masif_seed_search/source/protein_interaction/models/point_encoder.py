import tensorflow as tf
from tensorflow.keras import layers


class PointEncoder(layers.Layer):
    """Shared point-wise encoder (1x1 MLP over vertices)."""

    def __init__(self, hidden_dim: int = 128, dropout: float = 0.2, name: str = "point_encoder"):
        super().__init__(name=name)
        self.d1 = layers.Dense(hidden_dim, activation="relu")
        self.d2 = layers.Dense(hidden_dim, activation="relu")
        self.d3 = layers.Dense(hidden_dim, activation="relu")
        self.dropout = layers.Dropout(dropout)

    def call(self, x, training=False):
        x = self.d1(x)
        x = self.dropout(x, training=training)
        x = self.d2(x)
        x = self.dropout(x, training=training)
        x = self.d3(x)
        return x

