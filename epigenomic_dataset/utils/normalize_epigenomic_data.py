"""Module providing normalization for epigenomic data."""
from typing import Tuple
import numpy as np
from sklearn.impute import KNNImputer
from sklearn.preprocessing import StandardScaler, MinMaxScaler


def normalize_epigenomic_data(
    train_x: np.ndarray,
    test_x: np.ndarray = None
) -> Tuple[np.ndarray]:
    """Return imputed and normalized epigenomic data.

    We fit the imputation and normalization on the training data and
    apply it to both the training data and the test data.

    Parameters
    -------------------------
    train_x: np.ndarray,
        Training data to use to fit the imputer and scaled.
    test_x: np.ndarray = None,
        Test data to be normalized.

    Returns
    -------------------------
    Tuple with imputed and scaled train and test data.
    """
    # Create the imputer and scaler object
    imputer = KNNImputer()
    scaler = MinMaxScaler()
    # Fit the imputer object
    imputer.fit(train_x)
    # Impute the train and test data
    imputed_train_x = imputer.transform(train_x)
    if test_x is not None:
        imputed_test_x = imputer.transform(test_x)
    # Fit the scaler object
    scaler.fit(imputed_train_x)
    # Scale the train and test data
    scaled_train_x = scaler.transform(imputed_train_x)
    if test_x is not None:
        scaled_test_x = scaler.transform(imputed_test_x)
    if test_x is not None:
        # Return the normalized data
        return scaled_train_x, scaled_test_x
    return scaled_train_x