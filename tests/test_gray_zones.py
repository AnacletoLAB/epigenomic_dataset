from epigenomic_dataset import (
    active_enhancers_vs_active_promoters,
    active_enhancers_vs_inactive_enhancers,
    active_promoters_vs_inactive_promoters
)

def test_gray_zones():
    """Test that preprocessed epigenomic data are loadable."""
    X, y = active_enhancers_vs_inactive_enhancers()
    assert X.shape[0] > 0
    assert X.shape[0] == y.shape[0]
    assert all(X.index == y.index)
    assert set(y.values.flatten()) != set([0.0, 1.0])
    X, y = active_enhancers_vs_inactive_enhancers(binarize=True)
    assert X.shape[0] > 0
    assert X.shape[0] == y.shape[0]
    assert all(X.index == y.index)
    assert set(y.values.flatten()) == set([0.0, 1.0])

    X, y = active_promoters_vs_inactive_promoters()
    assert X.shape[0] > 0
    assert X.shape[0] == y.shape[0]
    assert all(X.index == y.index)
    assert set(y.values.flatten()) != set([0.0, 1.0])
    X, y = active_promoters_vs_inactive_promoters(binarize=True)
    assert X.shape[0] > 0
    assert X.shape[0] == y.shape[0]
    assert all(X.index == y.index)
    assert set(y.values.flatten()) == set([0.0, 1.0])