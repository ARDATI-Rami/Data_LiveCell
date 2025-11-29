# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025  Rami Ardati

"""Topology module for cell shape analysis and metrics."""

from .metrics import (
    infer_support,
    tally_roi_counts,
    dirichlet_posterior_mean,
    js_divergence,
    wasserstein1_discrete,
    multinomial_residuals,
    topological_charge,
    compute_topology_deviation_metrics
)

from .distributions import (
    SHAPE_NAMES,
    get_shape_name,
    compute_global_distribution,
    compute_frame_distribution,
    mean_sidedness_per_roi,
    compute_roi_deviations
)

__all__ = [
    # Metrics
    'infer_support',
    'tally_roi_counts',
    'dirichlet_posterior_mean',
    'js_divergence',
    'wasserstein1_discrete',
    'multinomial_residuals',
    'topological_charge',
    'compute_topology_deviation_metrics',
    # Distributions
    'SHAPE_NAMES',
    'get_shape_name',
    'compute_global_distribution',
    'compute_frame_distribution',
    'mean_sidedness_per_roi',
    'compute_roi_deviations',
]

