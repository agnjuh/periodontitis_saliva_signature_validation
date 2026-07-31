# Central path registry. Keep all filesystem conventions here.


# Shared project paths

CONFIG_DIR = Path("config")
DATA_DIR = Path("data")
RESULTS_DIR = Path("results")
LOG_DIR = Path("logs")
SCRIPT_DIR = Path("scripts")

SCREENING_DIR = RESULTS_DIR / "screening"
INTEGRATION_DIR = RESULTS_DIR / "integration"

DATASET_REGISTRY = CONFIG_DIR / "dataset_registry.csv"
MODULE_DEFINITION = CONFIG_DIR / "biological_modules.csv"
DATASET_PRIORITISATION = SCREENING_DIR / "dataset_prioritisation.csv"
DATASET_SCREENING_SUMMARY = SCREENING_DIR / "dataset_screening_summary.txt"


# PXD043901 proteomics

PXD_DATASET = "PXD043901"

PXD_RAW_DIR = DATA_DIR / "raw" / PXD_DATASET
PXD_PROCESSED_DIR = DATA_DIR / "processed" / PXD_DATASET

PXD_RESULTS_DIR = RESULTS_DIR / PXD_DATASET
PXD_QC_DIR = PXD_RESULTS_DIR / "qc"
PXD_QC_FIGURE_DIR = PXD_RESULTS_DIR / "figures" / "qc"
PXD_TRAJECTORY_DIR = PXD_RESULTS_DIR / "trajectories"
PXD_MODULE_DIR = PXD_RESULTS_DIR / "modules"
PXD_SENSITIVITY_DIR = PXD_MODULE_DIR / "sensitivity"
PXD_CLINICAL_DIR = PXD_RESULTS_DIR / "clinical_associations"
PXD_MARKER_DIR = PXD_RESULTS_DIR / "workflow_markers"

PXD_LOG_DIR = LOG_DIR / PXD_DATASET
PXD_SCRIPT_DIR = SCRIPT_DIR / "pxd043901"

RAW_DOWNLOAD_MARKER = PXD_MARKER_DIR / f"{PXD_DATASET}_download_complete.done"

METADATA = PXD_PROCESSED_DIR / f"{PXD_DATASET}_sample_metadata.csv"
METADATA_CLEAN = PXD_PROCESSED_DIR / "sample_metadata_clean.csv"
LFQ_RAW = PXD_PROCESSED_DIR / f"{PXD_DATASET}_lfq_matrix_raw.csv"
LFQ_LOG2 = PXD_PROCESSED_DIR / f"{PXD_DATASET}_lfq_matrix_log2.csv"
PROTEIN_ANNOTATION = PXD_PROCESSED_DIR / f"{PXD_DATASET}_protein_annotation.csv"
PROTEIN_ANNOTATION_EXTENDED = (
    PXD_PROCESSED_DIR / f"{PXD_DATASET}_protein_annotation_extended.csv"
)
MATRIX_SUMMARY = PXD_PROCESSED_DIR / f"{PXD_DATASET}_matrix_summary.txt"
ORGANISM_SUMMARY = PXD_PROCESSED_DIR / f"{PXD_DATASET}_organism_summary.txt"
DELTA_T1 = PXD_PROCESSED_DIR / f"{PXD_DATASET}_delta_T1_vs_T0.csv"
DELTA_T6 = PXD_PROCESSED_DIR / f"{PXD_DATASET}_delta_T6_vs_T0.csv"

QC_SAMPLE_METRICS = PXD_QC_DIR / "sample_qc_metrics.csv"
QC_PROTEIN_METRICS = PXD_QC_DIR / "protein_qc_metrics.csv"
QC_CORRELATION = PXD_QC_DIR / "sample_spearman_correlation.csv"
QC_PCA_COORDINATES = PXD_QC_DIR / "sample_pca_coordinates.csv"
QC_PCA_VARIANCE = PXD_QC_DIR / "pca_explained_variance.csv"
QC_SUMMARY = PXD_QC_DIR / "qc_summary.txt"

PROTEIN_RESPONSE_SUMMARY = PXD_TRAJECTORY_DIR / "protein_response_summary.csv"
ROBUST_PRELIMINARY = (
    PXD_TRAJECTORY_DIR / "robust_human_responses_preliminary.csv"
)
RESPONSE_SUMMARY = PXD_TRAJECTORY_DIR / "response_summary.txt"
LONGITUDINAL_TESTS = PXD_TRAJECTORY_DIR / "longitudinal_statistical_tests.csv"
STRINGENT_RESPONSES = PXD_TRAJECTORY_DIR / "stringent_human_responses.csv"
RANKED_RESPONSES = PXD_TRAJECTORY_DIR / "ranked_human_responses.csv"
STATISTICAL_TEST_SUMMARY = PXD_TRAJECTORY_DIR / "statistical_test_summary.txt"
TRAJECTORY_CLASSES = PXD_TRAJECTORY_DIR / "human_protein_trajectory_classes.csv"
PERSISTENT_RESPONSES = PXD_TRAJECTORY_DIR / "persistent_human_responses.csv"
TRAJECTORY_COUNTS = PXD_TRAJECTORY_DIR / "trajectory_class_counts.csv"
TRAJECTORY_SUMMARY = PXD_TRAJECTORY_DIR / "trajectory_summary.txt"

MODULE_SCORES = PXD_MODULE_DIR / "sample_module_scores.csv"
MODULE_COVERAGE = PXD_MODULE_DIR / "module_coverage.csv"
MODULE_PATIENT_DELTAS = PXD_MODULE_DIR / "patient_module_response_deltas.csv"
MODULE_STATISTICS = PXD_MODULE_DIR / "module_longitudinal_statistics.csv"
MODULE_SUMMARY = PXD_MODULE_DIR / "module_analysis_summary.txt"
MODULE_LOO_RESULTS = PXD_SENSITIVITY_DIR / "leave_one_protein_out_results.csv"
MODULE_LOO_SUMMARY = PXD_SENSITIVITY_DIR / "leave_one_out_summary.csv"
MODULE_DELTA_CORRELATIONS = PXD_SENSITIVITY_DIR / "module_delta_correlations.csv"
MODULE_HIGH_CORRELATIONS = PXD_SENSITIVITY_DIR / "high_module_correlations.csv"
MODULE_SENSITIVITY_SUMMARY = (
    PXD_SENSITIVITY_DIR / "module_sensitivity_summary.txt"
)


# GSE312585 salivary extracellular-vesicle miRNA

GSE_DATASET = "GSE312585"

GSE_RAW_DIR = DATA_DIR / "raw" / GSE_DATASET
GSE_PROCESSED_DIR = DATA_DIR / "processed" / GSE_DATASET

GSE_RESULTS_DIR = RESULTS_DIR / GSE_DATASET
GSE_TABLE_DIR = GSE_RESULTS_DIR / "tables"
GSE_FIGURE_DIR = GSE_RESULTS_DIR / "figures"
GSE_TEXT_DIR = GSE_RESULTS_DIR / "text"

GSE_LOG_DIR = LOG_DIR / GSE_DATASET
GSE_SCRIPT_DIR = SCRIPT_DIR / "gse312585"

GSE_METADATA = GSE_PROCESSED_DIR / "sample_metadata.tsv"
GSE_RAW_COUNT_MATRIX = GSE_PROCESSED_DIR / "raw_count_matrix.tsv"
GSE_MIRNA_RAW_COUNTS = GSE_PROCESSED_DIR / "mirna_raw_counts.tsv"
GSE_MIRNA_FILTERED_COUNTS = GSE_PROCESSED_DIR / "mirna_filtered_counts.tsv"
GSE_MIRNA_LOG2_CPM = GSE_PROCESSED_DIR / "mirna_log2_cpm.tsv"

GSE_TECHNICAL_SUMMARY = GSE_TABLE_DIR / "technical_count_summary.tsv"
GSE_FEATURE_AUDIT = GSE_TABLE_DIR / "feature_audit.tsv"
GSE_SAMPLE_AUDIT = GSE_TABLE_DIR / "sample_count_audit.tsv"
GSE_MIRNA_AUDIT = GSE_TABLE_DIR / "mirna_feature_audit.tsv"
GSE_SAMPLE_QC = GSE_TABLE_DIR / "mirna_sample_qc.tsv"
GSE_PCA_COORDINATES = GSE_TABLE_DIR / "pca_sample_coordinates.tsv"
GSE_PCA_VARIANCE = GSE_TABLE_DIR / "pca_explained_variance.tsv"
GSE_SAMPLE_CORRELATION = GSE_TABLE_DIR / "sample_pearson_correlation.tsv"
GSE_CLUSTERING_ORDER = GSE_TABLE_DIR / "sample_clustering_order.tsv"
GSE_DONOR_CENTERED = GSE_TABLE_DIR / "donor_centered_log2_cpm.tsv"
GSE_DONOR_PCA_COORDINATES = (
    GSE_TABLE_DIR / "donor_centered_pca_coordinates.tsv"
)
GSE_DONOR_PCA_VARIANCE = GSE_TABLE_DIR / "donor_centered_pca_variance.tsv"
GSE_DONOR_PCA_LOADINGS = GSE_TABLE_DIR / "donor_centered_pca_loadings.tsv"
GSE_WITHIN_DONOR_DISTANCES = (
    GSE_TABLE_DIR / "within_donor_method_distances.tsv"
)
GSE_WITHIN_DONOR_DISTANCE_SUMMARY = (
    GSE_TABLE_DIR / "within_donor_method_distance_summary.tsv"
)
GSE_METHOD_CENTROIDS = GSE_TABLE_DIR / "donor_centered_method_centroids.tsv"
GSE_ALL_PAIRED_EFFECTS = GSE_TABLE_DIR / "all_paired_feature_effects.tsv"
GSE_SELECTED_EFFECT_MATRIX = (
    GSE_TABLE_DIR / "selected_feature_effect_matrix.tsv"
)
GSE_PAIRED_OVERVIEW = GSE_TABLE_DIR / "paired_effect_contrast_overview.tsv"
GSE_CONSENSUS_SIGNATURE = (
    GSE_TABLE_DIR / "consensus_method_sensitive_signature.tsv"
)
GSE_METHOD_SENSITIVITY = (
    GSE_TABLE_DIR / "feature_method_sensitivity_summary.tsv"
)
GSE_EFFECT_MATRIX = GSE_TABLE_DIR / "feature_contrast_effect_matrix.tsv"
GSE_ROBUSTNESS_MATRIX = (
    GSE_TABLE_DIR / "feature_contrast_robustness_matrix.tsv"
)
GSE_EFFECTS_WITH_ROBUSTNESS = (
    GSE_TABLE_DIR / "all_paired_feature_effects_with_robustness.tsv"
)
GSE_ROBUST_DISTRIBUTION = (
    GSE_TABLE_DIR / "consensus_robust_contrast_distribution.tsv"
)
GSE_FEATURE_CLASSIFICATION = GSE_TABLE_DIR / "feature_classification.tsv"
GSE_FEATURE_CLASS_SUMMARY = GSE_TABLE_DIR / "feature_class_summary.tsv"
GSE_SENSITIVITY_TIER_SUMMARY = (
    GSE_TABLE_DIR / "method_sensitivity_tier_summary.tsv"
)
GSE_SUMMARY_COORDINATES = (
    GSE_TABLE_DIR / "summary_signature_coordinates.tsv"
)
GSE_SUMMARY_PNG = GSE_FIGURE_DIR / "summary_signature_scatter.png"
GSE_SUMMARY_PDF = GSE_FIGURE_DIR / "summary_signature_scatter.pdf"
GSE_FIGURE_LEGEND = GSE_TEXT_DIR / "Figure7_legend.md"


# GSE121910 salivary miRNA

GSE121910_DATASET = "GSE121910"

GSE121910_RAW_DIR = DATA_DIR / "raw" / GSE121910_DATASET
GSE121910_PROCESSED_DIR = DATA_DIR / "processed" / GSE121910_DATASET
GSE121910_RESULTS_DIR = RESULTS_DIR / GSE121910_DATASET
GSE121910_TABLE_DIR = GSE121910_RESULTS_DIR / "tables"
GSE121910_FIGURE_DIR = GSE121910_RESULTS_DIR / "figures"
GSE121910_LOG_DIR = LOG_DIR / GSE121910_DATASET
GSE121910_SCRIPT_DIR = SCRIPT_DIR / "gse121910"


# Cross-cohort integration

CROSS_COHORT_DIR = INTEGRATION_DIR / "cross_cohort_validation"
INTEGRATION_SCRIPT_DIR = SCRIPT_DIR / "integration"
