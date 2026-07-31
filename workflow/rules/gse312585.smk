# GSE312585 salivary extracellular-vesicle miRNA workflow.
# All commands are executed from the repository root.

rule gse312585_build_count_matrix:
    input:
        metadata=GSE_METADATA,
        script=GSE_SCRIPT_DIR / "00_build_count_matrix.py"
    output:
        counts=GSE_RAW_COUNT_MATRIX,
        technical=GSE_TECHNICAL_SUMMARY,
        feature_audit=GSE_FEATURE_AUDIT,
        sample_audit=GSE_SAMPLE_AUDIT
    log:
        GSE_LOG_DIR / "00_build_count_matrix.log"
    shell:
        """
        mkdir -p \
            {GSE_PROCESSED_DIR} \
            {GSE_TABLE_DIR} \
            {GSE_LOG_DIR}

        python {input.script} > {log} 2>&1
        """


rule gse312585_prepare_mirna_matrix:
    input:
        counts=GSE_RAW_COUNT_MATRIX,
        metadata=GSE_METADATA,
        feature_audit=GSE_FEATURE_AUDIT,
        script=GSE_SCRIPT_DIR / "01_prepare_mirna_matrix.py"
    output:
        raw=GSE_MIRNA_RAW_COUNTS,
        filtered=GSE_MIRNA_FILTERED_COUNTS,
        cpm=GSE_MIRNA_LOG2_CPM,
        audit=GSE_MIRNA_AUDIT,
        sample_qc=GSE_SAMPLE_QC
    log:
        GSE_LOG_DIR / "01_prepare_mirna_matrix.log"
    shell:
        """
        mkdir -p \
            {GSE_PROCESSED_DIR} \
            {GSE_TABLE_DIR} \
            {GSE_LOG_DIR}

        python {input.script} > {log} 2>&1
        """


rule gse312585_sample_structure_qc:
    input:
        expression=GSE_MIRNA_LOG2_CPM,
        metadata=GSE_METADATA,
        sample_qc=GSE_SAMPLE_QC,
        script=GSE_SCRIPT_DIR / "02_sample_structure_qc.py"
    output:
        coordinates=GSE_PCA_COORDINATES,
        variance=GSE_PCA_VARIANCE,
        correlation=GSE_SAMPLE_CORRELATION,
        clustering=GSE_CLUSTERING_ORDER,
        pca_method=GSE_FIGURE_DIR / "pca_method_and_donor.png",
        pca_donor=GSE_FIGURE_DIR / "pca_by_donor.png",
        correlation_plot=GSE_FIGURE_DIR / "sample_correlation_heatmap.png",
        clustering_plot=GSE_FIGURE_DIR / "sample_hierarchical_clustering.png",
        library_qc=GSE_FIGURE_DIR / "mirna_library_size_qc.png",
        feature_qc=GSE_FIGURE_DIR / "mirna_detected_features_qc.png"
    log:
        GSE_LOG_DIR / "02_sample_structure_qc.log"
    shell:
        """
        mkdir -p \
            {GSE_TABLE_DIR} \
            {GSE_FIGURE_DIR} \
            {GSE_LOG_DIR}

        python {input.script} > {log} 2>&1
        """


rule gse312585_donor_centered_method_effects:
    input:
        expression=GSE_MIRNA_LOG2_CPM,
        metadata=GSE_METADATA,
        script=GSE_SCRIPT_DIR / "03_donor_centered_method_effects.py"
    output:
        centered=GSE_DONOR_CENTERED,
        coordinates=GSE_DONOR_PCA_COORDINATES,
        variance=GSE_DONOR_PCA_VARIANCE,
        loadings=GSE_DONOR_PCA_LOADINGS,
        distances=GSE_WITHIN_DONOR_DISTANCES,
        distance_summary=GSE_WITHIN_DONOR_DISTANCE_SUMMARY,
        centroids=GSE_METHOD_CENTROIDS,
        pca=GSE_FIGURE_DIR / "donor_centered_pca.png",
        distance_plot=GSE_FIGURE_DIR / "within_donor_method_distances.png"
    log:
        GSE_LOG_DIR / "03_donor_centered_method_effects.log"
    shell:
        """
        mkdir -p \
            {GSE_TABLE_DIR} \
            {GSE_FIGURE_DIR} \
            {GSE_LOG_DIR}

        python {input.script} > {log} 2>&1
        """


rule gse312585_paired_feature_effects:
    input:
        expression=GSE_MIRNA_LOG2_CPM,
        metadata=GSE_METADATA,
        script=GSE_SCRIPT_DIR / "04_paired_feature_effects.py"
    output:
        effects=GSE_ALL_PAIRED_EFFECTS,
        selected=GSE_SELECTED_EFFECT_MATRIX,
        overview=GSE_PAIRED_OVERVIEW,
        heatmap=GSE_FIGURE_DIR / "paired_method_effect_heatmap.png"
    log:
        GSE_LOG_DIR / "04_paired_feature_effects.log"
    shell:
        """
        mkdir -p \
            {GSE_TABLE_DIR} \
            {GSE_FIGURE_DIR} \
            {GSE_LOG_DIR}

        python {input.script} > {log} 2>&1
        """


rule gse312585_consensus_method_signature:
    input:
        effects=GSE_ALL_PAIRED_EFFECTS,
        script=GSE_SCRIPT_DIR / "05_consensus_method_signature.py"
    output:
        signature=GSE_CONSENSUS_SIGNATURE,
        sensitivity=GSE_METHOD_SENSITIVITY,
        effects=GSE_EFFECT_MATRIX,
        robustness=GSE_ROBUSTNESS_MATRIX,
        annotated=GSE_EFFECTS_WITH_ROBUSTNESS,
        distribution=GSE_ROBUST_DISTRIBUTION,
        heatmap=GSE_FIGURE_DIR / "consensus_method_signature_heatmap.png",
        rank_plot=GSE_FIGURE_DIR / "consensus_signature_rank.png"
    log:
        GSE_LOG_DIR / "05_consensus_method_signature.log"
    shell:
        """
        mkdir -p \
            {GSE_TABLE_DIR} \
            {GSE_FIGURE_DIR} \
            {GSE_LOG_DIR}

        python {input.script} > {log} 2>&1
        """


rule gse312585_feature_classification:
    input:
        signature=GSE_CONSENSUS_SIGNATURE,
        script=GSE_SCRIPT_DIR / "06A_feature_classification.py"
    output:
        classification=GSE_FEATURE_CLASSIFICATION,
        class_summary=GSE_FEATURE_CLASS_SUMMARY,
        tier_summary=GSE_SENSITIVITY_TIER_SUMMARY,
        class_plot=GSE_FIGURE_DIR / "feature_class_distribution.png",
        top_plot=GSE_FIGURE_DIR / "top_consensus_features_by_class.png"
    log:
        GSE_LOG_DIR / "06A_feature_classification.log"
    shell:
        """
        mkdir -p \
            {GSE_TABLE_DIR} \
            {GSE_FIGURE_DIR} \
            {GSE_LOG_DIR}

        python {input.script} > {log} 2>&1
        """


rule gse312585_summary_signature_figure:
    input:
        classification=GSE_FEATURE_CLASSIFICATION,
        script=GSE_SCRIPT_DIR / "07_summary_signature_figure.py"
    output:
        png=GSE_SUMMARY_PNG,
        pdf=GSE_SUMMARY_PDF,
        coordinates=GSE_SUMMARY_COORDINATES,
        legend=GSE_FIGURE_LEGEND
    log:
        GSE_LOG_DIR / "07_summary_signature_figure.log"
    shell:
        """
        mkdir -p \
            {GSE_TABLE_DIR} \
            {GSE_FIGURE_DIR} \
            {GSE_TEXT_DIR} \
            {GSE_LOG_DIR}

        python {input.script} > {log} 2>&1
        """
