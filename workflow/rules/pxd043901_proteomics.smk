rule download_pxd043901:
    input:
        script=PXD_SCRIPT_DIR / "download_pxd043901.py"
    output:
        marker=RAW_DOWNLOAD_MARKER
    log:
        PXD_LOG_DIR / "00_download_pxd043901.log"
    shell:
        """
        mkdir -p {PXD_RAW_DIR} {PXD_LOG_DIR} {MARKER_DIR}
        python {input.script} > {log} 2>&1
        test -d {PXD_RAW_DIR}
        test "$(find {PXD_RAW_DIR} -type f | wc -l)" -gt 0
        touch {output.marker}
        """

rule build_sample_metadata:
    input:
        download_complete=RAW_DOWNLOAD_MARKER,
        script=PXD_SCRIPT_DIR / "01_build_sample_metadata.py"
    output:
        metadata=METADATA
    log:
        PXD_LOG_DIR / "01_build_sample_metadata.log"
    shell:
        """
        mkdir -p {PXD_PROCESSED_DIR} {PXD_LOG_DIR}
        python {input.script} > {log} 2>&1
        """

rule create_clean_sample_metadata:
    input:
        metadata=METADATA
    output:
        metadata_clean=METADATA_CLEAN
    log:
        PXD_LOG_DIR / "01b_create_clean_sample_metadata.log"
    shell:
        """
        mkdir -p {PXD_PROCESSED_DIR} {PXD_LOG_DIR}
        cp {input.metadata} {output.metadata_clean}
        printf 'Compatibility copy created from %s\n' "{input.metadata}" > {log}
        """

rule prepare_protein_matrix:
    input:
        download_complete=RAW_DOWNLOAD_MARKER,
        metadata=METADATA,
        script=PXD_SCRIPT_DIR / "02_prepare_protein_matrix.py"
    output:
        matrix_raw=LFQ_RAW,
        matrix_log2=LFQ_LOG2,
        annotation=PROTEIN_ANNOTATION,
        summary=MATRIX_SUMMARY
    log:
        PXD_LOG_DIR / "02_prepare_protein_matrix.log"
    shell:
        """
        mkdir -p {PXD_PROCESSED_DIR} {PXD_LOG_DIR}
        python {input.script} > {log} 2>&1
        """

rule annotate_protein_origin:
    input:
        annotation=PROTEIN_ANNOTATION,
        matrix=LFQ_LOG2,
        script=PXD_SCRIPT_DIR / "02b_annotate_protein_origin.py"
    output:
        extended_annotation=PROTEIN_ANNOTATION_EXTENDED,
        organism_summary=ORGANISM_SUMMARY
    log:
        PXD_LOG_DIR / "02b_annotate_protein_origin.log"
    shell:
        """
        mkdir -p {PXD_PROCESSED_DIR} {PXD_LOG_DIR}
        python {input.script} > {log} 2>&1
        """

rule run_quality_control:
    input:
        matrix=LFQ_LOG2,
        metadata=METADATA,
        annotation=PROTEIN_ANNOTATION_EXTENDED,
        script=PXD_SCRIPT_DIR / "03_run_quality_control.py"
    output:
        sample_metrics=QC_SAMPLE_METRICS,
        protein_metrics=QC_PROTEIN_METRICS,
        correlation=QC_CORRELATION,
        pca_coordinates=QC_PCA_COORDINATES,
        pca_variance=QC_PCA_VARIANCE,
        summary=QC_SUMMARY,
        detected=PXD_QC_FIGURE_DIR / "detected_proteins_per_sample.png",
        missing=PXD_QC_FIGURE_DIR / "missing_fraction_per_sample.png",
        boxplots=PXD_QC_FIGURE_DIR / "lfq_intensity_boxplots.png",
        heatmap=PXD_QC_FIGURE_DIR / "sample_correlation_heatmap.png",
        pca=PXD_QC_FIGURE_DIR / "pca_samples.png"
    log:
        PXD_LOG_DIR / "03_run_quality_control.log"
    shell:
        """
        mkdir -p {PXD_QC_DIR} {PXD_QC_FIGURE_DIR} {PXD_LOG_DIR}
        python {input.script} > {log} 2>&1
        """

rule calculate_patient_responses:
    input:
        matrix=LFQ_LOG2,
        metadata=METADATA,
        annotation=PROTEIN_ANNOTATION_EXTENDED,
        qc_summary=QC_SUMMARY,
        script=PXD_SCRIPT_DIR / "04_calculate_patient_responses.py"
    output:
        delta_t1=DELTA_T1,
        delta_t6=DELTA_T6,
        response_table=PROTEIN_RESPONSE_SUMMARY,
        preliminary=ROBUST_PRELIMINARY,
        summary=RESPONSE_SUMMARY
    log:
        PXD_LOG_DIR / "04_calculate_patient_responses.log"
    shell:
        """
        mkdir -p {PXD_TRAJECTORY_DIR} {PXD_LOG_DIR}
        python {input.script} > {log} 2>&1
        """

rule test_longitudinal_responses:
    input:
        delta_t1=DELTA_T1,
        delta_t6=DELTA_T6,
        annotation=PROTEIN_ANNOTATION_EXTENDED,
        preliminary=ROBUST_PRELIMINARY,
        script=PXD_SCRIPT_DIR / "04b_test_longitudinal_responses.py"
    output:
        tests=LONGITUDINAL_TESTS,
        stringent=STRINGENT_RESPONSES,
        ranked=RANKED_RESPONSES,
        summary=STATISTICAL_TEST_SUMMARY
    log:
        PXD_LOG_DIR / "04b_test_longitudinal_responses.log"
    shell:
        """
        mkdir -p {PXD_TRAJECTORY_DIR} {PXD_LOG_DIR}
        python {input.script} > {log} 2>&1
        """

rule cluster_protein_trajectories:
    input:
        tests=LONGITUDINAL_TESTS,
        delta_t1=DELTA_T1,
        delta_t6=DELTA_T6,
        annotation=PROTEIN_ANNOTATION_EXTENDED,
        script=PXD_SCRIPT_DIR / "05_cluster_protein_trajectories.py"
    output:
        classes=TRAJECTORY_CLASSES,
        persistent=PERSISTENT_RESPONSES,
        counts=TRAJECTORY_COUNTS,
        summary=TRAJECTORY_SUMMARY
    log:
        PXD_LOG_DIR / "05_cluster_protein_trajectories.log"
    shell:
        """
        mkdir -p {PXD_TRAJECTORY_DIR} {PXD_LOG_DIR}
        python {input.script} > {log} 2>&1
        """
