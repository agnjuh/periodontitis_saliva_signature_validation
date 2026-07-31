rule calculate_module_scores:
    input:
        matrix=LFQ_LOG2,
        annotation=PROTEIN_ANNOTATION_EXTENDED,
        metadata=METADATA,
        modules=MODULE_DEFINITION,
        trajectories=TRAJECTORY_CLASSES,
        script=PXD_SCRIPT_DIR / "06_calculate_module_scores.py"
    output:
        scores=MODULE_SCORES,
        coverage=MODULE_COVERAGE,
        deltas=MODULE_PATIENT_DELTAS,
        statistics=MODULE_STATISTICS,
        summary=MODULE_SUMMARY
    log:
        PXD_LOG_DIR / "06_calculate_module_scores.log"
    shell:
        """
        mkdir -p {PXD_MODULE_DIR} {PXD_LOG_DIR}
        python {input.script} > {log} 2>&1
        """

rule module_sensitivity_analysis:
    input:
        matrix=LFQ_LOG2,
        annotation=PROTEIN_ANNOTATION_EXTENDED,
        metadata=METADATA,
        modules=MODULE_DEFINITION,
        module_scores=MODULE_SCORES,
        module_statistics=MODULE_STATISTICS,
        script=PXD_SCRIPT_DIR / "06b_module_sensitivity_analysis.py"
    output:
        results=MODULE_LOO_RESULTS,
        summary_table=MODULE_LOO_SUMMARY,
        correlations=MODULE_DELTA_CORRELATIONS,
        high_correlations=MODULE_HIGH_CORRELATIONS,
        summary=MODULE_SENSITIVITY_SUMMARY
    log:
        PXD_LOG_DIR / "06b_module_sensitivity_analysis.log"
    shell:
        """
        mkdir -p {PXD_SENSITIVITY_DIR} {PXD_LOG_DIR}
        python {input.script} > {log} 2>&1
        """
