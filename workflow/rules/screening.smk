rule score_candidate_datasets:
    input:
        registry=DATASET_REGISTRY,
        script=SCRIPT_DIR / "00_score_candidate_datasets.py"
    output:
        prioritisation=DATASET_PRIORITISATION,
        summary=DATASET_SCREENING_SUMMARY
    log:
        LOG_DIR / "00_score_candidate_datasets.log"
    shell:
        """
        mkdir -p {SCREENING_DIR} {LOG_DIR}
        python {input.script} > {log} 2>&1
        """
