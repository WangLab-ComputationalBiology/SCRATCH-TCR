process SCIRPY_QUALITY {

    tag "scTCR repertoire quality"
    label 'process_medium'

    container "oandrefonseca/scratch-tcr:main"

    input:
        path(notebook)
        path(path_vdj_folder)
        path(exp_table)
        path(config)

    output:
        path("data/${params.project_name}_tcr_repertoire_object.h5ad")   , emit: anndata
        path("report/${notebook.baseName}.html")                         , emit: html
        path("_freeze/**/figure-html/*.png")                             , emit: figures

    when:
        task.ext.when == null || task.ext.when

    shell:
        def param_file = task.ext.args ? "-P path_vdj_folder:'${path_vdj_folder.join(';')}' -P meta_data_csv:${exp_table} -P ${task.ext.args}" : ""
        """
        quarto render ${notebook} -P ${param_file}
        """

    stub:
        def param_file = task.ext.args ? "-P path_vdj_folder:'${path_vdj_folder.join(';')}' -P meta_data_csv:${exp_table} -P ${task.ext.args}" : ""
        """
        mkdir -p data _freeze/${notebook.baseName}
        mkdir -p _freeze/DUMMY/figure-html

        touch _freeze/DUMMY/figure-html/FILE.png

        touch data/${params.project_name}_tcr_repertoire_object.h5ad

        mkdir -p report
        touch report/${notebook.baseName}.html
        echo ${param_file} > _freeze/params.yml

        """
}