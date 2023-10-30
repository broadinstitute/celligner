sparkles sub -n celligner-run \
    -u celligner:celligner \
    -u mnnpy:mnnpy \
    -u install_submodules_and_run.sh \
    -u run_celligner.py \
    -u run_celligner_rollback.py \
    -u inputs.json \
    -u $HOME/.taiga/token:.taiga-token \
    sh install_submodules_and_run.sh "run_celligner.py"