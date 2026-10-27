#!/bin/zsh
cd /Users/leighpierce/Git/246
export NL_WARP2=1 NL_BK=1 NL_BK_ALPHA=${1:-1} NL_BK_FILE=/Users/leighpierce/Git/246/analysis/nl/branch_kernel_warp2.json
.venv/bin/python tools/memguard.py --limit-gb 6 -- .venv/bin/python tools/nl_model.py --stage bench --workers 4 > analysis/nl_hfdrift/bench_warp2_bk_kf_a$NL_BK_ALPHA.log 2>&1
echo finished > analysis/nl_hfdrift/bench_kf_a$NL_BK_ALPHA.done
