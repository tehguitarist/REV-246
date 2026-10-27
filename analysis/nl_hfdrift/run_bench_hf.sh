#!/bin/zsh
cd /Users/leighpierce/Git/246
export NL_BK=1 NL_HF=1
.venv/bin/python tools/memguard.py --limit-gb 6 -- .venv/bin/python tools/nl_model.py --stage bench --workers 4 > analysis/nl_hfdrift/bench_hf.log 2>&1
echo finished > analysis/nl_hfdrift/bench_hf.done
