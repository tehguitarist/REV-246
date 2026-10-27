#!/bin/zsh
# Full warp2 benches, one after the other (memguard 6 GB, 4 workers). Logs in analysis/nl_hfdrift/.
cd /Users/leighpierce/Git/246
export NL_WARP2=1
unset NL_BK NL_BK_ALPHA NL_BK_FILE
.venv/bin/python tools/memguard.py --limit-gb 6 -- .venv/bin/python tools/nl_model.py --stage bench --workers 4 > analysis/nl_hfdrift/bench_warp2.log 2>&1
NL_BK=1 .venv/bin/python tools/memguard.py --limit-gb 6 -- .venv/bin/python tools/nl_model.py --stage bench --workers 4 > analysis/nl_hfdrift/bench_warp2_bk.log 2>&1
echo finished > analysis/nl_hfdrift/benches.done
