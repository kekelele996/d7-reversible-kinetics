#!/usr/bin/env bash
# 参考解：用 solution/app/ 下的完整实现替换 /app 里同名文件。
# 每个难点的改法、依据和独立复算写在 NOTES.md；这里只负责安装和留一法开关。
#
# 留一法：HTB_SKIP=<难点 id>[,<难点 id>...] 时把这些难点的修复退回原写法，其余照装。
#   harbor run -p <包> -a oracle --ae HTB_SKIP=difficulty_1
# 正式 Oracle 不设这个变量。退回用精确片段替换，每段必须恰好命中一次，否则直接失败，
# 免得参考解改了而留一法悄悄失效。
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
install -m 0644 "$HERE/app/main.py" /app/main.py

SKIP="${HTB_SKIP:-}" /usr/local/bin/python3 - /app/main.py <<'PY'
import os
import sys
from pathlib import Path

REVERT = {
    "difficulty_1": (
        '    response_base = case["blank_absorbance"] + path * concentration * epsilon_b\n'
        "    response_span = path * concentration * (epsilon_a - epsilon_b)",
        '    response_base = case["blank_absorbance"]\n'
        "    response_span = path * concentration * epsilon_a",
    ),
    "difficulty_2": (
        "    reverse_rate = relaxation * equilibrium\n"
        "    forward_rate = relaxation - reverse_rate",
        "    reverse_rate = relaxation - relaxation * equilibrium\n"
        "    forward_rate = relaxation - reverse_rate",
    ),
    "difficulty_3": (
        "    reported_dead_time = dead_time",
        "    reported_dead_time = 0.0",
    ),
}
skip = [item for item in os.environ.get("SKIP", "").split(",") if item]
unknown = sorted(set(skip) - set(REVERT))
if unknown:
    sys.exit(f"HTB_SKIP 里有未知难点：{unknown}")
path = Path(sys.argv[1])
text = path.read_text(encoding="utf-8")
for key in skip:
    fixed, original = REVERT[key]
    if text.count(fixed) != 1:
        sys.exit(f"{key} 的修复片段在参考实现里出现 {text.count(fixed)} 次，应为 1")
    text = text.replace(fixed, original)
    print(f"留一法：{key} 退回原写法", file=sys.stderr)
path.write_text(text, encoding="utf-8")
PY
chown runner:runner /app/main.py 2>/dev/null || true
