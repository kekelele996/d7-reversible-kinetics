# D.7 可逆一级反应参数反演

题目包入口为 `python3 /app/main.py <batch.json>`。环境只使用 Python 3.12 标准库；
测试依赖以 wheel 形式放在 `tests/vendor/`，不进入运行镜像。

本地验证命令：

```bash
harbor run -p . -a oracle
harbor run -p . -a nop
```

交付提交、推送远端和登记表格属于对外操作，需由提交人明确授权。
