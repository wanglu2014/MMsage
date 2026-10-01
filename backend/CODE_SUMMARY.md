# TrajMM Pipeline - Code Summary

## 本次更新内容

### 1. Bug 修复

| 文件 | 问题 | 修复 |
|------|------|------|
| `agents/pubmed_query.py` | `DEFAULT_SYNONYMS` 未定义导致 import 失败 | 添加空字典 `DEFAULT_SYNONYMS: Dict[str, List[str]] = {}` 作为兼容导出 |
| `step2_scoresp_energy.py` | 硬编码酸/盐同义词仅覆盖5对 metabolite | 改为动态 `_acid_ate_variants()` 规则，自动处理所有 `-ic_acid↔-ate`、`-ous_acid↔-ite` 转换 |
| `step2_chain_novelty.py`<br>`agents/pubmed_query.py` | IBD 硬编码关键词 (`"ibd"`, `"colitis"`, `"crohn"`) | 改为通用疾病后缀检测 (`"disease"`, `"syndrome"`, `"itis"`, `"osis"`, `"emia"` 等) |
| `db_checkers/binded_checker.py` | SapBERT 每次调用都尝试导入 torch 然后失败，浪费大量时间 | 添加模块级标志 `_sapbert_available`，失败一次后永久跳过 |
| `data/knowledge_graph/auto_built_kg.gml` | 节点缺少 `label` 属性 | 后处理添加，所有 1320 个节点补全 label |

### 2. Pipeline 运行结果

```
Step 1 (TrajMM Signal):    146 candidates, 0.0s
Step 2 (Chain Novelty):    146 scored, 7 with KG paths, 109 dark matter, 0.5s
Step 2b (Multi-Agent):     146 agents, 7366.4s (主要耗时)
Step 3 (Quadrant):         146 results, 0.6s

Quadrant 分布:
  I (高优先级):  55
  II (暗物质):   54
  III (已知):    18
  IV (低优先):   19
  ────────────────
  总计:         146
```

### 3. 关键输出文件

| 文件 | 说明 |
|------|------|
| `outputs/step1_candidates.json` | TrajMM 信号提取的候选对 |
| `outputs/step2_chain_novelty.json` | Chain Novelty 评分结果 |
| `outputs/step2b_agent_evidence.json` | Multi-Agent 证据聚合 |
| `outputs/step3_quadrant.json` | 四象限分配结果 |
| `data/knowledge_graph/auto_built_kg.gml` | 1320节点, 2130边, 97%边有OpenAlex增强 |

### 4. 验证状态

- ✅ 数据完整性: 146个候选对全链路一致
- ✅ 四象限求和: 55+54+18+19=146
- ✅ 评分逻辑: chain_count=0 → chain_novelty=1.0
- ✅ KG路径: 7个has_path=True全部验证存在
- ✅ Agent公式: chain_novelty = 1 - log(1+total)/log(1+500) 精确匹配

### 5. 运行方式

```bash
cd E:/Onedrive/Apps/Research/TrajMM_nov/backend

# 完整 pipeline
python run_pipeline.py

# 或分步运行
python step1_trajmm_signal.py
python step2_chain_novelty.py
python step3_quadrant.py
```

### 6. 依赖

- Python: `biopython`, `openai`, `networkx`, `pandas`, `numpy`
- R: 4.4.1 (用于 TrajMM 信号计算)
- 可选: `torch` (SapBERT 语义编码，当前未安装，使用索引+回退)
