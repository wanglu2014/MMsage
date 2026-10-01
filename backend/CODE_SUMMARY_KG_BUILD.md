# KG Build Pipeline — Code Summary for Claude Code

## 目标

从 PubMed 下载 "Akkermansia muciniphila" + IBD 全部实验文献 → DeepSeek LLM 提取三元组 → 生成 GML 知识图谱 → OpenAlex 增强 → 输出到 TrajMM_nov pipeline。

## 当前状态

- **已完成**: 上一轮用 6 个分散查询跑了 279 篇（272 成功），生成了 1437 节点 2005 边的 GML，已验证能被 TrajMM_nov 的 knowledge_pump/ScoreSP/Chain Novelty 正确加载
- **未完成**: 需要用精确单条查询重跑 303 篇（上一轮覆盖不全），然后跑 OpenAlex 增强
- **问题**: 上次精确查询运行时 DeepSeek API 连续 Connection error，需要重试

## 关键文件

### 1. `E:/Onedrive/Apps/Research/TrajMM_nov/backend/build_kg.py`（主代码，已修改可用）

核心流程:
- `DeepSeekPool`: 多 key 轮转调用 DeepSeek API（从 `E:/Onedrive/LabYYc/chatgpt/keys_0605.csv` 读 key）
- `search_pubmed()` / `fetch_abstracts()`: Biopython Entrez 搜索和下载
- `recognize_entities()`: 在摘要中识别已知实体（细菌/代谢物/疾病）
- `extract_relations_llm()`: 调 DeepSeek 提取 ===BEGIN NODES/EDGES=== CSV 格式三元组
- `_parse_csv_sections()`: 解析 LLM 返回的 CSV
- `_add_to_graph()`: 合并到 NetworkX DiGraph
- `build_kg_from_abstracts()`: 批量处理，带缓存（`cache/kg_build/extractions/{pmid}.json`）

**本轮做的 2 个 bug fix**:
1. `_add_to_graph()` 补充了 `description=node.get("description", "")` —— 原来漏写了节点描述
2. `extract_relations_llm()` 门槛从 `entity_lines < 2` 改为 `< 1` —— 原来只识别到 1 种实体就跳过，导致 260/279 "失败"
3. `expand_entity_names()` 新增 IBD 疾病同义词 fallback —— 原来 SapBERT 无 torch 时只有 "IBD" 一个词

### 2. `E:/Onedrive/Apps/Research/TrajMM_nov/backend/_run_precise.py`（精确查询运行脚本）

用单条 PubMed 查询精确获取 303 篇文献:

```
"Akkermansia muciniphila"
AND ("IBD" OR "inflammatory bowel disease" OR "ulcerative colitis"
     OR "Crohn's disease" OR "Crohn disease" OR "colitis" OR "enteritis"
     OR "ileitis" OR "intestinal inflammation" OR "gut inflammation"
     OR "mucosal inflammation")
AND ("in vivo" OR "in vitro" OR "cell line" OR "clinical trial"
     OR "experiment" OR "mouse" OR "mice" OR "rat" OR "rats"
     OR "murine" OR "animal model" OR "patient" OR "patients"
     OR "cohort" OR "randomized" OR "biopsy" OR "fecal"
     OR "stool" OR "culture" OR "fermentation")
```

- 已有 ~56 篇缓存（上一轮重叠的），剩余 ~245 篇需 LLM 处理
- 输出到: `E:/Onedrive/Apps/Research/TrajMM_nov/data/knowledge_graph/auto_built_kg.gml`

### 3. OpenAlex 增强（尚未执行）

参考代码: `E:/Onedrive/LabYYc/chatgpt/AgentClin/rep1221_clin/add_metrics_to_gml.py` 中的 `GMLMetricsEnricher`

流程:
1. 读取 GML，收集所有边的 unique PMID（边的 pmid 字段可能是逗号分隔的多个 PMID）
2. 对每个 PMID 调 OpenAlex API:
   - `https://api.openalex.org/works/pmid:{pmid}` → `cited_by_count`, `publication_year`
   - `https://api.openalex.org/sources/{source_id}` → `summary_stats.2yr_mean_citedness`（即 impact factor）
3. 写入边属性

**关键**: 属性名必须用 `edge_` 前缀对齐 TrajMM_nov:
- `edge_impact_factor`（不是 `impact_factor`）
- `edge_citation_count`（不是 `citation_count`）
- `publication_year`（原版 add_metrics_to_gml.py 没取这个字段，需要补）

原版 `add_metrics_to_gml.py` 用 `G.edges(keys=True, data=True)` 遍历（MultiDiGraph），但当前 GML 是 DiGraph，需要改用 `G.edges(data=True)`。

## GML 输出格式要求（TrajMM_nov 全链路）

### 节点
```
node [
  id 0
  label "akkermansia_muciniphila"    ← snake_case，nx.read_gml 默认用 label 做节点 key
  type "microbe"                     ← microbe/metabolite/pathway/disease
  node_type "microbe"                ← 同 type（冗余但 api_server 需要）
  description "..."                  ← LLM 生成的描述文本
]
```

### 边（增强前）
```
edge [
  source 0
  target 1
  relation "produces"                ← produces/metabolizes/inhibits/promotes/protects/associated_with/degrades/transports/modulates
  description "..."                  ← LLM 生成的证据描述
  pmid "32015508,30115164"           ← 逗号分隔多个 PMID
  journal "Nature Medicine"
  year "2020"
]
```

### 边（增强后，ScoreSP 需要）
```
edge [
  ... 同上 ...
  edge_impact_factor 82.9            ← OpenAlex 2yr_mean_citedness
  edge_citation_count 1256           ← OpenAlex cited_by_count
  publication_year 2020              ← OpenAlex publication_year
]
```

## 下游消费端（TrajMM_nov）

| 文件 | 读什么 |
|---|---|
| `knowledge_pump.py` | label, type, description, relation, pmid |
| `step2_scoresp_energy.py` | type, label, **edge_impact_factor**, **edge_citation_count**, **publication_year** |
| `step2_chain_novelty.py` | 通过 knowledge_pump 路径查找 |
| `api_server.py` | label, type, relation, edge_impact_factor, edge_citation_count |
| `agents/master_agent.py` | 通过 knowledge_pump.build_llm_context() |

## 缓存位置

- LLM 提取结果: `E:/Onedrive/Apps/Research/TrajMM_nov/cache/kg_build/extractions/{pmid}.json`
- SapBERT 同义词: `E:/Onedrive/Apps/Research/TrajMM_nov/cache/kg_build/sapbert_synonym_cache.json`
- DeepSeek API keys: `E:/Onedrive/LabYYc/chatgpt/keys_0605.csv`

## 运行环境

- Python 依赖: `biopython`, `openai`, `networkx`（已安装）
- `torch` 未安装（SapBERT 编码不可用，但索引查找和 fallback 同义词可用）
- 工作目录: `E:/Onedrive/Apps/Research/TrajMM_nov/backend`
- Windows PowerShell，stdout 需要 line_buffering=True 才能实时看到输出
