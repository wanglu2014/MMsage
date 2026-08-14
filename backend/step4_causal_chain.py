"""
Step 4: Causal Chain Generation
================================
Generate causal chain text using LLM API.
Produce final results JSON and markdown report.
"""

import os
import json
from pathlib import Path
from typing import Optional, Dict, Any, List

try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False


# LLM Configuration (set OPENAI_API_KEY / OPENAI_API_BASE / OPENAI_MODEL in environment)
LLM_API_KEY = os.getenv("OPENAI_API_KEY", "") or os.getenv("AI2API_KEY", "")
LLM_API_BASE = os.getenv("OPENAI_API_BASE", "https://api.openai.com/v1")
LLM_MODEL = os.getenv("OPENAI_MODEL", "gpt-4o")


SYSTEM_PROMPT = """You are an expert microbiome researcher. Given a microbe-metabolite pair with trajectory analysis scores and knowledge graph paths, generate a concise causal chain explaining the biological mechanism.

Output JSON with:
{
  "causal_chain": "One-sentence causal chain: Bacteria (releases/absorbs) Metabolite, through X mechanism affecting IBD",
  "mechanism_detail": "2-3 sentences expanding the mechanism",
  "evidence_pmids": ["PMID1", "PMID2"],
  "kegg_pathways": ["hsa00640", "hsa00280"],
  "confidence": "high/medium/low"
}"""


def build_context(item: Dict[str, Any]) -> str:
    """Build context string from quadrant result for LLM prompt."""
    cand = item.get('scoresp_result', {}).get('candidate', {})
    bacteria = cand.get('bacteria', 'Unknown')
    metabolite = cand.get('metabolite', 'Unknown')
    mms = item.get('mmsage_norm', 0)
    scoresp = item.get('score_sp', 0)
    quadrant = item.get('quadrant', '?')
    direction = item.get('direction', 'releases')
    path_str = item.get('scoresp_result', {}).get('path_str', 'No path')
    is_dark = item.get('is_dark_matter', False)

    path_nodes = item.get('scoresp_result', {}).get('path_nodes', [])
    path_edges = item.get('scoresp_result', {}).get('path_edges', [])

    ctx = f"""Microbe: {bacteria}
Metabolite: {metabolite}
Direction: {direction}
MMSage Signal: {mms:.4f}
ScoreSP Energy: {scoresp:.4f}
Quadrant: {quadrant} ({'Dark Matter - novel discovery' if is_dark else item.get('quadrant_label', '')})
Knowledge Graph Path: {path_str}
"""
    if path_nodes:
        ctx += "\nPath Node Details:\n"
        for n in path_nodes:
            ctx += f"  - {n['name']} (type={n['type']}, novelty={n['novelty']}, cred={n['cred']})\n"
    if path_edges:
        ctx += "\nPath Edge Details:\n"
        for e in path_edges:
            ctx += f"  - {e['source']} -> {e['target']} (support={e['support']}, resistance={e['r_e']})\n"

    return ctx


def generate_causal_chain_llm(context: str, bacteria: str, metabolite: str) -> Dict[str, Any]:
    """Call LLM API to generate causal chain."""
    if not HAS_REQUESTS or not LLM_API_KEY:
        raise RuntimeError("LLM API not available")

    headers = {
        "Content-Type": "application/json",
        "Authorization": f"Bearer {LLM_API_KEY}",
    }
    payload = {
        "model": LLM_MODEL,
        "messages": [
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": f"Analyze this microbe-metabolite pair:\n\n{context}"},
        ],
        "temperature": 0.7,
        "max_tokens": 1000,
        "response_format": {"type": "json_object"},
    }

    resp = requests.post(
        f"{LLM_API_BASE}/chat/completions",
        headers=headers,
        json=payload,
        timeout=60,
    )
    if resp.status_code != 200:
        raise RuntimeError(f"LLM API error: {resp.status_code}")

    content = resp.json()['choices'][0]['message']['content']
    return json.loads(content)



def generate_report_markdown(results: List[Dict[str, Any]]) -> str:
    """Generate markdown report table."""
    lines = [
        "# MMSage + ScoreSP Dual-Axis Analysis Report\n",
        f"Total candidates analyzed: {len(results)}\n",
    ]

    # Quadrant summary
    q_counts = {}
    for r in results:
        q = r.get('quadrant', '?')
        q_counts[q] = q_counts.get(q, 0) + 1

    lines.append("## Quadrant Summary\n")
    lines.append("| Quadrant | Label | Count |")
    lines.append("|----------|-------|-------|")
    labels = {'I': 'High Priority', 'II': 'Dark Matter', 'III': 'Need More Data', 'IV': 'Low Priority'}
    for q in ['I', 'II', 'III', 'IV']:
        marker = ' *' if q == 'II' else ''
        lines.append(f"| {q}{marker} | {labels.get(q, '')} | {q_counts.get(q, 0)} |")
    lines.append("")

    # Results table
    lines.append("## Top Candidates\n")
    lines.append("| Rank | Bacteria | Direction | Metabolite | MMSage | ScoreSP | Quadrant | Causal Chain |")
    lines.append("|------|----------|-----------|------------|--------|---------|----------|-------------|")

    for i, r in enumerate(results, 1):
        cand = r.get('scoresp_result', {}).get('candidate', {})
        bacteria = cand.get('bacteria', '?')
        metabolite = cand.get('metabolite', '?')
        mms = r.get('mmsage_norm', 0)
        scoresp = r.get('score_sp', 0)
        quad = r.get('quadrant', '?')
        direction = r.get('direction', '?')
        chain = r.get('causal_chain_result', {}).get('causal_chain', 'N/A')

        # Add star marker for dark matter
        quad_display = f"{quad}" + ("*" if r.get('is_dark_matter') else "")

        # Truncate chain for table
        if len(chain) > 80:
            chain = chain[:77] + "..."

        lines.append(
            f"| {i} | {bacteria} | {direction} | {metabolite} | "
            f"{mms:.3f} | {scoresp:.3f} | {quad_display} | {chain} |"
        )

    lines.append("")
    return "\n".join(lines)


def run_step4(
    quadrant_path: str,
    output_path: str,
    report_path: str,
    use_llm: bool = True,
    max_llm_calls: int = 20,
) -> list:
    """
    Run step 4: generate causal chains for all candidates.

    Args:
        quadrant_path: Path to step3 output JSON
        output_path: Path to save step4 output JSON
        report_path: Path to save markdown report
        use_llm: Whether to attempt LLM API calls
        max_llm_calls: Maximum number of LLM API calls (to control cost)

    Returns:
        List of final result dicts
    """
    print("=" * 60)
    print("  STEP 4: Causal Chain Generation")
    print("=" * 60)

    with open(quadrant_path, 'r', encoding='utf-8') as f:
        quadrant_results = json.load(f)

    if not HAS_REQUESTS:
        raise RuntimeError("requests library required. Install with: pip install requests")

    results = []
    llm_calls = 0
    llm_success = 0
    llm_failed = 0
    skipped = 0

    for i, item in enumerate(quadrant_results):
        cand = item.get('scoresp_result', {}).get('candidate', {})
        bacteria = cand.get('bacteria', 'Unknown')
        metabolite = cand.get('metabolite', 'Unknown')

        causal_result = None

        if use_llm and llm_calls < max_llm_calls:
            quadrant = item.get('quadrant', 'IV')
            if quadrant in ('I', 'II'):
                try:
                    context = build_context(item)
                    causal_result = generate_causal_chain_llm(context, bacteria, metabolite)
                    causal_result['generated_by'] = 'llm'
                    llm_calls += 1
                    llm_success += 1
                except Exception as e:
                    print(f"  [ERROR] LLM failed for {bacteria} x {metabolite}: {e}")
                    llm_calls += 1
                    llm_failed += 1

        if causal_result is None:
            skipped += 1
            item['causal_chain_result'] = {
                'causal_chain': None,
                'generated_by': 'skipped',
            }
        else:
            item['causal_chain_result'] = causal_result

        results.append(item)

        if (i + 1) % 50 == 0:
            print(f"  Processed {i + 1}/{len(quadrant_results)}")

    print(f"\nCausal chain generation: LLM OK={llm_success}, "
          f"LLM failed={llm_failed}, skipped={skipped}")
    if llm_failed > 0 and llm_success == 0:
        print(f"  [ERROR] All LLM calls failed. Check API key/proxy.")

    # Save results JSON
    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"Saved {len(results)} results to {output_path}")

    # Generate and save report
    report = generate_report_markdown(results)
    report_file = Path(report_path)
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write(report)
    print(f"Saved report to {report_path}")

    return results


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Step 4: Causal Chain Generation')
    parser.add_argument('--quadrant', default='outputs/step3_quadrant.json')
    parser.add_argument('--output', default='outputs/step4_final.json')
    parser.add_argument('--report', default='outputs/step4_report.md')
    parser.add_argument('--no-llm', action='store_true', help='Skip LLM API calls')
    parser.add_argument('--max-llm-calls', type=int, default=20)
    args = parser.parse_args()

    run_step4(args.quadrant, args.output, args.report,
              use_llm=not args.no_llm, max_llm_calls=args.max_llm_calls)
