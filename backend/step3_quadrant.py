"""
Step 3: Dual-Axis Quadrant Assignment
=======================================
Assign each candidate to a quadrant based on
MMSage signal strength (X) and Chain Novelty (Y).

Quadrant layout (Y axis = Chain Novelty, high = top):

  Chain Novelty (high)
     |
     | II. New but weak data     I. Dark Matter ★
     | MMSage weak + Nov high    MMSage strong + Nov high
     | -> needs more experiments  -> HIGHEST priority
     |─────────────────────────────────────────────
     | IV. Low priority           III. Known relationship
     | Both weak                  MMSage strong + Nov low
     | -> ignore                  -> validate / review
     |________________________________________________ MMSage (strong)
   weak                                             strong

Thresholds are adaptive: median of each axis within current batch.
"""

import json
import math
from pathlib import Path

QUADRANT_LABELS = {
    'I': 'Dark Matter',
    'II': 'Novel but Weak Signal',
    'III': 'Known Relationship',
    'IV': 'Low Priority',
}


def assign_quadrant(mmsage_norm: float, chain_novelty: float,
                    mmsage_threshold: float, novelty_threshold: float) -> dict:
    """
    Assign a candidate to a quadrant.

    Quadrant logic:
      I:   mmsage_norm > threshold AND chain_novelty >= threshold  (Dark Matter ★)
      II:  mmsage_norm <= threshold AND chain_novelty >= threshold  (Novel but Weak Signal)
      III: mmsage_norm > threshold AND chain_novelty < threshold    (Known Relationship)
      IV:  mmsage_norm <= threshold AND chain_novelty < threshold   (Low Priority)
    """
    high_mmsage = mmsage_norm > mmsage_threshold
    high_novelty = chain_novelty >= novelty_threshold

    if high_mmsage and high_novelty:
        quadrant = 'I'
    elif not high_mmsage and high_novelty:
        quadrant = 'II'
    elif high_mmsage and not high_novelty:
        quadrant = 'III'
    else:
        quadrant = 'IV'

    return {
        'quadrant': quadrant,
        'quadrant_label': QUADRANT_LABELS[quadrant],
        'is_dark_matter': quadrant == 'I',
    }


def run_step3(
    chain_novelty_path: str,
    output_path: str,
    step1_path: str = None,
    disease: str = "IBD",
) -> list:
    """
    Run step 3: assign quadrants to all candidates.

    Args:
        chain_novelty_path: Path to step2 output JSON
        output_path: Path to save step3 output JSON
        step1_path: Optional path to step1 output JSON (for latest mmsage_norm)

    Returns:
        List of quadrant-assigned candidate dicts
    """
    print("=" * 60)
    print("  STEP 3: Dual-Axis Quadrant Assignment")
    print("=" * 60)

    with open(chain_novelty_path, 'r', encoding='utf-8') as f:
        novelty_results = json.load(f)

    # Build lookup from latest step1 output if available
    step1_lookup = {}
    if step1_path:
        p = Path(step1_path)
        if p.exists():
            with open(p, 'r', encoding='utf-8') as f:
                for c in json.load(f):
                    key = (c.get('bacteria', ''), c.get('metabolite', ''))
                    step1_lookup[key] = c.get('mmsage_norm', 0)
            print(f"  Loaded {len(step1_lookup)} mmsage_norm values from step1")

    # Resolve mmsage_norm: prefer step1 latest, fallback to step2 embedded
    for item in novelty_results:
        cand = item.get('candidate', {})
        key = (item.get('bacteria', cand.get('bacteria', '')),
               item.get('metabolite', cand.get('metabolite', '')))
        if key in step1_lookup:
            cand['mmsage_norm'] = step1_lookup[key]

    # Compute adaptive thresholds using median of each axis
    mmsage_values = sorted(
        item.get('candidate', {}).get('mmsage_norm', 0)
        for item in novelty_results)
    novelty_values = sorted(item.get('chain_novelty', 1.0)
                            for item in novelty_results)
    n = len(mmsage_values)
    mmsage_threshold = (mmsage_values[n // 2] + mmsage_values[(n - 1) // 2]) / 2 if n else 0.5
    novelty_threshold = (novelty_values[n // 2] + novelty_values[(n - 1) // 2]) / 2 if n else 0.5

    print(f"\nAdaptive thresholds (median):")
    print(f"  MMSage: {mmsage_threshold:.4f}  (range: {min(mmsage_values):.4f} - {max(mmsage_values):.4f})")
    print(f"  Chain Novelty: {novelty_threshold:.4f}  (range: {min(novelty_values):.4f} - {max(novelty_values):.4f})")

    results = []
    quadrant_counts = {'I': 0, 'II': 0, 'III': 0, 'IV': 0}

    # Auto-detect disease from step2 data if not explicitly provided
    if disease == "IBD" and novelty_results:
        detected = novelty_results[0].get('disease', '')
        if detected and detected != "IBD":
            disease = detected
            print(f"  Auto-detected disease from step2 data: {disease}")

    for item in novelty_results:
        cand = item.get('candidate', {})
        mmsage_norm = cand.get('mmsage_norm', 0)
        chain_novelty = item.get('chain_novelty', 1.0)
        chain_count = item.get('chain_count', 0)
        pairwise = item.get('pairwise_counts', {})
        # Evidence foundation: how many experimental papers link this metabolite to the disease
        # Higher = more feasible to design validation experiments
        evidence_foundation = pairwise.get('pair_md_exp', 0)
        pair_bm_exp = pairwise.get('pair_bm_exp', 0)
        pair_md_exp = pairwise.get('pair_md_exp', 0)

        quad_info = assign_quadrant(mmsage_norm, chain_novelty,
                                     mmsage_threshold, novelty_threshold)

        results.append({
            # Core fields
            'bacteria': item.get('bacteria', cand.get('bacteria', '')),
            'metabolite': item.get('metabolite', cand.get('metabolite', '')),
            'disease': disease,
            'mmsage_norm': mmsage_norm,
            'chain_count': chain_count,
            'chain_novelty': chain_novelty,
            'evidence_foundation': evidence_foundation,
            'pair_bm_exp': pair_bm_exp,
            'pair_md_exp': pair_md_exp,
            # Quadrant
            **quad_info,
            # Chain details
            'chain_path_str': item.get('chain_path_str', ''),
            'chain_path': item.get('chain_path', []),
            'has_path': item.get('has_path', False),
            'chain_query': item.get('chain_query', ''),
            'edge_cooccurrences': item.get('edge_cooccurrences', []),
            'bottleneck_edge': item.get('bottleneck_edge'),
            # Original candidate data
            'candidate': cand,
        })

        quadrant_counts[quad_info['quadrant']] += 1

    # Print summary
    print(f"\nQuadrant distribution:")
    for q, label in QUADRANT_LABELS.items():
        count = quadrant_counts[q]
        marker = ' ★' if q == 'I' else ''
        print(f"  {q} ({label}): {count}{marker}")

    # Sort by quadrant priority (I > II > III > IV) then composite_score descending
    _quadrant_priority = {'I': 0, 'II': 1, 'III': 2, 'IV': 3}
    for r in results:
        bm = r.get('pair_bm_exp', 0)
        md = r.get('pair_md_exp', 0)
        mms = r.get('mmsage_norm', 0)
        r['composite_score'] = round(
            mms * (1 + math.log2(1 + bm)) * (1 + math.log2(1 + md)), 4)
    results.sort(key=lambda r: (
        _quadrant_priority.get(r['quadrant'], 9),
        -r['composite_score'],
    ))

    # Highlight dark matter candidates
    dark_matter = [r for r in results if r['is_dark_matter']]
    if dark_matter:
        print(f"\nDark matter candidates (Quadrant I, {len(dark_matter)}):")
        for dm in dark_matter[:10]:
            print(f"  {dm['bacteria']} x {dm['metabolite']} "
                  f"(MMSage={dm['mmsage_norm']:.3f}, "
                  f"chain_count={dm['chain_count']}, "
                  f"ChainNov={dm['chain_novelty']:.3f})")
        if len(dark_matter) > 10:
            print(f"  ... and {len(dark_matter) - 10} more")

    # Known relationships (Quadrant III)
    known = [r for r in results if r['quadrant'] == 'III']
    if known:
        print(f"\nKnown relationships (Quadrant III, {len(known)}):")
        for k in known[:5]:
            print(f"  {k['bacteria']} x {k['metabolite']} "
                  f"(MMSage={k['mmsage_norm']:.3f}, "
                  f"chain_count={k['chain_count']}, "
                  f"ChainNov={k['chain_novelty']:.3f})")

    # Save
    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    print(f"\nSaved {len(results)} quadrant results to {output_path}")
    return results


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Step 3: Quadrant Assignment')
    parser.add_argument('--input', default='outputs/step2_chain_novelty.json')
    parser.add_argument('--output', default='outputs/step3_quadrant.json')
    args = parser.parse_args()

    run_step3(args.input, args.output)
