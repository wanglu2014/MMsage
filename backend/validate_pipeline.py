"""
Quick validation: manually inject known pairs into candidates to verify
all four quadrants work correctly with real PubMed queries.
"""
import json
import sys
import io
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
if sys.stdout.encoding and sys.stdout.encoding.lower() != 'utf-8':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

from step2_chain_novelty import run_step2
from step3_quadrant import run_step3

# Create synthetic candidates that cover all scenarios
test_candidates = [
    # Should find KG path: Akk -> Propionate (known, high MMSage)
    {"bacteria": "Akkermansia_muciniphila", "metabolite": "Propionate",
     "pseudotime": 55.0, "mmsage_norm": 0.90, "rank_in_microbe": 1,
     "total_in_microbe": 100, "source_file": "test"},
    # Should find KG path: Faecalibacterium -> Butyrate (known, high MMSage)
    {"bacteria": "Faecalibacterium_prausnitzii", "metabolite": "Butyrate",
     "pseudotime": 58.0, "mmsage_norm": 0.95, "rank_in_microbe": 1,
     "total_in_microbe": 100, "source_file": "test"},
    # No KG path: Akk -> Isobutyric acid (dark matter, high MMSage)
    {"bacteria": "Akkermansia_muciniphila", "metabolite": "Isobutyric acid",
     "pseudotime": 57.0, "mmsage_norm": 0.91, "rank_in_microbe": 2,
     "total_in_microbe": 100, "source_file": "test"},
    # No KG path, low MMSage -> Q-II
    {"bacteria": "Akkermansia_muciniphila", "metabolite": "Betaine",
     "pseudotime": 20.0, "mmsage_norm": 0.30, "rank_in_microbe": 50,
     "total_in_microbe": 100, "source_file": "test"},
    # Should find KG path: E.coli -> Succinate (known, low MMSage) -> Q-IV
    {"bacteria": "Escherichia_coli", "metabolite": "Succinate",
     "pseudotime": 25.0, "mmsage_norm": 0.40, "rank_in_microbe": 30,
     "total_in_microbe": 100, "source_file": "test"},
]

out_dir = Path(__file__).parent.parent / "outputs" / "validation"
out_dir.mkdir(parents=True, exist_ok=True)

# Save test candidates
cand_path = str(out_dir / "test_candidates.json")
with open(cand_path, 'w', encoding='utf-8') as f:
    json.dump(test_candidates, f, indent=2)

gml_path = str(Path(__file__).parent.parent / "data" / "knowledge_graph" / "ibd_test_kg.gml")

# Run step2
step2_out = str(out_dir / "test_chain_novelty.json")
results2 = run_step2(cand_path, gml_path, step2_out, disease="IBD")

# Run step3
step3_out = str(out_dir / "test_quadrant.json")
results3 = run_step3(step2_out, step3_out)

# Summary
print("\n" + "=" * 70)
print("  VALIDATION RESULTS")
print("=" * 70)
print(f"\n{'Bacteria':<30} {'Metabolite':<20} {'MMSage':>7} {'ChainCt':>8} {'ChainNov':>9} {'Quad':>5} {'HasPath':>8}")
print("-" * 95)
for r in results3:
    print(f"{r['bacteria']:<30} {r['metabolite']:<20} {r['mmsage_norm']:>7.3f} "
          f"{r['chain_count']:>8} {r['chain_novelty']:>9.3f} {r['quadrant']:>5} "
          f"{'Yes' if r['has_path'] else 'No':>8}")

# Check expected quadrants
print("\n=== Expected vs Actual ===")
expected = {
    "Propionate": ("III", "Known: Akk produces Propionate, high PubMed co-occurrence"),
    "Butyrate": ("III", "Known: Faecal produces Butyrate, high PubMed co-occurrence"),
    "Isobutyric acid": ("I", "Dark matter: no KG path Akk->IBA"),
    "Betaine": ("II", "No path + low MMSage"),
    "Succinate": ("IV or III", "E.coli->Succinate path exists, low MMSage"),
}
for r in results3:
    met = r['metabolite']
    exp_q, reason = expected.get(met, ("?", ""))
    actual_q = r['quadrant']
    match = "OK" if exp_q in actual_q else "CHECK"
    print(f"  {met:<20} Expected={exp_q:<8} Actual={actual_q:<5} [{match}] {reason}")
