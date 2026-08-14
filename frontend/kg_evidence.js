/**
 * Knowledge graph node click → literature evidence panel (shared by Dashboard & browse/record).
 */
(function (global) {
  const HOP_LABELS = {
    microbe_metabolite: "Microbe–Metabolite",
    metabolite_disease: "Metabolite–Disease",
    microbe_disease: "Microbe–Disease",
  };

  let _boundCy = null;
  let _getContext = () => ({});

  function escapeHtml(s) {
    if (s == null) return "";
    return String(s)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;");
  }

  function formatKgLabel(s) {
    const acronyms = new Set(["ibd", "ibs", "crc", "uc", "cd", "nafld", "t2d"]);
    return String(s || "")
      .replace(/_/g, " ")
      .split(/\s+/)
      .filter(Boolean)
      .map((word) => {
        const lower = word.toLowerCase();
        if (acronyms.has(lower)) return lower.toUpperCase();
        return word.charAt(0).toUpperCase() + word.slice(1);
      })
      .join(" ");
  }

  function pmidLink(pmid) {
    const id = escapeHtml(pmid);
    return `<a href="https://pubmed.ncbi.nlm.nih.gov/${id}/" target="_blank" rel="noopener" class="text-indigo-600 hover:underline font-mono text-xs">PMID ${id}</a>`;
  }

  function renderEdgeEvidence(edges) {
    if (!edges || !edges.length) {
      return '<p class="text-gray-400 text-xs">No KG edge literature for this node.</p>';
    }
    return edges
      .map((e) => {
        const arrow = e.direction === "outgoing" ? "→" : "←";
        const rel = escapeHtml(e.relation || "related_to");
        const desc = e.description ? `<p class="text-gray-600 mt-1">${escapeHtml(e.description)}</p>` : "";
        const pmids = (e.pmids || []).map(pmidLink).join(" · ");
        const meta = [];
        if (e.impact_factor) meta.push(`IF ${e.impact_factor}`);
        if (e.citation_count) meta.push(`${e.citation_count} citations`);
        const metaHtml = meta.length
          ? `<p class="text-gray-400 text-xs mt-1">${meta.join(" · ")}</p>`
          : "";
        return `<div class="border border-gray-100 rounded-lg p-3 mb-2 bg-gray-50/80">
          <p class="font-medium text-gray-800 text-xs">${escapeHtml(formatKgLabel(e.neighbor_label))} <span class="text-gray-400">${arrow}</span></p>
          <p class="text-xs text-gray-500 mt-0.5">${rel}</p>
          ${desc}
          ${pmids ? `<p class="mt-2 flex flex-wrap gap-2">${pmids}</p>` : '<p class="text-gray-400 text-xs mt-1">No PMID on this edge</p>'}
          ${metaHtml}
        </div>`;
      })
      .join("");
  }

  function renderHopEvidence(hops) {
    if (!hops || !hops.length) return "";
    return `<div class="mb-4">
      <h4 class="text-xs font-semibold text-gray-500 uppercase tracking-wide mb-2">PubMed co-occurrence (multi-agent)</h4>
      ${hops
        .map((h) => {
          const title = HOP_LABELS[h.hop_type] || h.hop_type;
          const sources = (h.sources || [])
            .map((s) => `<li class="text-gray-700">${escapeHtml(s)}</li>`)
            .join("");
          const dbHits = h.db_hits
            ? Object.entries(h.db_hits)
                .filter(([, v]) => v)
                .map(([k]) => k)
                .join(", ")
            : "";
          return `<div class="border-l-4 border-indigo-400 pl-3 mb-3">
            <p class="font-medium text-gray-800 text-sm">${escapeHtml(title)}</p>
            <p class="text-xs text-gray-500 mt-0.5">PubMed hits: <span class="font-mono">${h.pubmed_count ?? 0}</span>${dbHits ? ` · DB: ${escapeHtml(dbHits)}` : ""}</p>
            ${h.query_used ? `<p class="text-xs text-gray-400 mt-1 break-all">Query: ${escapeHtml(h.query_used)}</p>` : ""}
            ${sources ? `<ul class="list-disc list-inside text-xs mt-2 space-y-1">${sources}</ul>` : ""}
          </div>`;
        })
        .join("")}
    </div>`;
  }

  function renderArticles(articles) {
    if (!articles || !articles.length) return "";
    return `<div>
      <h4 class="text-xs font-semibold text-gray-500 uppercase tracking-wide mb-2">Linked abstracts</h4>
      ${articles
        .map((a) => {
          const title = escapeHtml(a.title || "(no title)");
          const journal = [a.journal, a.year].filter(Boolean).join(", ");
          const abs = a.abstract
            ? `<p class="text-gray-600 text-xs mt-1 line-clamp-3">${escapeHtml(a.abstract)}</p>`
            : "";
          return `<div class="mb-3 pb-3 border-b border-gray-100 last:border-0">
            <p class="text-sm font-medium text-gray-900">${title}</p>
            <p class="text-xs mt-1">${pmidLink(a.pmid)}${journal ? ` · <span class="text-gray-500">${escapeHtml(journal)}</span>` : ""}</p>
            ${abs}
          </div>`;
        })
        .join("")}
    </div>`;
  }

  function renderPayload(data) {
    const node = data.node || {};
    const typeLabel = escapeHtml(node.type || "?");
    const role = node.core_role ? `<span class="text-emerald-700 font-medium"> · ${escapeHtml(node.core_role)}</span>` : "";
    let html = `<div class="mb-3">
      <p class="text-base font-semibold text-gray-900">${escapeHtml(formatKgLabel(node.label || node.id))}</p>
      <p class="text-xs text-gray-500">${typeLabel}${role}</p>`;
    if (node.description) {
      html += `<p class="text-sm text-gray-600 mt-2">${escapeHtml(node.description)}</p>`;
    }
    if (node.kegg_id) {
      html += `<p class="mt-1"><a href="https://www.kegg.jp/entry/${escapeHtml(node.kegg_id)}" target="_blank" class="text-indigo-600 hover:underline text-xs">KEGG ${escapeHtml(node.kegg_id)}</a></p>`;
    }
    html += `</div>`;

    html += renderHopEvidence(data.hop_evidence);
    html += `<div class="mb-4">
      <h4 class="text-xs font-semibold text-gray-500 uppercase tracking-wide mb-2">KG edge literature (${(data.edge_evidence || []).length})</h4>
      ${renderEdgeEvidence(data.edge_evidence)}
    </div>`;
    html += renderArticles(data.articles);
    return html;
  }

  async function showNode(nodeId, nodeData) {
    const panel = document.getElementById("kg-evidence-panel");
    const content = document.getElementById("kg-evidence-content");
    if (!panel || !content) return;

    panel.classList.remove("hidden");
    content.innerHTML = '<p class="text-gray-400 text-sm py-2">Loading evidence...</p>';

    const ctx = _getContext() || {};
    const params = new URLSearchParams({
      node_id: nodeId,
      bacteria: ctx.bacteria || "",
      metabolite: ctx.metabolite || "",
      disease: ctx.disease || "IBD",
      node_type: nodeData.type || "",
      core_role: nodeData.core_role || "",
    });
    if (ctx.run_id) params.set("run_id", ctx.run_id);

    try {
      const base = ctx.api_base || "";
      const res = await fetch(`${base}/api/graph/node-evidence?${params}`);
      const data = await res.json();
      if (!res.ok || data.error) throw new Error(data.error || res.statusText);
      content.innerHTML = renderPayload(data);
    } catch (err) {
      content.innerHTML = `<p class="text-red-600 text-sm">${escapeHtml(err.message || "Failed to load evidence")}</p>`;
    }
  }

  function onNodeTap(evt) {
    const node = evt.target;
    if (!node.isNode()) return;
    showNode(node.id(), node.data());
  }

  function bind(cy, getContext) {
    if (!cy) return;
    if (_boundCy && _boundCy !== cy) {
      _boundCy.off("tap", "node", onNodeTap);
    }
    _boundCy = cy;
    _getContext = typeof getContext === "function" ? getContext : () => ({});
    cy.off("tap", "node", onNodeTap);
    cy.on("tap", "node", onNodeTap);
  }

  function clearPanel() {
    const panel = document.getElementById("kg-evidence-panel");
    const content = document.getElementById("kg-evidence-content");
    if (panel) panel.classList.add("hidden");
    if (content) content.innerHTML = "";
  }

  global.KgEvidence = { bind, showNode, clearPanel, renderPayload };
})(typeof window !== "undefined" ? window : global);
