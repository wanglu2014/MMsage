const {
  useEffect,
  useState,
  useMemo,
  useCallback,
  useRef
} = React;
const API_BASE = window.location.origin;
const MM_TO_PX_150 = 150 / 25.4;

/** Step 1: full A4 portrait width, height = 1/3 page (210 × 99 mm) */
const A4_PORTRAIT_THIRD = {
  widthMm: 210,
  heightMm: 297 / 3,
  widthPx: Math.round(210 * MM_TO_PX_150),
  heightPx: Math.round(297 / 3 * MM_TO_PX_150)
};

/** A4 landscape — 297 × 210 mm @ ~150 dpi (Step 2 KG export/print) */
const A4_LANDSCAPE = {
  widthMm: 297,
  heightMm: 210,
  widthPx: 1754,
  heightPx: 1240
};

/** Publication figure settings */
const FIGURE = {
  plot: {
    format: "svg",
    width: A4_PORTRAIT_THIRD.widthPx,
    height: A4_PORTRAIT_THIRD.heightPx
  },
  plotPng: {
    format: "png",
    width: A4_PORTRAIT_THIRD.widthPx,
    height: A4_PORTRAIT_THIRD.heightPx,
    scale: 3
  },
  kgSvg: {
    full: true,
    scale: 2,
    bg: "#ffffff"
  },
  kgPngScale: 4
};
function normalizeSvgSize(svg, widthMm, heightMm) {
  if (!svg || typeof svg !== "string") return svg;
  return svg.replace(/<svg([^>]*)>/i, (match, attrs) => {
    const cleaned = attrs.replace(/\s*width\s*=\s*["'][^"']*["']/gi, "").replace(/\s*height\s*=\s*["'][^"']*["']/gi, "").replace(/\s*preserveAspectRatio\s*=\s*["'][^"']*["']/gi, "");
    return `<svg${cleaned} width="${widthMm}mm" height="${heightMm}mm" preserveAspectRatio="xMidYMid meet">`;
  });
}
function normalizeSvgA4PortraitThird(svg) {
  return normalizeSvgSize(svg, A4_PORTRAIT_THIRD.widthMm, A4_PORTRAIT_THIRD.heightMm);
}
function normalizeSvgA4Landscape(svg) {
  return normalizeSvgSize(svg, A4_LANDSCAPE.widthMm, A4_LANDSCAPE.heightMm);
}
function plotlyToPrintImage(plotEl) {
  if (!plotEl || !window.Plotly?.toImage) return Promise.resolve(null);
  return window.Plotly.toImage(plotEl, FIGURE.plot).then(dataUri => {
    const svg = dataUriToSvg(dataUri);
    return svg ? svgToDataUri(normalizeSvgA4PortraitThird(svg)) : dataUri;
  }).catch(() => null);
}
function cytoscapeSvgString(cy) {
  if (!cy || typeof cy.svg !== "function") return null;
  try {
    cy.resize();
    cy.fit(undefined, 30);
    return normalizeSvgA4Landscape(cy.svg(FIGURE.kgSvg));
  } catch (e) {
    console.warn("KG SVG export failed", e);
    return null;
  }
}
function svgToDataUri(svg) {
  return "data:image/svg+xml;charset=utf-8," + encodeURIComponent(svg);
}
function dataUriToSvg(dataUri) {
  if (!dataUri || !dataUri.startsWith("data:")) return null;
  const comma = dataUri.indexOf(",");
  if (comma < 0) return null;
  const meta = dataUri.slice(0, comma);
  const payload = dataUri.slice(comma + 1);
  if (meta.includes(";base64")) {
    try {
      return atob(payload);
    } catch (e) {
      return null;
    }
  }
  return payload.startsWith("<") ? payload : decodeURIComponent(payload);
}
function downloadPlotSvg(plotEl) {
  if (!plotEl || !window.Plotly?.toImage) return;
  plotlyToPrintImage(plotEl).then(dataUri => {
    const svg = dataUriToSvg(dataUri);
    if (!svg) return;
    const blob = new Blob([svg], {
      type: "image/svg+xml;charset=utf-8"
    });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = "quadrant_plot.svg";
    a.click();
    URL.revokeObjectURL(a.href);
  }).catch(() => null);
}
function cytoscapeToPrintSvgUri(cy) {
  const svg = cytoscapeSvgString(cy);
  if (svg) return svgToDataUri(svg);
  try {
    return cy.png({
      output: "base64uri",
      bg: "#ffffff",
      scale: FIGURE.kgPngScale,
      full: true
    });
  } catch (e) {
    return null;
  }
}
function downloadKgSvg(cy) {
  const svg = cytoscapeSvgString(cy);
  if (!svg) return;
  const blob = new Blob([svg], {
    type: "image/svg+xml;charset=utf-8"
  });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = "kg_chain_view.svg";
  a.click();
  URL.revokeObjectURL(a.href);
}
const QUADRANT_ORDER = ["I", "II", "III", "IV"];
const QUADRANT_PRIORITY = {
  I: 0,
  II: 1,
  III: 2,
  IV: 3
};
const QUADRANT_HEADERS = {
  I: "I · Dark Matter ★",
  II: "II · Novel, Weak Signal",
  III: "III · Known Relationship",
  IV: "IV · Low Priority"
};
function sortCandidatesGrouped(candidates) {
  return [...candidates].sort((a, b) => {
    const pa = QUADRANT_PRIORITY[a.quadrant] ?? 9;
    const pb = QUADRANT_PRIORITY[b.quadrant] ?? 9;
    if (pa !== pb) return pa - pb;
    return (b.mmsage_norm ?? 0) - (a.mmsage_norm ?? 0);
  });
}
function buildCandidateGroups(candidates) {
  const sorted = sortCandidatesGrouped(candidates);
  const groups = [];
  for (const quadrant of QUADRANT_ORDER) {
    const items = sorted.filter(c => c.quadrant === quadrant);
    if (!items.length) continue;
    groups.push({
      quadrant,
      items
    });
  }
  const other = sorted.filter(c => !QUADRANT_ORDER.includes(c.quadrant));
  if (other.length) {
    groups.push({
      quadrant: "?",
      items: other
    });
  }
  return groups;
}
function formatLabel(value) {
  const acronyms = new Set(["ibd", "ibs", "crc", "uc", "cd", "nafld", "t2d"]);
  return String(value || "").replace(/_/g, " ").split(/\s+/).filter(Boolean).map(word => {
    const lower = word.toLowerCase();
    if (acronyms.has(lower)) return lower.toUpperCase();
    return word.charAt(0).toUpperCase() + word.slice(1);
  }).join(" ");
}
function supportLevelClass(level) {
  switch ((level || "").toLowerCase()) {
    case "direct":
      return "bg-emerald-100 text-emerald-800";
    case "indirect":
      return "bg-amber-100 text-amber-800";
    default:
      return "bg-gray-100 text-gray-700";
  }
}
function App() {
  const initialRunId = new URLSearchParams(window.location.search).get("run_id") || "";
  const [candidates, setCandidates] = useState([]);
  const [thresholds, setThresholds] = useState({
    mmsage: 0.5,
    evidence: 0
  });
  const [selectedPair, setSelectedPair] = useState(null);
  const [analysisResult, setAnalysisResult] = useState(null);
  const [protocol, setProtocol] = useState("");
  const [protocolError, setProtocolError] = useState("");
  const [validationPlan, setValidationPlan] = useState(null);
  const [validationPlanError, setValidationPlanError] = useState("");
  const [validationPlanLoading, setValidationPlanLoading] = useState(false);
  const [validationPlanJobId, setValidationPlanJobId] = useState("");
  const [validationPlanStatus, setValidationPlanStatus] = useState(null);
  const [questionPlan, setQuestionPlan] = useState(null);
  const [questionPlanError, setQuestionPlanError] = useState("");
  const [questionPlanLoading, setQuestionPlanLoading] = useState(false);
  const [questionPlanJobId, setQuestionPlanJobId] = useState("");
  const [questionPlanStatus, setQuestionPlanStatus] = useState(null);
  const [questionPlanResearchQuestion, setQuestionPlanResearchQuestion] = useState("");
  const [questionPlanPromptConstraints, setQuestionPlanPromptConstraints] = useState("");
  const [loading, setLoading] = useState(false);
  const [uploading, setUploading] = useState(false);
  const [uploadMsg, setUploadMsg] = useState("");
  const [disease, setDisease] = useState("IBD");
  const [maxQueries, setMaxQueries] = useState(15);
  const [maxArticlesPerQuery, setMaxArticlesPerQuery] = useState(20);
  const [currentRunId, setCurrentRunId] = useState(initialRunId);
  const [defaultRunId, setDefaultRunId] = useState("");
  const [currentJobId, setCurrentJobId] = useState("");
  const [plotPrintImg, setPlotPrintImg] = useState("");
  const [kgPrintImg, setKgPrintImg] = useState("");
  const [unchartedExpanded, setUnchartedExpanded] = useState(false);
  const fileRef = useRef(null);
  const UNCHARTED_PREVIEW = 3;
  const candidateGroups = useMemo(() => buildCandidateGroups(candidates), [candidates]);
  const unchartedDarkMatter = useMemo(() => candidates.filter(c => (c.pair_bm_exp ?? 0) + (c.pair_md_exp ?? 0) === 0).sort((a, b) => (b.mmsage_norm ?? 0) - (a.mmsage_norm ?? 0)), [candidates]);
  const buildMechanismSummary = useCallback(pair => {
    if (!pair) return "";
    return `Composite: ${(pair.composite_score ?? 0).toFixed(2)}, MMSage: ${(pair.mmsage_norm ?? 0).toFixed(3)}, bm_exp: ${pair.pair_bm_exp ?? 0}, md_exp: ${pair.pair_md_exp ?? 0}, Novelty: ${(pair.novelty_score ?? 0).toFixed(3)}`;
  }, []);
  const updateDashboardUrl = useCallback(runId => {
    const url = new URL(window.location.href);
    if (runId) {
      url.searchParams.set("run_id", runId);
    } else {
      url.searchParams.delete("run_id");
    }
    window.history.replaceState({}, "", `${url.pathname}${url.search}`);
  }, []);
  const formatValidationStage = useCallback(stage => {
    const labels = {
      queued: "Queued",
      started: "Worker Started",
      prepare_support: "Load Candidate Support",
      round1_search: "Round 1 Literature Search",
      audit_protocol: "Protocol Audit",
      audit_scope: "Question Scope Audit",
      followup_queries: "Generate Follow-up Queries",
      round2_search: "Follow-up Literature Search",
      fulltext_methods: "Extract Full-Text Methods",
      protocols_io_search: "Search protocols.io",
      evidence_adjudication: "Adjudicate Direct Evidence & Hypotheses",
      in_vitro_design: "Design In Vitro Package",
      in_vivo_design: "Design In Vivo Package",
      self_reflection: "Self-Reflection & Revision",
      prepare_question: "Build Question Focus",
      assemble: "Assemble Final Plan",
      completed: "Completed",
      error: "Error"
    };
    return labels[stage] || formatLabel(stage || "running");
  }, []);
  const loadCandidates = useCallback(async (targetRunId = "") => {
    const params = new URLSearchParams();
    if (targetRunId) params.set("run_id", targetRunId);
    const qs = params.toString();
    const res = await fetch(`${API_BASE}/api/dual-axis-candidates${qs ? `?${qs}` : ""}`);
    const data = await res.json();
    setCandidates(data?.candidates || []);
    if (data?.thresholds) setThresholds(data.thresholds);
    setSelectedPair(null);
    setAnalysisResult(null);
    const nextRunId = data?.run_id || targetRunId || "";
    setCurrentRunId(prevRunId => {
      if (prevRunId !== nextRunId) {
        setProtocol("");
        setProtocolError("");
        setValidationPlan(null);
        setValidationPlanError("");
        setValidationPlanLoading(false);
        setValidationPlanJobId("");
        setValidationPlanStatus(null);
      }
      return nextRunId;
    });
    setDefaultRunId(data?.default_dashboard_run_id || "");
    setDisease(data?.disease || "IBD");
  }, []);

  // Load candidates
  useEffect(() => {
    loadCandidates(initialRunId).catch(console.error);
  }, [initialRunId, loadCandidates]);

  // Render quadrant plot
  useEffect(() => {
    if (candidates.length > 0) {
      setTimeout(() => renderQuadrantPlot(candidates, thresholds), 100);
    }
  }, [candidates, thresholds]);
  useEffect(() => {
    setUnchartedExpanded(false);
  }, [candidates]);
  useEffect(() => {
    setProtocol("");
    setProtocolError("");
    setValidationPlan(null);
    setValidationPlanError("");
    setValidationPlanLoading(false);
    setValidationPlanJobId("");
    setValidationPlanStatus(null);
  }, [selectedPair?.bacteria, selectedPair?.metabolite, disease]);
  useEffect(() => {
    if (!validationPlanJobId || !validationPlanLoading) return undefined;
    let cancelled = false;
    const pollStatus = async () => {
      try {
        const res = await fetch(`${API_BASE}/api/validation-plan/status/${encodeURIComponent(validationPlanJobId)}`);
        const data = await res.json().catch(() => ({}));
        if (cancelled) return;
        if (!res.ok) {
          throw new Error(parseApiError(data, `Validation plan status failed (${res.status})`));
        }
        setValidationPlanStatus(data);
        if (data.status === "completed") {
          setValidationPlan(data.result || null);
          setProtocol(data?.result?.validation_protocol_text ? String(data.result.validation_protocol_text).trim() : "");
          setProtocolError("");
          setLoading(false);
          setValidationPlanLoading(false);
          return;
        }
        if (data.status === "error") {
          setValidationPlanError(data.error || data.current_message || "Validation plan generation failed.");
          setProtocol("");
          setProtocolError(data.error || data.current_message || "Validation plan generation failed.");
          setLoading(false);
          setValidationPlanLoading(false);
          return;
        }
        window.setTimeout(pollStatus, 1800);
      } catch (err) {
        if (cancelled) return;
        setValidationPlanError(err.message || "Validation plan status polling failed.");
        setProtocol("");
        setProtocolError(err.message || "Validation plan status polling failed.");
        setLoading(false);
        setValidationPlanLoading(false);
      }
    };
    pollStatus();
    return () => {
      cancelled = true;
    };
  }, [validationPlanJobId, validationPlanLoading]);
  useEffect(() => {
    if (!questionPlanJobId || !questionPlanLoading) return undefined;
    let cancelled = false;
    const pollStatus = async () => {
      try {
        const res = await fetch(`${API_BASE}/api/validation-plan/status/${encodeURIComponent(questionPlanJobId)}`);
        const data = await res.json().catch(() => ({}));
        if (cancelled) return;
        if (!res.ok) {
          throw new Error(parseApiError(data, `Standalone validation status failed (${res.status})`));
        }
        setQuestionPlanStatus(data);
        if (data.status === "completed") {
          setQuestionPlan(data.result || null);
          setQuestionPlanLoading(false);
          return;
        }
        if (data.status === "error") {
          setQuestionPlanError(data.error || data.current_message || "Standalone validation generation failed.");
          setQuestionPlanLoading(false);
          return;
        }
        window.setTimeout(pollStatus, 1800);
      } catch (err) {
        if (cancelled) return;
        setQuestionPlanError(err.message || "Standalone validation status polling failed.");
        setQuestionPlanLoading(false);
      }
    };
    pollStatus();
    return () => {
      cancelled = true;
    };
  }, [questionPlanJobId, questionPlanLoading]);

  // Reflow charts before printing to avoid canvas/SVG offset.
  useEffect(() => {
    const handleBeforePrint = () => {
      const plotEl = document.getElementById("main-plot");
      if (plotEl && window.Plotly) {
        plotlyToPrintImage(plotEl).then(url => {
          if (url) setPlotPrintImg(url);
        });
      }
      if (window._kgCy) {
        window._kgCy.resize();
        window._kgCy.fit(undefined, 30);
        const dataUri = cytoscapeToPrintSvgUri(window._kgCy);
        if (dataUri) setKgPrintImg(dataUri);
      }
    };
    window.addEventListener("beforeprint", handleBeforePrint);
    // Some browsers fire print layout after beforeprint; run once more shortly after.
    const mediaQuery = window.matchMedia("print");
    const onMqChange = e => {
      if (e.matches) {
        setTimeout(handleBeforePrint, 80);
      }
    };
    if (mediaQuery.addEventListener) {
      mediaQuery.addEventListener("change", onMqChange);
    } else if (mediaQuery.addListener) {
      mediaQuery.addListener(onMqChange);
    }
    return () => {
      window.removeEventListener("beforeprint", handleBeforePrint);
      if (mediaQuery.removeEventListener) {
        mediaQuery.removeEventListener("change", onMqChange);
      } else if (mediaQuery.removeListener) {
        mediaQuery.removeListener(onMqChange);
      }
    };
  }, []);
  const renderQuadrantPlot = useCallback((data, thresh) => {
    const el = document.getElementById("main-plot");
    if (!el) return;

    // Only show points with evidence > 0 (bm+md=0 → see Uncharted Dark Matter panel below)
    const withEvidence = data.filter(d => (d.pair_bm_exp ?? 0) + (d.pair_md_exp ?? 0) > 0);

    // Three groups:
    // 1. Known (novelty < 1.0) → cross marker
    // 2. Top 5% novel by composite_score → star marker + auto-label metabolite name
    // 3. Other novel → circle marker
    const known = withEvidence.filter(d => (d.novelty_score ?? 1) < 1.0);
    const allNovel = withEvidence.filter(d => (d.novelty_score ?? 1) >= 1.0);
    const topN = 10;
    const novelSorted = [...allNovel].sort((a, b) => (b.composite_score ?? 0) - (a.composite_score ?? 0));
    const topNovel = novelSorted.slice(0, topN);
    const topSet = new Set(topNovel.map(d => d.bacteria + "|" + d.metabolite));
    const otherNovel = allNovel.filter(d => !topSet.has(d.bacteria + "|" + d.metabolite));
    const makeTip = d => {
      const bact = d.bacteria.replace(/_/g, " ");
      let tip = `<b>${d.metabolite}</b><br>${bact}`;
      tip += `<br>MMSage: ${d.mmsage_norm.toFixed(3)}`;
      tip += `<br>bm_exp: ${d.pair_bm_exp ?? 0}  |  md_exp: ${d.pair_md_exp ?? 0}`;
      tip += `<br>Composite: ${(d.composite_score ?? 0).toFixed(2)}`;
      tip += `<br>Novelty: ${(d.novelty_score ?? 0).toFixed(3)} ${(d.novelty_score ?? 1) >= 1.0 ? "★ novel" : ""}`;
      return tip;
    };

    // Point size based on composite_score (continuous, good spread)
    // novelty is almost always 1.0 so useless for sizing
    const compScores = withEvidence.map(d => d.composite_score ?? 0);
    const compMax = Math.max(...compScores, 1);
    const compSize = (d, minS, maxS) => {
      const ratio = (d.composite_score ?? 0) / compMax;
      return minS + ratio * (maxS - minS);
    };
    const traces = [
    // Known → cross (×)
    {
      x: known.map(d => d.mmsage_norm),
      y: known.map(d => (d.pair_bm_exp ?? 0) + (d.pair_md_exp ?? 0)),
      mode: "markers",
      type: "scatter",
      name: "× Known",
      marker: {
        color: "#94a3b8",
        size: known.map(d => compSize(d, 7, 16)),
        opacity: 0.6,
        symbol: "x",
        line: {
          width: 2,
          color: "#94a3b8"
        }
      },
      text: known.map(makeTip),
      hovertemplate: "%{text}<extra></extra>",
      _items: known
    },
    // Other novel → circle (●)
    {
      x: otherNovel.map(d => d.mmsage_norm),
      y: otherNovel.map(d => (d.pair_bm_exp ?? 0) + (d.pair_md_exp ?? 0)),
      mode: "markers",
      type: "scatter",
      name: "● Novel",
      marker: {
        color: "#6366f1",
        size: otherNovel.map(d => compSize(d, 7, 18)),
        opacity: 0.7,
        symbol: "circle",
        line: {
          width: 1,
          color: "#ffffff"
        }
      },
      text: otherNovel.map(makeTip),
      hovertemplate: "%{text}<extra></extra>",
      _items: otherNovel
    },
    // Top 5% novel → star (★) with text labels
    {
      x: topNovel.map(d => d.mmsage_norm),
      y: topNovel.map(d => (d.pair_bm_exp ?? 0) + (d.pair_md_exp ?? 0)),
      mode: "markers+text",
      type: "scatter",
      name: "★ Top 10 Novel",
      marker: {
        color: "#f59e0b",
        size: topNovel.map(d => compSize(d, 16, 24)),
        opacity: 1,
        symbol: "star",
        line: {
          width: 1.5,
          color: "#ffffff"
        }
      },
      text: topNovel.map(d => `<b>${d.metabolite.replace(/_/g, " ")}</b>`),
      textposition: "top center",
      textfont: {
        size: 17,
        color: "#78350f",
        family: "Inter, system-ui, sans-serif"
      },
      cliponaxis: false,
      hovertext: topNovel.map(makeTip),
      hovertemplate: "%{hovertext}<extra></extra>",
      _items: topNovel
    }];

    // Size legend: dummy traces to show composite_score → size mapping
    const sizeLegendValues = [0.2, 0.5, 1.0];
    const sizeLegendLabels = ["Composite Low", "Composite Mid", "Composite High"];
    sizeLegendValues.forEach((ratio, i) => {
      const sz = 7 + ratio * 11;
      traces.push({
        x: [null],
        y: [null],
        mode: "markers",
        type: "scatter",
        name: sizeLegendLabels[i],
        marker: {
          color: "#d1d5db",
          size: sz,
          symbol: "circle",
          line: {
            width: 1,
            color: "#9ca3af"
          }
        },
        showlegend: true,
        hoverinfo: "skip",
        legendgroup: "size"
      });
    });
    const allEF = data.map(d => (d.pair_bm_exp ?? 0) + (d.pair_md_exp ?? 0));

    // Use log(1+y) transform for Y axis to spread out the skewed distribution
    // Most candidates have evidence=0, a few have 100+, log scale makes all 4 quadrants visible
    const logY = v => Math.log2(1 + v);
    const yRaw = allEF;
    const yTransformed = yRaw.map(logY);
    const yMaxT = Math.max(...yTransformed, 1) * 1.1;
    const tx = thresh.mmsage;
    // Transform the threshold too; if threshold=0, use a small visual offset so bottom quadrants are visible
    const tyRaw = Math.max(thresh.evidence, 1); // at least 1 paper as threshold if median=0
    const ty = logY(tyRaw);

    // Transform the trace Y values
    traces.forEach(tr => {
      tr.y = tr.y.map(v => logY(v));
    });
    const layout = {
      paper_bgcolor: "#ffffff",
      plot_bgcolor: "#ffffff",
      font: {
        color: "#374151",
        family: "Inter, system-ui, sans-serif",
        size: 15
      },
      margin: {
        t: 62,
        l: 92,
        r: 38,
        b: 88
      },
      xaxis: {
        title: {
          text: "MMSage Signal",
          font: {
            size: 19,
            color: "#4b5563"
          }
        },
        range: [-0.05, 1.05],
        showgrid: false,
        zeroline: false,
        linecolor: "#e5e7eb",
        linewidth: 1,
        tickfont: {
          size: 15,
          color: "#374151"
        },
        ticks: "outside",
        ticklen: 6
      },
      yaxis: {
        title: {
          text: "Evidence Foundation  log₂(1 + bm + md)",
          font: {
            size: 19,
            color: "#4b5563"
          }
        },
        range: [-0.5, yMaxT],
        showgrid: true,
        gridcolor: "#f3f4f6",
        zeroline: false,
        linecolor: "#e5e7eb",
        linewidth: 1,
        tickfont: {
          size: 15,
          color: "#374151"
        },
        ticks: "outside",
        ticklen: 6
      },
      title: {
        text: "MMSage Signal × Evidence Foundation  (★ = novel, no triple evidence)",
        font: {
          size: 21,
          color: "#1f2937"
        }
      },
      legend: {
        x: 0.01,
        y: 0.99,
        bgcolor: "rgba(255,255,255,0.95)",
        bordercolor: "#e5e7eb",
        borderwidth: 1,
        font: {
          size: 14
        }
      },
      shapes: [],
      annotations: []
    };
    Plotly.newPlot("main-plot", traces, layout, {
      responsive: true
    });
    plotlyToPrintImage(el).then(url => {
      if (url) setPlotPrintImg(url);
    });
    el.on("plotly_click", event => {
      const pt = event.points[0];
      const trace = traces[pt.curveNumber];
      if (trace && trace._items) {
        setSelectedPair(trace._items[pt.pointIndex]);
      }
    });
  }, []);

  // Upload coordinates CSV
  const handleUpload = async () => {
    const file = fileRef.current?.files?.[0];
    if (!file) return;
    if (!disease.trim()) {
      setUploadMsg("Please enter a disease name.");
      return;
    }
    if (!Number.isInteger(maxQueries) || maxQueries < 1 || !Number.isInteger(maxArticlesPerQuery) || maxArticlesPerQuery < 1) {
      setUploadMsg("PubMed limits must be positive integers.");
      return;
    }
    setUploading(true);
    setUploadMsg("Uploading...");
    const formData = new FormData();
    formData.append("file", file);
    try {
      const params = new URLSearchParams({
        disease: disease.trim(),
        max_queries: String(maxQueries),
        max_articles_per_query: String(maxArticlesPerQuery)
      });
      const res = await fetch(`${API_BASE}/api/upload-coordinates?${params.toString()}`, {
        method: "POST",
        body: formData
      });
      const data = await res.json();
      if (data.error) {
        setUploadMsg("Error: " + data.error);
        setUploading(false);
        return;
      }
      const jobId = data.job_id || "";
      setCurrentJobId(jobId);
      // Pipeline started in background — poll for progress
      setUploadMsg(`Pipeline started for ${disease}... (step1: MMSage signal)`);
      const poll = setInterval(async () => {
        try {
          const pr = await fetch(`${API_BASE}/api/jobs/${encodeURIComponent(jobId)}/status`);
          const pg = await pr.json();
          if (pg.status === "done") {
            clearInterval(poll);
            setUploadMsg(`Done! ${pg.candidates || 0} candidates for ${pg.disease}. Loading your isolated run...`);
            setUploading(false);
            if (pg.run_id) {
              updateDashboardUrl(pg.run_id);
              await loadCandidates(pg.run_id);
            }
          } else if (pg.status === "error") {
            clearInterval(poll);
            setUploadMsg("Pipeline error: " + (pg.error || "unknown"));
            setUploading(false);
          } else {
            setUploadMsg(`Pipeline running... ${pg.step || ""} (${pg.disease})`);
          }
        } catch (e) {/* ignore poll errors */}
      }, 3000);
    } catch (err) {
      setUploadMsg("Upload failed: " + err.message);
      setUploading(false);
    }
  };
  const parseApiError = (data, fallback) => {
    const d = data?.detail;
    if (Array.isArray(d)) return d.map(x => x.msg || x).join("; ");
    if (typeof d === "string") return d;
    return fallback;
  };

  // Analyze selected pair
  const handleAnalyze = async () => {
    if (!selectedPair) return;
    setLoading(true);
    setProtocol("");
    setProtocolError("");
    setValidationPlan(null);
    setValidationPlanError("");
    setValidationPlanLoading(false);
    setValidationPlanJobId("");
    setValidationPlanStatus(null);
    if (window.KgEvidence) window.KgEvidence.clearPanel();
    const {
      bacteria,
      metabolite
    } = selectedPair;
    const kgParams = new URLSearchParams({
      bacteria,
      metabolite,
      disease
    });
    if (currentRunId) kgParams.set("run_id", currentRunId);
    try {
      const kgRes = await fetch(`${API_BASE}/api/graph/chain?${kgParams.toString()}`);
      const kg = await kgRes.json().catch(() => ({}));
      if (!kgRes.ok) {
        throw new Error(parseApiError(kg, `Knowledge graph request failed (${kgRes.status})`));
      }
      renderKG(kg);
      setAnalysisResult({
        ...selectedPair,
        prioritization: `${bacteria.replace(/_/g, " ")} → ${metabolite} → ${disease} | Composite: ${(selectedPair.composite_score ?? 0).toFixed(2)}`
      });
    } catch (err) {
      console.error("KG render failed:", err);
      setProtocolError(prev => prev ? `${prev} Knowledge graph: ${err.message}` : `Knowledge graph: ${err.message}`);
    }
    try {
      const protoRes = await fetch(`${API_BASE}/api/validation-plan/start`, {
        method: "POST",
        headers: {
          "Content-Type": "application/json"
        },
        body: JSON.stringify({
          bacteria,
          metabolite,
          disease,
          run_id: currentRunId,
          mechanism_summary: buildMechanismSummary(selectedPair),
          protocol_text: "",
          research_question: "",
          prompt_constraints: ""
        })
      });
      const protoData = await protoRes.json().catch(() => ({}));
      if (!protoRes.ok) {
        throw new Error(parseApiError(protoData, `Step 3 validation generation failed (${protoRes.status})`));
      }
      if (!protoData || typeof protoData !== "object" || !protoData.job_id) {
        throw new Error("Server returned an invalid Step 3 job payload.");
      }
      setValidationPlanJobId(protoData.job_id);
      setValidationPlanStatus({
        job_id: protoData.job_id,
        status: protoData.status || "started",
        progress_percent: 0,
        current_stage: "queued",
        current_message: protoData.message || "Step 3 validation protocol job started.",
        logs: []
      });
      setValidationPlanLoading(true);
    } catch (err) {
      console.error("Step 3 validation generation failed:", err);
      setProtocol("");
      setProtocolError(prev => prev ? `${prev} Step 3: ${err.message || "Validation protocol generation failed."}` : err.message || "Validation protocol generation failed.");
      setLoading(false);
    }
  };
  const handleGenerateValidationPlan = async () => {
    if (!questionPlanResearchQuestion.trim()) return;
    setQuestionPlanLoading(true);
    setQuestionPlan(null);
    setQuestionPlanError("");
    setQuestionPlanStatus(null);
    setQuestionPlanJobId("");
    try {
      const planRes = await fetch(`${API_BASE}/api/validation-plan/start`, {
        method: "POST",
        headers: {
          "Content-Type": "application/json"
        },
        body: JSON.stringify({
          bacteria: "",
          metabolite: "",
          disease: "",
          run_id: "",
          mechanism_summary: "",
          protocol_text: "",
          mode: "question_driven",
          research_question: questionPlanResearchQuestion.trim(),
          prompt_constraints: questionPlanPromptConstraints.trim()
        })
      });
      const planData = await planRes.json().catch(() => ({}));
      if (!planRes.ok) {
        throw new Error(parseApiError(planData, `Standalone validation generation failed (${planRes.status})`));
      }
      if (!planData || typeof planData !== "object" || !planData.job_id) {
        throw new Error("Server returned an invalid standalone validation job payload.");
      }
      setQuestionPlanJobId(planData.job_id);
      setQuestionPlanStatus({
        job_id: planData.job_id,
        status: planData.status || "started",
        progress_percent: 0,
        current_stage: "queued",
        current_message: planData.message || "Standalone validation job started.",
        logs: []
      });
    } catch (err) {
      console.error("Standalone validation generation failed:", err);
      setQuestionPlan(null);
      setQuestionPlanError(err.message || "Standalone validation generation failed.");
      setQuestionPlanStatus(null);
      setQuestionPlanJobId("");
      setQuestionPlanLoading(false);
    }
  };
  const renderKG = kg => {
    const container = document.getElementById("kg");
    if (!container) return;
    if (window._kgCy) window._kgCy.destroy();
    const KG_COLORS = {
      microbe: "#3b82f6",
      metabolite: "#f59e0b",
      pathway: "#10b981",
      disease: "#ef4444",
      enzyme: "#f97316",
      receptor: "#06b6d4",
      unknown: "#9ca3af"
    };
    const KG_SHAPES = {
      microbe: "ellipse",
      metabolite: "round-rectangle",
      pathway: "diamond",
      disease: "hexagon",
      enzyme: "triangle",
      receptor: "vee"
    };
    const nodes = (kg.nodes || []).map(node => ({
      ...node,
      data: {
        ...node.data,
        label: formatLabel(node?.data?.label || node?.data?.id || "")
      }
    }));
    const edges = kg.edges || [];
    window._kgCy = cytoscape({
      container,
      elements: [...nodes, ...edges],
      // Enable mouse-wheel zooming and panning.
      userZoomingEnabled: true,
      userPanningEnabled: true,
      style: [
      // Base nodes use larger, high-contrast labels.
      {
        selector: "node",
        style: {
          "label": "data(label)",
          "background-color": function (ele) {
            return KG_COLORS[ele.data("type")] || "#cbd5e1";
          },
          "shape": function (ele) {
            return KG_SHAPES[ele.data("type")] || "ellipse";
          },
          "color": "#111827",
          "font-size": 18,
          // Improve readability in scientific figures.
          "font-weight": "700",
          "text-wrap": "wrap",
          "text-max-width": 180,
          "text-valign": "bottom",
          "text-halign": "center",
          "text-margin-y": 10,
          "text-background-color": "#ffffff",
          "text-background-opacity": 0.96,
          "text-background-padding": 5,
          "text-background-shape": "roundrectangle",
          "text-outline-color": "#ffffff",
          "text-outline-width": 2,
          "width": 66,
          // Enlarge base nodes.
          "height": 66,
          "border-width": 2,
          "border-color": "#ffffff"
        }
      },
      // Emphasize the three core target nodes.
      {
        selector: "node[?is_core]",
        style: {
          "width": 116,
          // Enlarge core nodes.
          "height": 116,
          "font-size": 24,
          // Increase core label size.
          "font-weight": "900",
          "text-max-width": 220,
          "border-width": 5,
          "border-color": "#1e3a8a",
          "text-background-opacity": 1
        }
      },
      // Base edges use thicker strokes and white text outlines.
      {
        selector: "edge",
        style: {
          "width": 3,
          "line-color": "#cbd5e1",
          "target-arrow-color": "#cbd5e1",
          "curve-style": "bezier",
          "target-arrow-shape": "triangle",
          "label": "data(relation)",
          "font-size": 14,
          "color": "#475569",
          "text-rotation": "autorotate",
          "text-outline-color": "#ffffff",
          "text-outline-width": 4
        }
      },
      // Highlight the core inference chain with a thick dark-blue stroke.
      {
        selector: "edge[?is_chain]",
        style: {
          "width": 7,
          "line-color": "#2563eb",
          "target-arrow-color": "#2563eb",
          "font-size": 17,
          "font-weight": "bold",
          "color": "#1e40af",
          "text-outline-width": 4,
          "z-index": 10
        }
      }, {
        selector: "node:selected",
        style: {
          "border-width": 5,
          "border-color": "#eab308"
        }
      }],
      // Tighten layout spacing so the graph fills the viewport.
      layout: {
        name: "dagre",
        rankDir: "LR",
        nodeSep: 30,
        // Reduce vertical spacing.
        rankSep: 90,
        // Reduce rank spacing.
        padding: 15 // Remove excess outer padding.
      }
    });
    const captureKgImage = () => {
      const dataUri = cytoscapeToPrintSvgUri(window._kgCy);
      if (dataUri) setKgPrintImg(dataUri);
    };
    window._kgCy.one("layoutstop", () => setTimeout(captureKgImage, 60));
    setTimeout(captureKgImage, 250);
    window._kgEvidenceContext = () => ({
      bacteria: selectedPair?.bacteria || "",
      metabolite: selectedPair?.metabolite || "",
      disease: disease || "IBD",
      run_id: currentRunId || "",
      api_base: API_BASE
    });
    if (window.KgEvidence) {
      window.KgEvidence.bind(window._kgCy, () => window._kgEvidenceContext());
    }
  };
  return /*#__PURE__*/React.createElement("div", {
    className: "min-h-screen p-6"
  }, /*#__PURE__*/React.createElement("header", {
    className: "mb-6 flex items-center justify-between flex-wrap gap-4"
  }, /*#__PURE__*/React.createElement("div", null, /*#__PURE__*/React.createElement("h1", {
    className: "text-3xl font-bold text-gray-900"
  }, "MMSage-Insight"), /*#__PURE__*/React.createElement("p", {
    className: "text-gray-500"
  }, "MMSage Signal × Evidence Foundation | ★ = Novel (no triple evidence)")), /*#__PURE__*/React.createElement("a", {
    href: "/browse",
    className: "inline-flex items-center gap-2 px-4 py-2 bg-emerald-600 hover:bg-emerald-700 text-white rounded-lg text-sm font-medium transition shadow-sm"
  }, /*#__PURE__*/React.createElement("svg", {
    className: "w-4 h-4",
    fill: "none",
    stroke: "currentColor",
    viewBox: "0 0 24 24",
    "aria-hidden": "true"
  }, /*#__PURE__*/React.createElement("path", {
    strokeLinecap: "round",
    strokeLinejoin: "round",
    strokeWidth: 2,
    d: "M4 7v10c0 2 1 3 3 3h10c2 0 3-1 3-3V7c0-2-1-3-3-3H7C5 4 4 5 4 7zm16 0H4m16 4H4m5 4h6"
  })), "Results Catalog"), currentRunId && defaultRunId && currentRunId !== defaultRunId && /*#__PURE__*/React.createElement("button", {
    onClick: () => {
      updateDashboardUrl("");
      loadCandidates("").catch(console.error);
    },
    className: "inline-flex items-center gap-2 px-4 py-2 bg-white hover:bg-gray-50 text-gray-700 border border-gray-300 rounded-lg text-sm font-medium transition shadow-sm"
  }, "View Home Result")), /*#__PURE__*/React.createElement("div", {
    className: "mb-6 p-4 bg-gray-50 border border-gray-200 rounded-xl flex items-center gap-4 flex-wrap"
  }, /*#__PURE__*/React.createElement("div", null, /*#__PURE__*/React.createElement("span", {
    className: "text-sm font-medium text-gray-700"
  }, "Disease:"), /*#__PURE__*/React.createElement("input", {
    value: disease,
    onChange: e => setDisease(e.target.value),
    className: "ml-2 border border-gray-300 rounded-lg px-3 py-1.5 text-sm w-32",
    placeholder: "e.g. IBD"
  })), /*#__PURE__*/React.createElement("div", null, /*#__PURE__*/React.createElement("span", {
    className: "text-sm font-medium text-gray-700"
  }, "Max PubMed queries:"), /*#__PURE__*/React.createElement("input", {
    type: "number",
    min: "1",
    step: "1",
    value: maxQueries,
    onChange: e => setMaxQueries(Number(e.target.value)),
    className: "ml-2 border border-gray-300 rounded-lg px-3 py-1.5 text-sm w-20"
  })), /*#__PURE__*/React.createElement("div", null, /*#__PURE__*/React.createElement("span", {
    className: "text-sm font-medium text-gray-700"
  }, "Articles per query:"), /*#__PURE__*/React.createElement("input", {
    type: "number",
    min: "1",
    step: "1",
    value: maxArticlesPerQuery,
    onChange: e => setMaxArticlesPerQuery(Number(e.target.value)),
    className: "ml-2 border border-gray-300 rounded-lg px-3 py-1.5 text-sm w-20"
  })), /*#__PURE__*/React.createElement("div", {
    className: "border-l border-gray-300 pl-4"
  }, /*#__PURE__*/React.createElement("span", {
    className: "text-sm font-medium text-gray-700"
  }, "Upload Coordinates CSV:")), /*#__PURE__*/React.createElement("input", {
    ref: fileRef,
    type: "file",
    accept: ".csv",
    className: "text-sm text-gray-600 file:mr-2 file:py-1.5 file:px-4 file:rounded-lg file:border file:border-gray-300 file:bg-white file:text-gray-700 file:cursor-pointer hover:file:bg-gray-50"
  }), /*#__PURE__*/React.createElement("button", {
    onClick: handleUpload,
    disabled: uploading,
    className: "px-5 py-2 bg-indigo-600 hover:bg-indigo-700 disabled:bg-gray-400 text-white rounded-lg text-sm font-medium transition"
  }, uploading ? "Uploading..." : "Upload & Process"), uploadMsg && /*#__PURE__*/React.createElement("span", {
    className: "text-sm text-amber-600"
  }, uploadMsg), currentRunId && /*#__PURE__*/React.createElement("span", {
    className: "text-xs text-gray-500"
  }, "Viewing run: ", currentRunId === defaultRunId ? "home default" : currentRunId), currentJobId && uploading && /*#__PURE__*/React.createElement("span", {
    className: "text-xs text-gray-400 font-mono"
  }, "Job: ", currentJobId)), /*#__PURE__*/React.createElement("div", {
    className: "flex flex-col gap-6"
  }, /*#__PURE__*/React.createElement("div", {
    className: "grid grid-cols-12 gap-6 items-start"
  }, /*#__PURE__*/React.createElement("div", {
    id: "step1-card",
    className: "col-span-12 lg:col-span-8 bg-white rounded-xl p-4 border border-gray-200 shadow-sm flex flex-col"
  }, /*#__PURE__*/React.createElement("h2", {
    className: "text-2xl mb-2 flex items-center gap-2 text-gray-900 font-semibold"
  }, "Step 1: Signal Discovery", /*#__PURE__*/React.createElement("span", {
    className: "text-base text-gray-500 font-normal"
  }, "Click a point to select"), /*#__PURE__*/React.createElement("div", {
    className: "ml-auto flex gap-2 font-normal"
  }, /*#__PURE__*/React.createElement("button", {
    onClick: () => {
      const plotEl = document.getElementById("main-plot");
      if (plotEl) downloadPlotSvg(plotEl);
    },
    className: "text-xs px-3 py-1 bg-emerald-600 hover:bg-emerald-700 text-white rounded-lg transition",
    title: "Vector SVG - 210 x 99 mm (A4 page width x one-third page height)"
  }, "Export SVG"), /*#__PURE__*/React.createElement("button", {
    onClick: () => {
      Plotly.downloadImage("main-plot", {
        ...FIGURE.plotPng,
        filename: "quadrant_plot"
      });
    },
    className: "text-xs px-3 py-1 bg-amber-500 hover:bg-amber-600 text-white rounded-lg transition",
    title: "High-resolution raster preview"
  }, "Export PNG"))), /*#__PURE__*/React.createElement("div", {
    id: "main-plot",
    className: "h-[480px] min-h-[480px] shrink-0 screen-only"
  }), plotPrintImg && /*#__PURE__*/React.createElement("div", {
    className: "print-only"
  }, /*#__PURE__*/React.createElement("img", {
    src: plotPrintImg,
    alt: "Step 1 Signal Discovery Plot",
    className: "print-figure print-figure--vector print-figure--a4-portrait-third"
  }))), /*#__PURE__*/React.createElement("div", {
    className: "col-span-12 lg:col-span-4 bg-white rounded-xl p-4 border border-gray-200 shadow-sm flex flex-col"
  }, /*#__PURE__*/React.createElement("h3", {
    className: "text-lg mb-2 text-gray-900"
  }, "Candidates (", candidates.length, " total)"), /*#__PURE__*/React.createElement("div", {
    className: "space-y-2 max-h-[480px] overflow-y-auto pr-1 flex-1"
  }, (() => {
    let rank = 0;
    return candidateGroups.map(group => /*#__PURE__*/React.createElement("div", {
      key: group.quadrant,
      className: "space-y-2"
    }, /*#__PURE__*/React.createElement("div", {
      className: "sticky top-0 z-10 bg-white/95 backdrop-blur px-1 py-1 text-xs font-semibold text-gray-600 border-b border-gray-100"
    }, QUADRANT_HEADERS[group.quadrant] || `Quadrant ${group.quadrant}`, /*#__PURE__*/React.createElement("span", {
      className: "ml-1 font-normal text-gray-400"
    }, "(", group.items.length, ")")), group.items.map(c => {
      rank += 1;
      const isSelected = selectedPair?.bacteria === c.bacteria && selectedPair?.metabolite === c.metabolite;
      return /*#__PURE__*/React.createElement("div", {
        key: `${c.bacteria}-${c.metabolite}`,
        onClick: () => setSelectedPair(c),
        className: `p-2 rounded-lg border cursor-pointer transition ${isSelected ? "border-indigo-400 bg-indigo-50" : "border-gray-200 hover:border-gray-400 hover:bg-gray-50"}`
      }, /*#__PURE__*/React.createElement("div", {
        className: "flex items-center justify-between"
      }, /*#__PURE__*/React.createElement("div", {
        className: "text-sm font-medium text-gray-800 truncate"
      }, c.metabolite), /*#__PURE__*/React.createElement("span", {
        className: "text-xs text-gray-400"
      }, "#", rank)), /*#__PURE__*/React.createElement("div", {
        className: "text-xs text-gray-500 flex gap-3 mt-1 flex-wrap"
      }, /*#__PURE__*/React.createElement("span", null, "MMSage: ", c.mmsage_norm.toFixed(3)), /*#__PURE__*/React.createElement("span", null, "bm: ", c.pair_bm_exp ?? 0), /*#__PURE__*/React.createElement("span", null, "md: ", c.pair_md_exp ?? 0), /*#__PURE__*/React.createElement("span", {
        className: `font-medium ${(c.novelty_score ?? 1) >= 1.0 ? "text-amber-500" : "text-emerald-600"}`
      }, (c.novelty_score ?? 1) >= 1.0 ? "★ novel" : "known")));
    })));
  })()), /*#__PURE__*/React.createElement("button", {
    className: "mt-4 w-full bg-indigo-600 hover:bg-indigo-700 disabled:bg-gray-400 text-white py-2 rounded-lg font-medium shadow-sm transition mt-auto",
    onClick: handleAnalyze,
    disabled: !selectedPair || loading
  }, loading ? "Generating..." : "Generate Step 3 Protocol"))), /*#__PURE__*/React.createElement("div", {
    id: "uncharted-dark-matter",
    className: "col-span-12 bg-gray-50/80 rounded-lg px-3 py-2.5 border border-gray-100"
  }, /*#__PURE__*/React.createElement("div", {
    className: "flex flex-wrap items-center gap-x-3 gap-y-1"
  }, /*#__PURE__*/React.createElement("h3", {
    className: "text-xs font-medium text-gray-500"
  }, "Uncharted Dark Matter", /*#__PURE__*/React.createElement("span", {
    className: "ml-1 text-gray-400"
  }, "(", unchartedDarkMatter.length, ")")), /*#__PURE__*/React.createElement("span", {
    className: "text-[11px] text-gray-400"
  }, "bm = 0, md = 0 · excluded from Step 1 plot"), /*#__PURE__*/React.createElement("details", {
    className: "text-[11px] text-gray-400 ml-auto"
  }, /*#__PURE__*/React.createElement("summary", {
    className: "cursor-pointer hover:text-gray-600 select-none list-none [&::-webkit-details-marker]:hidden"
  }, "Why separate?"), /*#__PURE__*/React.createElement("p", {
    className: "mt-1.5 text-gray-500 max-w-2xl leading-relaxed"
  }, "No PubMed evidence for microbe–metabolite (bm) or metabolite–disease (md) links; zero-literature pairs are omitted from the scatter plot but remain high-priority MMSage signals."))), unchartedDarkMatter.length === 0 ? /*#__PURE__*/React.createElement("p", {
    className: "text-xs text-gray-400 mt-2"
  }, "No uncharted candidates in this run.") : /*#__PURE__*/React.createElement(React.Fragment, null, /*#__PURE__*/React.createElement("div", {
    className: "uncharted-grid grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-1.5 mt-2"
  }, (unchartedExpanded ? unchartedDarkMatter : unchartedDarkMatter.slice(0, UNCHARTED_PREVIEW)).map((c, i) => {
    const isSelected = selectedPair?.bacteria === c.bacteria && selectedPair?.metabolite === c.metabolite;
    return /*#__PURE__*/React.createElement("div", {
      key: `${c.bacteria}-${c.metabolite}`,
      onClick: () => setSelectedPair(c),
      className: `px-2 py-1.5 rounded border cursor-pointer transition text-xs ${isSelected ? "border-gray-400 bg-white" : "border-gray-200 bg-white/60 hover:border-gray-300 hover:bg-white"}`
    }, /*#__PURE__*/React.createElement("div", {
      className: "flex items-center justify-between gap-2"
    }, /*#__PURE__*/React.createElement("div", {
      className: "font-medium text-gray-700 truncate"
    }, c.metabolite), /*#__PURE__*/React.createElement("span", {
      className: "text-gray-400 shrink-0"
    }, "#", i + 1)), /*#__PURE__*/React.createElement("div", {
      className: "text-gray-400 flex gap-2 mt-0.5"
    }, /*#__PURE__*/React.createElement("span", null, "MMSage ", c.mmsage_norm.toFixed(3)), /*#__PURE__*/React.createElement("span", null, "uncharted")));
  })), unchartedDarkMatter.length > UNCHARTED_PREVIEW && /*#__PURE__*/React.createElement("button", {
    type: "button",
    onClick: () => setUnchartedExpanded(v => !v),
    className: "mt-2 text-[11px] text-gray-400 hover:text-gray-600 transition"
  }, unchartedExpanded ? "Show less" : `Show all ${unchartedDarkMatter.length} uncharted pairs`))), /*#__PURE__*/React.createElement("div", {
    className: "grid grid-cols-12 gap-6"
  }, /*#__PURE__*/React.createElement("div", {
    id: "step2-card",
    className: "col-span-12 lg:col-span-8 bg-white rounded-xl p-4 border border-gray-200 shadow-sm flex flex-col"
  }, /*#__PURE__*/React.createElement("div", {
    className: "flex items-center justify-between mb-2"
  }, /*#__PURE__*/React.createElement("h2", {
    className: "text-2xl text-gray-900 font-semibold"
  }, "Step 2: Knowledge Graph"), /*#__PURE__*/React.createElement("div", {
    className: "flex gap-2"
  }, /*#__PURE__*/React.createElement("button", {
    onClick: () => downloadKgSvg(window._kgCy),
    className: "text-xs px-3 py-1 bg-emerald-600 hover:bg-emerald-700 text-white rounded-lg transition",
    title: "Vector SVG for journal figures"
  }, "Export SVG"), /*#__PURE__*/React.createElement("button", {
    onClick: () => {
      if (!window._kgCy) return;
      const png = window._kgCy.png({
        output: "blob",
        bg: "#ffffff",
        scale: FIGURE.kgPngScale,
        full: true
      });
      const a = document.createElement("a");
      a.href = URL.createObjectURL(png);
      a.download = "kg_chain_view.png";
      a.click();
      URL.revokeObjectURL(a.href);
    },
    className: "text-xs px-3 py-1 bg-amber-500 hover:bg-amber-600 text-white rounded-lg transition",
    title: "Raster fallback if SVG preview differs"
  }, "Export PNG"))), /*#__PURE__*/React.createElement("div", {
    id: "kg",
    className: "flex-1 min-h-[480px] bg-gray-50 rounded-lg border border-gray-100 screen-only"
  }), /*#__PURE__*/React.createElement("div", {
    className: "mt-3 flex flex-wrap gap-x-4 gap-y-2 text-xs text-gray-600 screen-only"
  }, /*#__PURE__*/React.createElement("div", {
    className: "flex items-center gap-2"
  }, /*#__PURE__*/React.createElement("span", {
    className: "inline-block w-3.5 h-3.5 rounded-full bg-blue-500 border border-white shadow-sm"
  }), /*#__PURE__*/React.createElement("span", null, "Microbe")), /*#__PURE__*/React.createElement("div", {
    className: "flex items-center gap-2"
  }, /*#__PURE__*/React.createElement("span", {
    className: "inline-block w-4 h-3 rounded-sm bg-amber-500 border border-white shadow-sm"
  }), /*#__PURE__*/React.createElement("span", null, "Metabolite")), /*#__PURE__*/React.createElement("div", {
    className: "flex items-center gap-2"
  }, /*#__PURE__*/React.createElement("span", {
    className: "inline-block w-3.5 h-3.5 bg-emerald-500 border border-white shadow-sm rotate-45"
  }), /*#__PURE__*/React.createElement("span", null, "Pathway")), /*#__PURE__*/React.createElement("div", {
    className: "flex items-center gap-2"
  }, /*#__PURE__*/React.createElement("span", {
    className: "inline-flex items-center justify-center w-4 h-4 bg-red-500 text-white text-[9px] leading-none shadow-sm",
    style: {
      clipPath: "polygon(25% 0%, 75% 0%, 100% 50%, 75% 100%, 25% 100%, 0% 50%)"
    }
  }), /*#__PURE__*/React.createElement("span", null, "Disease")), /*#__PURE__*/React.createElement("div", {
    className: "flex items-center gap-2"
  }, /*#__PURE__*/React.createElement("span", {
    className: "inline-block w-0 h-0 border-l-[8px] border-r-[8px] border-b-[14px] border-l-transparent border-r-transparent border-b-orange-500 drop-shadow-sm"
  }), /*#__PURE__*/React.createElement("span", null, "Enzyme")), /*#__PURE__*/React.createElement("div", {
    className: "flex items-center gap-2"
  }, /*#__PURE__*/React.createElement("span", {
    className: "inline-block w-0 h-0 border-l-[8px] border-r-[8px] border-t-[14px] border-l-transparent border-r-transparent border-t-cyan-500 drop-shadow-sm"
  }), /*#__PURE__*/React.createElement("span", null, "Receptor")), /*#__PURE__*/React.createElement("div", {
    className: "flex items-center gap-2"
  }, /*#__PURE__*/React.createElement("span", {
    className: "inline-block w-3.5 h-3.5 rounded-full bg-gray-400 border border-white shadow-sm"
  }), /*#__PURE__*/React.createElement("span", null, "Unknown"))), /*#__PURE__*/React.createElement("p", {
    className: "text-xs text-gray-400 mt-2 screen-only"
  }, "Click a node to view literature evidence below."), /*#__PURE__*/React.createElement("div", {
    id: "kg-evidence-panel",
    className: "hidden mt-4 border-t border-gray-100 pt-4 screen-only"
  }, /*#__PURE__*/React.createElement("h3", {
    className: "text-sm font-semibold text-gray-800 mb-2"
  }, "Literature evidence"), /*#__PURE__*/React.createElement("div", {
    id: "kg-evidence-content",
    className: "text-sm max-h-[360px] overflow-y-auto pr-1"
  })), kgPrintImg && /*#__PURE__*/React.createElement("div", {
    className: "print-only"
  }, /*#__PURE__*/React.createElement("img", {
    src: kgPrintImg,
    alt: "Step 2 Knowledge Graph",
    className: "print-figure print-figure--vector print-figure--a4-landscape"
  }))), /*#__PURE__*/React.createElement("div", {
    className: "col-span-12 lg:col-span-4 bg-white rounded-xl p-4 border border-gray-200 shadow-sm flex flex-col"
  }, /*#__PURE__*/React.createElement("h2", {
    className: "text-xl mb-2 text-gray-900"
  }, "Analysis"), selectedPair ? /*#__PURE__*/React.createElement("div", {
    className: "flex-1 flex flex-col space-y-4"
  }, /*#__PURE__*/React.createElement("div", {
    className: "text-sm"
  }, /*#__PURE__*/React.createElement("div", {
    className: "flex justify-between mb-2"
  }, /*#__PURE__*/React.createElement("span", {
    className: "text-gray-500"
  }, "MMSage Signal:"), /*#__PURE__*/React.createElement("span", {
    className: "font-mono text-gray-800"
  }, selectedPair.mmsage_norm.toFixed(3))), /*#__PURE__*/React.createElement("div", {
    className: "flex justify-between mb-2"
  }, /*#__PURE__*/React.createElement("span", {
    className: "text-gray-500"
  }, "Microbe-Metabolite (bm):"), /*#__PURE__*/React.createElement("span", {
    className: "font-mono text-gray-800"
  }, selectedPair.pair_bm_exp ?? 0, " papers")), /*#__PURE__*/React.createElement("div", {
    className: "flex justify-between mb-2"
  }, /*#__PURE__*/React.createElement("span", {
    className: "text-gray-500"
  }, "Metabolite-Disease (md):"), /*#__PURE__*/React.createElement("span", {
    className: "font-mono text-gray-800"
  }, selectedPair.pair_md_exp ?? 0, " papers")), /*#__PURE__*/React.createElement("div", {
    className: "flex justify-between font-semibold border-t border-gray-200 pt-2 mt-2"
  }, /*#__PURE__*/React.createElement("span", {
    className: "text-gray-700"
  }, "Composite Score:"), /*#__PURE__*/React.createElement("span", {
    className: "font-mono text-indigo-600"
  }, (selectedPair.composite_score ?? 0).toFixed(2))), /*#__PURE__*/React.createElement("div", {
    className: "flex justify-between mt-2"
  }, /*#__PURE__*/React.createElement("span", {
    className: "text-gray-500"
  }, "Novelty:"), /*#__PURE__*/React.createElement("span", {
    className: "font-mono",
    style: {
      color: (selectedPair.novelty_score ?? 1) >= 1.0 ? "#d97706" : "#059669"
    }
  }, (selectedPair.novelty_score ?? 0).toFixed(3), " ", (selectedPair.novelty_score ?? 1) >= 1.0 ? "★ novel" : "known"))), analysisResult && /*#__PURE__*/React.createElement("div", {
    className: "bg-indigo-50 rounded-lg p-3 border-l-4 border-indigo-500 mt-auto"
  }, /*#__PURE__*/React.createElement("div", {
    className: "text-xs font-semibold uppercase text-indigo-600 mb-1"
  }, "Prioritization"), /*#__PURE__*/React.createElement("div", {
    className: "text-sm text-indigo-900 font-medium"
  }, analysisResult.prioritization))) : /*#__PURE__*/React.createElement("div", {
    className: "text-gray-400 text-center flex-1 flex items-center justify-center"
  }, "Select a candidate to view analysis")))), /*#__PURE__*/React.createElement("div", {
    id: "step3-card",
    className: "mt-6 bg-white rounded-xl p-5 border border-gray-200 shadow-sm"
  }, /*#__PURE__*/React.createElement("div", {
    className: "border-b border-gray-100 pb-2 mb-3"
  }, /*#__PURE__*/React.createElement("h2", {
    className: "text-xl text-gray-900 font-semibold"
  }, "Step 3: Evidence Gap and Follow-up Plan"), /*#__PURE__*/React.createElement("p", {
    className: "mt-1 text-sm text-gray-500"
  }, "First identify what is already known and what remains unknown, then propose the most useful follow-up experiment(s).")), loading && /*#__PURE__*/React.createElement("p", {
    className: "text-sm text-indigo-600 mb-2"
  }, "Generating validation protocol..."), protocolError && /*#__PURE__*/React.createElement("div", {
    className: "mb-3 rounded-lg border border-red-200 bg-red-50 px-3 py-2 text-sm text-red-700"
  }, protocolError), (loading || protocol) && /*#__PURE__*/React.createElement("pre", {
    className: "whitespace-pre-wrap text-base leading-relaxed text-gray-800 max-h-[600px] overflow-y-auto bg-gray-50/50 p-5 rounded-lg font-serif"
  }, loading ? "Please wait while the validation protocol is being generated..." : protocol), validationPlanStatus && (validationPlanLoading || validationPlanStatus.status === "error") && /*#__PURE__*/React.createElement("div", {
    className: "mt-4 rounded-xl border border-indigo-200 bg-indigo-50/70 p-4"
  }, /*#__PURE__*/React.createElement("div", {
    className: "flex items-start justify-between gap-3 flex-wrap"
  }, /*#__PURE__*/React.createElement("div", null, /*#__PURE__*/React.createElement("h3", {
    className: "text-sm font-semibold text-indigo-900"
  }, validationPlanStatus?.status === "error" ? "Step 3 Generation Failed" : "Step 3 Generation Running"), /*#__PURE__*/React.createElement("p", {
    className: "mt-1 text-sm text-indigo-800"
  }, validationPlanStatus?.current_message || "Generating Step 3 evidence-gap summary..."), validationPlanJobId && /*#__PURE__*/React.createElement("p", {
    className: "mt-1 text-xs text-indigo-600 font-mono"
  }, "Job: ", validationPlanJobId)), /*#__PURE__*/React.createElement("div", {
    className: "min-w-[120px] text-right"
  }, /*#__PURE__*/React.createElement("div", {
    className: "text-lg font-semibold text-indigo-900"
  }, validationPlanStatus?.progress_percent ?? 0, "%"), /*#__PURE__*/React.createElement("div", {
    className: "text-xs text-indigo-700"
  }, formatValidationStage(validationPlanStatus?.current_stage)))), /*#__PURE__*/React.createElement("div", {
    className: "mt-3 h-2 w-full overflow-hidden rounded-full bg-indigo-100"
  }, /*#__PURE__*/React.createElement("div", {
    className: "h-full rounded-full bg-indigo-600 transition-all",
    style: {
      width: `${Math.max(6, validationPlanStatus?.progress_percent ?? 6)}%`
    }
  })), /*#__PURE__*/React.createElement("div", {
    className: "mt-4 grid gap-4 lg:grid-cols-[220px_minmax(0,1fr)]"
  }, /*#__PURE__*/React.createElement("div", {
    className: "rounded-lg border border-indigo-200 bg-white/80 p-3"
  }, /*#__PURE__*/React.createElement("div", {
    className: "text-xs font-semibold uppercase tracking-wide text-indigo-800"
  }, "Current Stage"), /*#__PURE__*/React.createElement("div", {
    className: "mt-2 text-sm font-medium text-gray-900"
  }, formatValidationStage(validationPlanStatus?.current_stage)), validationPlanStatus?.candidate?.bacteria && /*#__PURE__*/React.createElement("div", {
    className: "mt-3 text-xs text-gray-600 space-y-1"
  }, /*#__PURE__*/React.createElement("div", null, "Microbe: ", validationPlanStatus.candidate.bacteria), /*#__PURE__*/React.createElement("div", null, "Metabolite: ", validationPlanStatus.candidate.metabolite), /*#__PURE__*/React.createElement("div", null, "Disease: ", validationPlanStatus.candidate.disease))), /*#__PURE__*/React.createElement("div", {
    className: "rounded-lg border border-indigo-200 bg-white/80 p-3"
  }, /*#__PURE__*/React.createElement("div", {
    className: "text-xs font-semibold uppercase tracking-wide text-indigo-800"
  }, "Live Logs"), Array.isArray(validationPlanStatus?.logs) && validationPlanStatus.logs.length > 0 ? /*#__PURE__*/React.createElement("div", {
    className: "mt-2 space-y-2 max-h-56 overflow-y-auto pr-1"
  }, validationPlanStatus.logs.slice(-8).reverse().map((item, idx) => /*#__PURE__*/React.createElement("div", {
    key: `validation-log-${idx}`,
    className: "rounded-lg border border-gray-200 bg-gray-50/80 p-2.5"
  }, /*#__PURE__*/React.createElement("div", {
    className: "flex items-center gap-2 flex-wrap"
  }, /*#__PURE__*/React.createElement("span", {
    className: "text-xs font-semibold text-gray-900"
  }, formatValidationStage(item.stage)), /*#__PURE__*/React.createElement("span", {
    className: "text-xs text-gray-500"
  }, item.progress_percent ?? 0, "%"), item.timestamp && /*#__PURE__*/React.createElement("span", {
    className: "text-xs text-gray-400"
  }, new Date(item.timestamp).toLocaleTimeString())), /*#__PURE__*/React.createElement("p", {
    className: "mt-1 text-sm text-gray-700"
  }, item.message), item.extra && Object.keys(item.extra).length > 0 && /*#__PURE__*/React.createElement("p", {
    className: "mt-1 text-xs text-gray-500"
  }, Object.entries(item.extra).map(([k, v]) => `${k}: ${v}`).join(" | "))))) : /*#__PURE__*/React.createElement("p", {
    className: "mt-2 text-sm text-gray-500"
  }, "Waiting for backend logs...")))), validationPlanError && /*#__PURE__*/React.createElement("div", {
    className: "mt-4 rounded-lg border border-red-200 bg-red-50 px-3 py-2 text-sm text-red-700"
  }, validationPlanError), /*#__PURE__*/React.createElement("div", {
    className: "mt-6 bg-slate-50/80 rounded-xl p-5 border border-slate-200 shadow-sm"
  }, /*#__PURE__*/React.createElement("div", {
    className: "border-b border-slate-200 pb-2 mb-3"
  }, /*#__PURE__*/React.createElement("h2", {
    className: "text-xl text-slate-900 font-semibold"
  }, "Specific Experiment Question"), /*#__PURE__*/React.createElement("p", {
    className: "mt-1 text-sm text-slate-600"
  }, "Use this independent module to ask about any concrete experiment, even if it is not related to the selected candidate above. Provide the claim, model system, or constraints if available.")), /*#__PURE__*/React.createElement("div", {
    className: "grid gap-3"
  }, /*#__PURE__*/React.createElement("div", null, /*#__PURE__*/React.createElement("label", {
    htmlFor: "standalone-validation-question",
    className: "block text-xs font-semibold uppercase tracking-wide text-slate-600"
  }, "EXPERIMENT QUESTION"), /*#__PURE__*/React.createElement("textarea", {
    id: "standalone-validation-question",
    value: questionPlanResearchQuestion,
    onChange: e => setQuestionPlanResearchQuestion(e.target.value),
    placeholder: "Example: How should we validate whether metabolite X suppresses epithelial inflammation in Caco-2 cells or intestinal organoids? Include rationale, groups, controls, readouts, expected outcomes, and decision rules.",
    className: "mt-1 min-h-[96px] w-full rounded-lg border border-slate-300 bg-white px-3 py-2 text-sm text-slate-800 shadow-sm outline-none transition focus:border-slate-500 focus:ring-2 focus:ring-slate-100"
  })), /*#__PURE__*/React.createElement("div", null, /*#__PURE__*/React.createElement("label", {
    htmlFor: "standalone-validation-constraints",
    className: "block text-xs font-semibold uppercase tracking-wide text-slate-600"
  }, "OPTIONAL CONTEXT / CONSTRAINTS"), /*#__PURE__*/React.createElement("textarea", {
    id: "standalone-validation-constraints",
    value: questionPlanPromptConstraints,
    onChange: e => setQuestionPlanPromptConstraints(e.target.value),
    placeholder: "Example: Keep only the most basic in vitro and in vivo experiments. Avoid mechanism-heavy studies and use literature-supported conditions when available.",
    className: "mt-1 min-h-[84px] w-full rounded-lg border border-slate-300 bg-white px-3 py-2 text-sm text-slate-800 shadow-sm outline-none transition focus:border-slate-500 focus:ring-2 focus:ring-slate-100"
  })), /*#__PURE__*/React.createElement("div", {
    className: "flex items-center justify-end"
  }, /*#__PURE__*/React.createElement("button", {
    onClick: handleGenerateValidationPlan,
    disabled: questionPlanLoading || !questionPlanResearchQuestion.trim(),
    className: `px-4 py-2 rounded-lg text-sm font-medium transition ${!questionPlanLoading && questionPlanResearchQuestion.trim() ? "bg-slate-800 hover:bg-slate-900 text-white" : "bg-slate-200 text-slate-400 cursor-not-allowed"}`
  }, questionPlanLoading ? "Generating..." : "Generate Standalone Protocol"))), questionPlanStatus && (questionPlanLoading || questionPlanStatus.status === "error") && /*#__PURE__*/React.createElement("div", {
    className: "mt-4 rounded-xl border border-slate-200 bg-white/80 p-4"
  }, /*#__PURE__*/React.createElement("div", {
    className: "flex items-start justify-between gap-3 flex-wrap"
  }, /*#__PURE__*/React.createElement("div", null, /*#__PURE__*/React.createElement("h3", {
    className: "text-sm font-semibold text-slate-900"
  }, questionPlanStatus?.status === "error" ? "Standalone Generation Failed" : "Standalone Generation Running"), /*#__PURE__*/React.createElement("p", {
    className: "mt-1 text-sm text-slate-700"
  }, questionPlanStatus?.current_message || "Generating a standalone evidence-driven protocol..."), questionPlanJobId && /*#__PURE__*/React.createElement("p", {
    className: "mt-1 text-xs text-slate-500 font-mono"
  }, "Job: ", questionPlanJobId)), /*#__PURE__*/React.createElement("div", {
    className: "min-w-[120px] text-right"
  }, /*#__PURE__*/React.createElement("div", {
    className: "text-lg font-semibold text-slate-900"
  }, questionPlanStatus?.progress_percent ?? 0, "%"), /*#__PURE__*/React.createElement("div", {
    className: "text-xs text-slate-600"
  }, formatValidationStage(questionPlanStatus?.current_stage)))), /*#__PURE__*/React.createElement("div", {
    className: "mt-3 h-2 w-full overflow-hidden rounded-full bg-slate-100"
  }, /*#__PURE__*/React.createElement("div", {
    className: "h-full rounded-full bg-slate-700 transition-all",
    style: {
      width: `${Math.max(6, questionPlanStatus?.progress_percent ?? 6)}%`
    }
  })), Array.isArray(questionPlanStatus?.logs) && questionPlanStatus.logs.length > 0 ? /*#__PURE__*/React.createElement("div", {
    className: "mt-3 space-y-2 max-h-56 overflow-y-auto pr-1"
  }, questionPlanStatus.logs.slice(-8).reverse().map((item, idx) => /*#__PURE__*/React.createElement("div", {
    key: `question-log-${idx}`,
    className: "rounded-lg border border-slate-200 bg-slate-50 p-2.5"
  }, /*#__PURE__*/React.createElement("div", {
    className: "flex items-center gap-2 flex-wrap"
  }, /*#__PURE__*/React.createElement("span", {
    className: "text-xs font-semibold text-slate-900"
  }, formatValidationStage(item.stage)), /*#__PURE__*/React.createElement("span", {
    className: "text-xs text-slate-500"
  }, item.progress_percent ?? 0, "%"), item.timestamp && /*#__PURE__*/React.createElement("span", {
    className: "text-xs text-slate-400"
  }, new Date(item.timestamp).toLocaleTimeString())), /*#__PURE__*/React.createElement("p", {
    className: "mt-1 text-sm text-slate-700"
  }, item.message), item.extra && Object.keys(item.extra).length > 0 && /*#__PURE__*/React.createElement("p", {
    className: "mt-1 text-xs text-slate-500"
  }, Object.entries(item.extra).map(([k, v]) => `${k}: ${v}`).join(" | "))))) : /*#__PURE__*/React.createElement("p", {
    className: "mt-3 text-sm text-slate-500"
  }, "Waiting for backend logs...")), questionPlanError && /*#__PURE__*/React.createElement("div", {
    className: "mt-4 rounded-lg border border-red-200 bg-red-50 px-3 py-2 text-sm text-red-700"
  }, questionPlanError), questionPlan && /*#__PURE__*/React.createElement("pre", {
    className: "mt-4 whitespace-pre-wrap text-base leading-relaxed text-slate-800 max-h-[720px] overflow-y-auto rounded-lg bg-white p-5 font-serif"
  }, questionPlan.validation_protocol_text || "Standalone validation protocol text was not returned."))), /*#__PURE__*/React.createElement("footer", {
    className: "mt-8 text-center text-gray-400 text-sm"
  }, "MMSage-Insight | X: MMSage Signal | Y: Evidence Foundation (bm_exp + md_exp) | ★ = Novel (chain_novelty=1, no triple experimental co-occurrence) | Composite = MMSage × (1+log₂(1+bm)) × (1+log₂(1+md))"));
}
ReactDOM.createRoot(document.getElementById("root")).render(/*#__PURE__*/React.createElement(App, null));
