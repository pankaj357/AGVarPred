
import { useEffect, useState } from "react";
import {
  ArrowDownRight,
  ArrowRight,
  CheckCircle2,
  CircleHelp,
  Dna,
  FileUp,
  FlaskConical,
  Github,
  LoaderCircle,
  Plus,
  ShieldCheck,
  Sparkles,
  Trash2,
  TriangleAlert,
  X,
} from "lucide-react";

type Variant = { chrom: string; pos: string; ref: string; alt: string; gene: string };
type Health = {
  status: string;
  model_version: string;
  model_type: string;
  feature_count: number;
  threshold: number;
  live_prediction_enabled: boolean;
  live_provider_ready: boolean;
};
type ApiError = { detail?: string };
type PredictionResult = {
  variant_id: string;
  chrom: string;
  pos: number;
  ref: string;
  alt: string;
  gene?: string | null;
  status: "ok" | "error";
  score?: number;
  predicted_class?: string;
  threshold?: number;
  gnomAD_AF?: number | null;
  consequence?: string | null;
  error?: string;
};
type PredictResponse = ApiError & {
  results?: PredictionResult[];
  requested_variants?: number;
  scored_variants?: number;
  model_version?: string;
  notice?: string;
};

const emptyVariant = (): Variant => ({ chrom: "chr", pos: "", ref: "", alt: "", gene: "" });
const STATIC_PREVIEW = import.meta.env.VITE_STATIC_PREVIEW === "true";
const API_BASE = STATIC_PREVIEW
  ? ""
  : (import.meta.env.VITE_API_BASE_URL ?? "http://127.0.0.1:8000").replace(/\/$/, "");

function App() {
  const [variants, setVariants] = useState<Variant[]>([{ ...emptyVariant(), chrom: "chr1" }]);
  const [busy, setBusy] = useState(false);
  const [notice, setNotice] = useState<{ type: "error" | "info"; message: string } | null>(null);
  const [health, setHealth] = useState<Health | null>(null);
  const [results, setResults] = useState<PredictionResult[]>([]);
  const [resultMeta, setResultMeta] = useState<{ model_version?: string; notice?: string } | null>(null);
  const [apiUnavailable, setApiUnavailable] = useState(false);

  useEffect(() => {
    let active = true;
    if (STATIC_PREVIEW) {
      if (active) {
        setHealth({
          status: "ok",
          model_version: "model_full",
          model_type: "full",
          feature_count: 120,
          threshold: 0.42,
          live_prediction_enabled: false,
          live_provider_ready: false,
        });
        setApiUnavailable(true);
      }
      return;
    }
    fetch(`${API_BASE}/api/health`)
      .then(async (response) => {
        if (!response.ok) throw new Error("API health endpoint returned an error.");
        return response.json() as Promise<Health>;
      })
      .then((data) => {
        if (active) { setHealth(data); setApiUnavailable(false); }
      })
      .catch(() => {
        if (active) { setApiUnavailable(true); setHealth(null); }
      });
    return () => { active = false; };
  }, []);

  const canSubmit = variants.length > 0 && variants.every(
    (v) => v.chrom.trim() && /^\d+$/.test(v.pos) && Number(v.pos) > 0 &&
      /^[ACGTN]+$/i.test(v.ref) && /^[ACGTN]+$/i.test(v.alt) && v.ref.toUpperCase() !== v.alt.toUpperCase(),
  );

  const updateVariant = (index: number, key: keyof Variant, value: string) => {
    setVariants((current) => current.map((variant, i) => i === index ? { ...variant, [key]: value } : variant));
  };

  const addVariant = () => {
    if (variants.length < 25) setVariants((current) => [...current, emptyVariant()]);
  };

  const removeVariant = (index: number) => {
    setVariants((current) => current.length === 1 ? current : current.filter((_, i) => i !== index));
  };

  const submit = async (event: React.FormEvent) => {
    event.preventDefault();
    if (!canSubmit) {
      setNotice({ type: "error", message: "Check every row: use a GRCh38 chromosome, a positive 1-based position, and different A/C/G/T/N alleles." });
      return;
    }
    if (STATIC_PREVIEW) {
      setNotice({ type: "info", message: "Online predictions are not yet available in this static preview." });
      return;
    }
    setBusy(true);
    setNotice(null);
    setResults([]);
    setResultMeta(null);
    try {
      const response = await fetch(`${API_BASE}/api/predict`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          genome_build: "GRCh38",
          variants: variants.map((v) => ({
            chrom: v.chrom,
            pos: Number(v.pos),
            ref: v.ref.toUpperCase(),
            alt: v.alt.toUpperCase(),
            ...(v.gene.trim() ? { gene: v.gene.trim() } : {}),
          })),
        }),
      });
      const payload = (await response.json()) as PredictResponse;
      if (!response.ok) throw new Error(payload.detail ?? `Request failed (${response.status}).`);
      setResults(payload.results ?? []);
      setResultMeta({ model_version: payload.model_version, notice: payload.notice });
      setNotice({
        type: "info",
        message: `${payload.scored_variants ?? 0} of ${payload.requested_variants ?? variants.length} variants scored.`,
      });
    } catch (error) {
      setNotice({
        type: "error",
        message: error instanceof Error ? error.message : "The request could not be completed.",
      });
    } finally {
      setBusy(false);
    }
  };

  return (
    <div className="min-h-screen bg-ink text-emerald-50">
      <div className="pointer-events-none fixed inset-0 grid-fade opacity-70" />
      <header className="relative z-10 border-b border-white/10">
        <div className="mx-auto flex max-w-7xl items-center justify-between px-5 py-5 sm:px-8">
          <a href="#" className="flex items-center gap-3 no-underline">
            <span className="flex h-10 w-10 items-center justify-center rounded-xl border border-mint/30 bg-mint/10 text-mint">
              <Dna size={22} />
            </span>
            <span>
              <span className="block text-sm font-semibold tracking-wide">AGVarPred</span>
              <span className="block text-xs text-emerald-100/45">Variant interpretation research</span>
            </span>
          </a>
          <a
            className="hidden items-center gap-2 rounded-full border border-white/15 px-4 py-2 text-sm text-emerald-50/75 transition hover:border-mint/40 hover:text-mint sm:inline-flex"
            href="https://github.com/pankaj357/AGVarPred"
            target="_blank"
            rel="noreferrer"
          >
            <Github size={15} /> Source code <ArrowRight size={14} />
          </a>
        </div>
      </header>

      <main className="relative z-10 mx-auto max-w-7xl px-5 pb-16 pt-12 sm:px-8 sm:pt-20">
        <section className="grid items-start gap-12 lg:grid-cols-[1.05fr_.95fr] lg:gap-16">
          <div className="pt-2">
            <div className="mb-6 inline-flex items-center gap-2 rounded-full border border-mint/20 bg-mint/[0.07] px-3 py-1.5 text-xs font-medium tracking-wide text-mint">
              <Sparkles size={13} /> ALPHAGENOME-INFORMED PREDICTION
            </div>
            <h1 className="max-w-2xl text-4xl font-semibold leading-[1.08] tracking-[-0.04em] text-white sm:text-5xl lg:text-6xl">
              Explore the functional context of a <span className="text-mint">variant.</span>
            </h1>
            <p className="mt-6 max-w-xl text-base leading-7 text-emerald-50/65 sm:text-lg">
              AGVarPred integrates functional-genomics scores, population allele frequency and variant annotations to support research on germline variant pathogenicity.
            </p>
            <div className="mt-8 flex flex-wrap items-center gap-3 text-sm text-emerald-50/65">
              <span className="inline-flex items-center gap-2 rounded-full border border-white/10 bg-white/[0.03] px-3 py-2"><ShieldCheck size={15} className="text-aqua" /> GRCh38 input</span>
              <span className="inline-flex items-center gap-2 rounded-full border border-white/10 bg-white/[0.03] px-3 py-2"><FlaskConical size={15} className="text-aqua" /> Research use</span>
              <span className="inline-flex items-center gap-2 rounded-full border border-white/10 bg-white/[0.03] px-3 py-2"><CheckCircle2 size={15} className="text-aqua" /> Frozen model</span>
            </div>

            <div className="mt-12 grid max-w-xl grid-cols-3 gap-3">
              <Metric label="Selected features" value={health ? String(health.feature_count) : "120"} />
              <Metric label="Model threshold" value={health ? health.threshold.toFixed(2) : "0.42"} />
              <Metric label="Genome build" value="GRCh38" />
            </div>
          </div>

          <section className="rounded-3xl border border-white/10 bg-panel/90 p-5 shadow-glow backdrop-blur sm:p-7">
            <div className="flex items-start justify-between gap-4">
              <div>
                <div className="text-xs font-semibold uppercase tracking-[0.16em] text-mint/80">Start an analysis</div>
                <h2 className="mt-2 text-2xl font-semibold tracking-tight text-white">Enter variants</h2>
                <p className="mt-2 text-sm leading-6 text-emerald-50/50">One-based reference coordinates. Small variants only; GRCh38.</p>
              </div>
              <div className="rounded-xl border border-white/10 bg-white/[0.04] p-3 text-mint"><FileUp size={21} /></div>
            </div>

            {STATIC_PREVIEW && (
              <div className="mt-5 flex items-start gap-3 rounded-xl border border-amber-200/20 bg-amber-100/[0.06] p-4 text-sm leading-6 text-amber-100/80">
                <TriangleAlert size={18} className="mt-0.5 shrink-0" />
                <span>
                  <strong>Static preview.</strong> Online predictions are not yet available from this GitHub Pages site.
                  The form below demonstrates the interface; submissions are disabled until a backend is deployed.
                </span>
              </div>
            )}

            <form className="mt-7 space-y-4" onSubmit={submit}>
              {variants.map((variant, index) => (
                <div key={index} className="rounded-2xl border border-white/10 bg-ink/45 p-4">
                  <div className="mb-3 flex items-center justify-between">
                    <span className="text-xs font-medium text-emerald-50/45">VARIANT {String(index + 1).padStart(2, "0")}</span>
                    {variants.length > 1 && (
                      <button aria-label={`Remove variant ${index + 1}`} type="button" onClick={() => removeVariant(index)} className="rounded-lg p-1.5 text-emerald-50/35 transition hover:bg-white/5 hover:text-red-300"><Trash2 size={15} /></button>
                    )}
                  </div>
                  <div className="grid grid-cols-2 gap-3 sm:grid-cols-4">
                    <Field label="Chromosome" value={variant.chrom} placeholder="chr17" onChange={(v) => updateVariant(index, "chrom", v)} />
                    <Field label="Position" value={variant.pos} placeholder="43000000" inputMode="numeric" onChange={(v) => updateVariant(index, "pos", v)} />
                    <Field label="Reference" value={variant.ref} placeholder="A" maxLength={50} onChange={(v) => updateVariant(index, "ref", v)} />
                    <Field label="Alternate" value={variant.alt} placeholder="G" maxLength={50} onChange={(v) => updateVariant(index, "alt", v)} />
                  </div>
                  <div className="mt-3">
                    <Field label="Gene (optional)" value={variant.gene} placeholder="Gene symbol, if known" onChange={(v) => updateVariant(index, "gene", v)} />
                  </div>
                </div>
              ))}

              <div className="flex flex-wrap items-center justify-between gap-3">
                <button type="button" onClick={addVariant} disabled={variants.length >= 25} className="inline-flex items-center gap-2 rounded-xl px-3 py-2 text-sm text-emerald-50/65 transition hover:bg-white/5 hover:text-mint disabled:opacity-30">
                  <Plus size={15} /> Add another <span className="text-xs text-emerald-50/35">({variants.length}/25)</span>
                </button>
                <button type="submit" disabled={busy || STATIC_PREVIEW} className="inline-flex min-w-40 items-center justify-center gap-2 rounded-xl bg-mint px-5 py-3 text-sm font-semibold text-ink transition hover:bg-[#b8f0d5] disabled:cursor-not-allowed disabled:opacity-70">
                  {busy ? <LoaderCircle size={16} className="animate-spin" /> : <ArrowDownRight size={16} />}
                  {STATIC_PREVIEW ? "Preview only" : busy ? "Submitting…" : "Analyse variants"}
                </button>
              </div>
            </form>

            {notice && (
              <div role="alert" className={`mt-5 flex items-start gap-3 rounded-xl border p-4 text-sm leading-6 ${notice.type === "error" ? "border-amber-300/20 bg-amber-300/[0.06] text-amber-100/80" : "border-aqua/20 bg-aqua/[0.06] text-aqua"}`}>
                {notice.type === "error" ? <TriangleAlert size={18} className="mt-0.5 shrink-0" /> : <CheckCircle2 size={18} className="mt-0.5 shrink-0" />}
                <span>{notice.message}</span>
                <button onClick={() => setNotice(null)} aria-label="Dismiss message" className="ml-auto rounded p-1 opacity-65 hover:opacity-100"><X size={15} /></button>
              </div>
            )}

            {apiUnavailable && (
              <div className="mt-5 flex items-start gap-3 rounded-xl border border-white/10 bg-white/[0.025] p-4 text-sm leading-6 text-emerald-50/55">
                <CircleHelp size={18} className="mt-0.5 shrink-0 text-mint/65" />
                <span>API not detected. Start FastAPI on port 8000 to check service and model status. Live prediction is intentionally not enabled in this prototype.</span>
              </div>
            )}
            {health && (
              <div className="mt-5 flex items-center gap-2 text-xs text-emerald-50/45">
                <span className={`h-1.5 w-1.5 rounded-full ${health.status === "ok" ? "bg-aqua" : "bg-amber-300"}`} />
                API {health.status} · {health.model_version} · {health.live_provider_ready ? "prediction provider ready" : "live provider not configured"}
              </div>
            )}

            {results.length > 0 && (
              <div className="mt-7 overflow-hidden rounded-2xl border border-white/10">
                <div className="border-b border-white/10 bg-amber-100/[0.06] p-4">
                  <div className="flex items-start gap-3">
                    <TriangleAlert size={18} className="mt-0.5 shrink-0 text-amber-200/90" />
                    <div>
                      <h4 className="text-sm font-semibold text-amber-50">Research use only — not for clinical diagnosis or treatment</h4>
                      <p className="mt-1 text-xs leading-5 text-amber-50/65">
                        These scores are produced by a frozen research prototype and are calibrated to the ClinVar prevalence used during model development. They may not transfer to populations or diseases with different prevalence. Do not interpret them as a diagnosis, prognosis, or treatment recommendation.
                      </p>
                    </div>
                  </div>
                </div>
                <div className="flex flex-col justify-between gap-3 border-b border-white/10 bg-white/[0.025] p-4 sm:flex-row sm:items-center">
                  <div>
                    <h3 className="font-semibold text-white">Analysis results</h3>
                    <p className="mt-1 text-xs text-emerald-50/45">
                      {resultMeta?.model_version ?? "AGVarPred"} · Research use only
                    </p>
                  </div>
                  <button type="button" onClick={() => downloadCsv(results)} className="inline-flex items-center justify-center gap-2 rounded-lg border border-white/15 px-3 py-2 text-xs font-medium text-emerald-50/70 transition hover:border-mint/40 hover:text-mint">
                    <FileUp size={14} /> Download CSV
                  </button>
                </div>
                <div className="max-h-[420px] overflow-auto">
                  <table className="w-full min-w-[800px] border-collapse text-left text-xs">
                    <thead className="sticky top-0 bg-[#153331] text-emerald-50/50">
                      <tr>
                        {["Variant", "Gene", "AGVarPred score", "Threshold class", "gnomAD AF", "Consequence / status"].map((label) => (
                          <th key={label} className="px-4 py-3 font-medium">{label}</th>
                        ))}
                      </tr>
                    </thead>
                    <tbody>
                      {results.map((result) => (
                        <tr key={result.variant_id} className="border-t border-white/[0.07]">
                          <td className="px-4 py-3 font-mono text-emerald-50/75">{result.variant_id}</td>
                          <td className="px-4 py-3 text-emerald-50/60">{result.gene || "—"}</td>
                          <td className="px-4 py-3 font-semibold text-white">{result.status === "ok" && typeof result.score === "number" ? result.score.toFixed(4) : "—"}</td>
                          <td className="px-4 py-3">
                            {result.status === "ok" ? (
                              <span className={`rounded-full px-2 py-1 ${result.predicted_class === "Pathogenic" ? "bg-amber-200/10 text-amber-100" : "bg-aqua/10 text-aqua"}`}>{result.predicted_class}</span>
                            ) : <span className="text-amber-200">Not scored</span>}
                          </td>
                          <td className="px-4 py-3 font-mono text-emerald-50/60">{typeof result.gnomAD_AF === "number" ? result.gnomAD_AF.toExponential(2) : "—"}</td>
                          <td className="max-w-xs px-4 py-3 text-emerald-50/55">{result.status === "error" ? result.error : result.consequence || "Annotation unavailable"}</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
                {resultMeta?.notice && <p className="border-t border-white/10 px-4 py-3 text-xs leading-5 text-emerald-50/45">{resultMeta.notice}</p>}
              </div>
            )}
          </section>
        </section>

        <section className="mt-14 grid gap-4 sm:grid-cols-3">
          <InfoCard number="01" title="Functional signals" text="AlphaGenome-derived genomic activity and splicing tracks contribute functional context." />
          <InfoCard number="02" title="Variant annotations" text="VEP-derived consequence information and population allele frequency complement the functional features." />
          <InfoCard number="03" title="Calibrated score" text="The frozen model returns a ClinVar-calibrated score; calibration may not transfer to populations with different prevalence." />
        </section>

        <section className="mt-12 flex flex-col justify-between gap-4 rounded-2xl border border-amber-200/15 bg-amber-100/[0.035] p-5 sm:flex-row sm:items-center sm:p-6">
          <div className="flex items-start gap-3">
            <TriangleAlert size={19} className="mt-0.5 shrink-0 text-amber-200/80" />
            <div>
              <h3 className="text-sm font-semibold text-amber-50">Research prototype — not for clinical decisions</h3>
              <p className="mt-1 max-w-3xl text-sm leading-6 text-amber-50/55">
                Predictions are not a clinical diagnosis, prognosis, or treatment recommendation. The operating threshold is empirical, scores are calibrated to the ClinVar prevalence used during model development, and live AlphaGenome-derived features may differ from the historical training snapshot. Do not use these results for patient care without independent expert review.
              </p>
            </div>
          </div>
          <a href="https://doi.org/10.5281/zenodo.21997464" target="_blank" rel="noreferrer" className="inline-flex shrink-0 items-center gap-2 text-sm text-amber-100/70 hover:text-amber-50">
            Archive <ArrowRight size={14} />
          </a>
        </section>

        <footer className="mt-10 flex flex-col justify-between gap-3 border-t border-white/10 pt-6 text-xs text-emerald-50/35 sm:flex-row sm:items-center">
          <span>AGVarPred · Open research software prototype</span>
        </footer>
      </main>
    </div>
  );
}

function downloadCsv(results: PredictionResult[]) {
  const columns: (keyof PredictionResult)[] = [
    "variant_id", "chrom", "pos", "ref", "alt", "gene", "status",
    "score", "predicted_class", "threshold", "gnomAD_AF", "consequence", "error",
  ];
  const escapeCell = (value: unknown) => {
    const text = value === null || value === undefined ? "" : String(value);
    return /[",\r\n]/.test(text) ? `"${text.replace(/"/g, '""')}"` : text;
  };
  const csv = [
    columns.join(","),
    ...results.map((row) => columns.map((column) => escapeCell(row[column])).join(",")),
  ].join("\r\n");
  const blob = new Blob([csv], { type: "text/csv;charset=utf-8;" });
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = "agvarpred_results.csv";
  link.click();
  URL.revokeObjectURL(url);
}

function Field({
  label, value, placeholder, onChange, inputMode, maxLength,
}: {
  label: string; value: string; placeholder: string; onChange: (v: string) => void;
  inputMode?: "numeric"; maxLength?: number;
}) {
  return (
    <label className="block min-w-0">
      <span className="mb-1.5 block text-xs text-emerald-50/50">{label}</span>
      <input
        className="w-full min-w-0 rounded-lg border border-white/10 bg-white/[0.035] px-3 py-2.5 text-sm text-white outline-none transition placeholder:text-emerald-50/20 focus:border-mint/50 focus:ring-2 focus:ring-mint/10"
        value={value}
        placeholder={placeholder}
        onChange={(event) => onChange(event.target.value)}
        inputMode={inputMode}
        maxLength={maxLength}
        autoCapitalize="off"
        autoComplete="off"
        spellCheck={false}
      />
    </label>
  );
}

function Metric({ label, value }: { label: string; value: string }) {
  return (
    <div className="rounded-xl border border-white/10 bg-white/[0.025] px-3 py-4">
      <div className="text-xl font-semibold tracking-tight text-white sm:text-2xl">{value}</div>
      <div className="mt-1 text-xs leading-5 text-emerald-50/45">{label}</div>
    </div>
  );
}

function InfoCard({ number, title, text }: { number: string; title: string; text: string }) {
  return (
    <div className="rounded-2xl border border-white/10 bg-panel/70 p-5 transition hover:border-mint/20">
      <div className="text-xs font-semibold tracking-[0.18em] text-mint/60">{number} / PLATFORM</div>
      <h3 className="mt-4 text-lg font-semibold text-white">{title}</h3>
      <p className="mt-2 text-sm leading-6 text-emerald-50/50">{text}</p>
    </div>
  );
}

export default App;
