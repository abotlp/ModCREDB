(function () {
  "use strict";

  const byId = (id) => document.getElementById(id);
  const num = (value) => (Number.isFinite(Number(value)) ? Number(value) : 0);
  let nodes = [];
  let nodeById = new Map();
  let selectedFamilyId = null;

  const childPriority = {
    "1": ["1.1", "1.2", "1.3"],
    "2": ["2.3", "2.1", "2.2", "2.5", "2.6", "2.7", "2.8", "2.9"],
    "3": ["3.1", "3.3", "3.2", "3.4", "3.5", "3.6", "3.7"],
    "4": ["4.1", "4.2"],
    "5": ["5.1", "5.3"],
    "6": ["6.2", "6.3", "6.1", "6.5", "6.4", "6.6", "6.7"],
    "7": ["7.1", "7.2"],
    "8": ["8.1", "8.2"],
    "9": ["9.1"],
    "0": ["0.0", "0.5", "0.4", "0.1", "0.2", "0.3", "0.6"]
  };

  function readEmbeddedTree() {
    const node = byId("family-tree-data");
    if (!node) return [];
    try { return JSON.parse(node.textContent || "[]"); }
    catch (error) { console.error("Could not parse family-tree data", error); return []; }
  }

  async function loadTreeData() {
    const fallback = readEmbeddedTree();
    try {
      const response = await fetch("/static/family_tree_data.generated.json", { cache: "no-store" });
      if (!response.ok) return fallback;
      const data = await response.json();
      if (Array.isArray(data.nodes) && data.nodes.length) return data.nodes;
      if (Array.isArray(data) && data.length) return data;
    } catch (error) {}
    return fallback;
  }

  function tfCount(node) { return num(node && (node.tf_count || node.count)); }
  function fmt(value) { return num(value).toLocaleString(); }
  function searchUrl(node) {
    const q = node && (node.search_query || node.short_label || node.label || "");
    return q ? `/search?q=${encodeURIComponent(q)}` : "/search";
  }
  function displayLabel(node) {
    if (!node) return "TF family";
    return node.short_label || String(node.label || "TF family")
      .replace(" DNA-binding domains", "")
      .replace(" domains", "")
      .replace(" factors", "");
  }
  function childrenOf(parentId) {
    const raw = nodes.filter((node) => node.parent_id === parentId);
    const priority = childPriority[parentId] || [];
    return raw.sort((a, b) => {
      const ia = priority.indexOf(a.family_id);
      const ib = priority.indexOf(b.family_id);
      if (ia !== -1 || ib !== -1) return (ia === -1 ? 999 : ia) - (ib === -1 ? 999 : ib);
      const ca = tfCount(a);
      const cb = tfCount(b);
      if (ca || cb) return cb - ca;
      return num(a.display_order) - num(b.display_order) || displayLabel(a).localeCompare(displayLabel(b));
    });
  }
  function setText(id, value) {
    const element = byId(id);
    if (element) element.textContent = value;
  }
  function showHome() {
    const home = byId("home-view");
    const browse = byId("family-browse-view");
    if (home) home.classList.remove("hidden");
    if (browse) browse.classList.add("hidden");
    window.scrollTo({ top: 0, behavior: "smooth" });
  }
  function showBrowse() {
    const home = byId("home-view");
    const browse = byId("family-browse-view");
    if (home) home.classList.add("hidden");
    if (browse) browse.classList.remove("hidden");
    window.scrollTo({ top: 0, behavior: "smooth" });
  }
  function metricSpan(label, value) {
    const count = num(value);
    if (!count) return "";
    return `<span>${label} ${fmt(count)}</span>`;
  }
  function updateSelectedSubfamily(node) {
    if (!node) return;
    setText("selected-subfamily-title", displayLabel(node));
    const metrics = [
      metricSpan("Records", tfCount(node)),
      metricSpan("Known", node.known_count),
      metricSpan("Nearest Neighbor", node.homologous_count),
      metricSpan("50-70%", node.relative_homologous_count),
      metricSpan("Generated PWM", node.generated_pwm_tf_count),
      metricSpan("3D models", node.model_tf_count || node.active_model_tf_count)
    ].filter(Boolean).join("") || "<span>Open matching records in search</span>";
    const metricsNode = byId("selected-subfamily-metrics");
    if (metricsNode) metricsNode.innerHTML = metrics;
    const open = byId("selected-subfamily-open");
    if (open) open.href = searchUrl(node);
    const input = byId("family-browse-q");
    if (input) input.value = node.search_query || node.short_label || node.label || "";
    document.querySelectorAll(".subfamily-card").forEach((card) => {
      card.classList.toggle("is-selected", card.getAttribute("data-family-id") === node.family_id);
    });
  }
  function renderSubfamilies(parentId) {
    const grid = byId("subfamily-grid");
    if (!grid) return;
    grid.textContent = "";
    const children = childrenOf(parentId);
    children.forEach((child) => {
      const button = document.createElement("button");
      button.type = "button";
      button.className = "subfamily-card";
      button.setAttribute("data-family-id", child.family_id);
      const count = tfCount(child);
      button.innerHTML = `<strong>${displayLabel(child)}</strong><span class="badge">${count ? fmt(count) : "Open"}</span>`;
      button.addEventListener("click", () => updateSelectedSubfamily(child));
      grid.appendChild(button);
    });
    updateSelectedSubfamily(children[0] || nodeById.get(parentId));
  }
  function openFamily(parentId) {
    const parent = nodeById.get(parentId);
    if (!parent) return;
    selectedFamilyId = parentId;
    setText("family-browse-title", displayLabel(parent));
    setText("family-browse-description", "Choose a subfamily, then filter matching records.");
    const input = byId("family-browse-q");
    if (input) {
      input.value = parent.search_query || "";
      input.placeholder = `Search within ${displayLabel(parent)}`;
    }
    renderSubfamilies(parentId);
    showBrowse();
  }
  function bindFamilyCards() {
    document.querySelectorAll(".home-family-card[data-family-id]").forEach((card) => {
      card.addEventListener("click", () => openFamily(card.getAttribute("data-family-id")));
    });
    const back = byId("back-to-home");
    if (back) back.addEventListener("click", showHome);
  }

  document.addEventListener("DOMContentLoaded", async () => {
    nodes = await loadTreeData();
    nodeById = new Map(nodes.map((node) => [node.family_id, node]));
    bindFamilyCards();
    if (selectedFamilyId) openFamily(selectedFamilyId);
  });
})();
