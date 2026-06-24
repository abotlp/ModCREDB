(function () {
  "use strict";

  const byId = (id) => document.getElementById(id);
  const num = (value) => (Number.isFinite(Number(value)) ? Number(value) : 0);
  const X = { root: 90, branch: 330, leaf: 650 };
  const TOP = 48;
  const ROW = 70;
  const RIGHT_PAD = 520;

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

  function tfCount(node) { return num(node.tf_count || node.count); }
  function countLabel(node) {
    const count = tfCount(node);
    return count ? `${count.toLocaleString()} records` : "count pending";
  }
  function childMap(nodes) {
    const out = new Map();
    nodes.forEach((node) => {
      const parent = node.parent_id || null;
      if (!out.has(parent)) out.set(parent, []);
      out.get(parent).push(node);
    });
    out.forEach((children) => children.sort((a, b) => num(a.display_order) - num(b.display_order) || String(a.label).localeCompare(String(b.label))));
    return out;
  }

  function computeLayout(nodes) {
    const children = childMap(nodes);
    const byFamily = new Map(nodes.map((node) => [node.family_id, node]));
    let y = TOP;

    function layoutNode(node, depth) {
      const nodeChildren = children.get(node.family_id) || [];
      node._depth = depth;
      node._x = depth === 0 ? X.root : depth === 1 ? X.branch : X.leaf;
      if (!nodeChildren.length) {
        node._y = y;
        y += ROW;
        return node._y;
      }
      const childYs = nodeChildren.map((child) => layoutNode(child, depth + 1));
      node._y = (Math.min(...childYs) + Math.max(...childYs)) / 2;
      return node._y;
    }

    const roots = children.get(null).length ? children.get(null) : [byFamily.get("root") || nodes[0]];
    roots.filter(Boolean).forEach((root) => layoutNode(root, 0));
    return { byFamily, height: Math.max(700, y + TOP), width: X.leaf + RIGHT_PAD };
  }

  function radius(node, maxCount) {
    if (node.family_id === "root") return 38;
    const count = tfCount(node);
    if (!count || !maxCount) return 9;
    return Math.max(8, Math.min(34, 7 + Math.sqrt(count / maxCount) * 31));
  }
  function innerRadius(node, outer) {
    const total = tfCount(node);
    const pwm = num(node.generated_pwm_tf_count);
    if (!total || !pwm) return Math.max(2, outer * 0.30);
    return outer * Math.sqrt(Math.min(pwm / total, 1));
  }
  function pathD(parent, child) {
    const mx = (parent._x + child._x) / 2;
    return `M${parent._x},${parent._y} C${mx},${parent._y} ${mx},${child._y} ${child._x},${child._y}`;
  }
  function setText(id, value) {
    const node = byId(id);
    if (node) node.textContent = value;
  }
  function setMetric(id, value) { setText(id, num(value).toLocaleString()); }
  function updatePanel(node) {
    setText("selected-family-eyebrow", node.parent_id ? "Selected family" : "Database root");
    setText("selected-family-label", node.label || node.short_label || "Family");
    setMetric("selected-family-count", tfCount(node));
    setMetric("selected-known-count", node.known_count);
    setMetric("selected-homologous-count", node.homologous_count);
    setMetric("selected-relative-count", node.relative_homologous_count);
    setMetric("selected-predicted-count", node.predicted_low_count || (num(node.modcre_count) + num(node.alphafold_count)));
    setMetric("selected-generated-count", node.generated_pwm_tf_count);
    setMetric("selected-model-count", node.model_tf_count || node.active_model_tf_count);
    setMetric("selected-monomer-count", node.monomer_model_tf_count);
    setMetric("selected-dimer-count", node.dimer_model_tf_count);
    const open = byId("selected-family-open");
    if (open) open.href = node.open_url || (node.search_query ? `/search?q=${encodeURIComponent(node.search_query)}` : "/search");
    const search = byId("selected-family-search");
    if (search) {
      const q = node.search_query || node.short_label || node.label || "";
      search.href = q ? `/search?q=${encodeURIComponent(q)}` : "/search";
      search.textContent = q ? `Search ${q}` : "Search records";
    }
  }
  function addText(svg, text, x, y, className, anchor) {
    const element = document.createElementNS("http://www.w3.org/2000/svg", "text");
    element.setAttribute("x", x);
    element.setAttribute("y", y);
    element.setAttribute("class", className);
    if (anchor) element.setAttribute("text-anchor", anchor);
    element.textContent = text;
    svg.appendChild(element);
  }
  function select(nodes, familyId) {
    const selected = nodes.find((node) => node.family_id === familyId);
    if (!selected) return;
    document.querySelectorAll(".family-tree-node").forEach((node) => {
      node.classList.toggle("is-selected", node.getAttribute("data-family-id") === familyId);
    });
    updatePanel(selected);
  }
  function render(nodes) {
    const svg = byId("family-tree-svg");
    if (!svg || !nodes.length) return;
    svg.textContent = "";
    const { byFamily, height, width } = computeLayout(nodes);
    svg.setAttribute("viewBox", `0 0 ${width} ${height}`);
    svg.style.height = `${height}px`;
    svg.style.width = `${width}px`;
    const maxCount = Math.max(...nodes.map(tfCount), 1);
    const tooltip = byId("family-tree-tooltip");

    nodes.forEach((node) => {
      const parent = byFamily.get(node.parent_id);
      if (!parent) return;
      const path = document.createElementNS("http://www.w3.org/2000/svg", "path");
      path.setAttribute("class", "family-tree-link");
      path.setAttribute("d", pathD(parent, node));
      svg.appendChild(path);
    });

    nodes.forEach((node) => {
      const group = document.createElementNS("http://www.w3.org/2000/svg", "g");
      group.setAttribute("class", "family-tree-node");
      group.setAttribute("data-family-id", node.family_id);
      group.setAttribute("tabindex", "0");
      group.setAttribute("role", "button");
      group.setAttribute("transform", `translate(${node._x},${node._y})`);
      const outerRadius = radius(node, maxCount);
      const outer = document.createElementNS("http://www.w3.org/2000/svg", "circle");
      outer.setAttribute("class", "family-tree-node-ring");
      outer.setAttribute("r", outerRadius);
      group.appendChild(outer);
      const inner = document.createElementNS("http://www.w3.org/2000/svg", "circle");
      inner.setAttribute("class", "family-tree-node-fill");
      inner.setAttribute("r", innerRadius(node, outerRadius));
      group.appendChild(inner);
      group.addEventListener("click", () => select(nodes, node.family_id));
      group.addEventListener("keydown", (event) => {
        if (event.key === "Enter" || event.key === " ") { event.preventDefault(); select(nodes, node.family_id); }
      });
      group.addEventListener("mousemove", (event) => {
        if (!tooltip) return;
        tooltip.hidden = false;
        tooltip.style.left = `${event.offsetX + 16}px`;
        tooltip.style.top = `${event.offsetY + 16}px`;
        tooltip.innerHTML = `<strong>${node.label}</strong><br>${countLabel(node)}<br>Generated PWM TFs: ${num(node.generated_pwm_tf_count).toLocaleString()}<br>3D model TFs: ${num(node.model_tf_count).toLocaleString()}`;
      });
      group.addEventListener("mouseleave", () => { if (tooltip) tooltip.hidden = true; });
      svg.appendChild(group);
      const anchor = node._depth < 2 ? "end" : "start";
      const labelX = node._x + (anchor === "end" ? -(outerRadius + 10) : outerRadius + 12);
      addText(svg, node.short_label || node.label, labelX, node._y - 2, "family-tree-node-label", anchor);
      addText(svg, countLabel(node), labelX, node._y + 15, "family-tree-node-count", anchor);
    });
    const first = nodes.find((node) => node.family_id === "2.3") || nodes.find((node) => tfCount(node) > 0) || nodes[0];
    if (first) select(nodes, first.family_id);
  }

  document.addEventListener("DOMContentLoaded", async () => render(await loadTreeData()));
})();
