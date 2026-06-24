(function () {
  "use strict";

  function byId(id) {
    return document.getElementById(id);
  }

  function asNumber(value, fallback) {
    const number = Number(value);
    return Number.isFinite(number) ? number : fallback;
  }

  function escapeText(value) {
    return String(value == null ? "" : value);
  }

  function readEmbeddedTree() {
    const node = byId("family-tree-data");
    if (!node) {
      return [];
    }
    try {
      return JSON.parse(node.textContent || "[]");
    } catch (error) {
      console.error("Could not parse embedded family tree data", error);
      return [];
    }
  }

  async function loadTreeData() {
    const fallback = readEmbeddedTree();
    try {
      const response = await fetch("/static/family_tree_data.generated.json", { cache: "no-store" });
      if (!response.ok) {
        return fallback;
      }
      const generated = await response.json();
      if (Array.isArray(generated) && generated.length) {
        return generated;
      }
      if (generated && Array.isArray(generated.nodes) && generated.nodes.length) {
        return generated.nodes;
      }
    } catch (error) {
      // Generated DB counts are optional. Fall back to embedded taxonomy/counts.
    }
    return fallback;
  }

  function descendantsByParent(nodes) {
    const byParent = new Map();
    nodes.forEach((node) => {
      if (!node.parent_id) {
        return;
      }
      if (!byParent.has(node.parent_id)) {
        byParent.set(node.parent_id, []);
      }
      byParent.get(node.parent_id).push(node);
    });
    return byParent;
  }

  function familyCount(node) {
    return asNumber(node.tf_count, asNumber(node.count, 0));
  }

  function nodeRadius(node, maxCount) {
    const count = familyCount(node);
    if (node.family_id === "root") {
      return 43;
    }
    if (!maxCount || count <= 0) {
      return node.children_count ? 18 : 10;
    }
    return Math.max(10, Math.min(40, 8 + Math.sqrt(count / maxCount) * 36));
  }

  function generatedRadius(node, outerRadius) {
    const count = asNumber(node.generated_pwm_tf_count, 0);
    const total = familyCount(node);
    if (!count || !total) {
      return Math.max(0, outerRadius * 0.35);
    }
    return Math.max(2, outerRadius * Math.sqrt(Math.min(count / total, 1)));
  }

  function pathD(parent, child) {
    const px = asNumber(parent.x, 0);
    const py = asNumber(parent.y, 0);
    const cx = asNumber(child.x, 0);
    const cy = asNumber(child.y, 0);
    const mx = (px + cx) / 2;
    return `M${px},${py} C${mx},${py} ${mx},${cy} ${cx},${cy}`;
  }

  function countLabel(node) {
    const count = familyCount(node);
    if (count > 0) {
      return `${count.toLocaleString()} records`;
    }
    if (node.placeholder_count) {
      return node.placeholder_count;
    }
    return "count pending";
  }

  function setText(id, value) {
    const node = byId(id);
    if (node) {
      node.textContent = value;
    }
  }

  function setMetric(id, value) {
    const node = byId(id);
    if (!node) {
      return;
    }
    const number = Number(value || 0);
    node.textContent = Number.isFinite(number) ? number.toLocaleString() : String(value || "0");
  }

  function updateSelectedPanel(node) {
    setText("selected-family-eyebrow", node.parent_id ? "Selected family" : "Database root");
    setText("selected-family-label", node.label || node.short_label || "Family");
    setMetric("selected-family-count", familyCount(node));
    setText("selected-family-count-label", node.family_id === "root" ? "TF sequence records in ModCREDB" : "TF sequence records");

    setMetric("selected-known-count", node.known_count || node.identical_count || 0);
    setMetric("selected-homologous-count", node.homologous_count || 0);
    setMetric("selected-relative-count", node.relative_homologous_count || 0);
    setMetric("selected-predicted-count", node.predicted_low_count || (asNumber(node.modcre_count, 0) + asNumber(node.alphafold_count, 0)));
    setMetric("selected-generated-count", node.generated_pwm_tf_count || 0);
    setMetric("selected-model-count", node.model_tf_count || node.active_model_tf_count || 0);
    setMetric("selected-monomer-count", node.monomer_model_tf_count || 0);
    setMetric("selected-dimer-count", node.dimer_model_tf_count || 0);

    const open = byId("selected-family-open");
    if (open) {
      open.href = node.open_url || (node.search_query ? `/search?q=${encodeURIComponent(node.search_query)}` : "/search");
    }
    const search = byId("selected-family-search");
    if (search) {
      const query = node.search_query || node.short_label || node.label || "";
      search.href = query ? `/search?q=${encodeURIComponent(query)}` : "/search";
      search.textContent = query ? `Search ${query}` : "Search records";
    }
  }

  function addSvgText(svg, text, x, y, className, anchor) {
    const element = document.createElementNS("http://www.w3.org/2000/svg", "text");
    element.setAttribute("x", x);
    element.setAttribute("y", y);
    element.setAttribute("class", className);
    if (anchor) {
      element.setAttribute("text-anchor", anchor);
    }
    element.textContent = text;
    svg.appendChild(element);
    return element;
  }

  function renderTree(nodes) {
    const svg = byId("family-tree-svg");
    if (!svg || !nodes.length) {
      return;
    }
    svg.textContent = "";

    const byFamily = new Map(nodes.map((node) => [node.family_id, node]));
    const maxCount = Math.max(...nodes.map(familyCount), 1);
    const tooltip = byId("family-tree-tooltip");

    nodes.forEach((node) => {
      if (!node.parent_id) {
        return;
      }
      const parent = byFamily.get(node.parent_id);
      if (!parent) {
        return;
      }
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
      group.setAttribute("aria-label", `${node.label}: ${countLabel(node)}`);
      group.setAttribute("transform", `translate(${asNumber(node.x, 0)},${asNumber(node.y, 0)})`);

      const outerRadius = nodeRadius(node, maxCount);
      const outer = document.createElementNS("http://www.w3.org/2000/svg", "circle");
      outer.setAttribute("class", "family-tree-node-ring");
      outer.setAttribute("r", outerRadius);
      group.appendChild(outer);

      const inner = document.createElementNS("http://www.w3.org/2000/svg", "circle");
      inner.setAttribute("class", "family-tree-node-fill");
      inner.setAttribute("r", generatedRadius(node, outerRadius));
      group.appendChild(inner);

      group.addEventListener("click", () => selectNode(nodes, node.family_id));
      group.addEventListener("keydown", (event) => {
        if (event.key === "Enter" || event.key === " ") {
          event.preventDefault();
          selectNode(nodes, node.family_id);
        }
      });
      group.addEventListener("mousemove", (event) => {
        if (!tooltip) {
          return;
        }
        tooltip.hidden = false;
        tooltip.style.left = `${event.offsetX + 16}px`;
        tooltip.style.top = `${event.offsetY + 16}px`;
        tooltip.innerHTML = `<strong>${escapeText(node.label)}</strong><br>${escapeText(countLabel(node))}<br>Generated PWM TFs: ${asNumber(node.generated_pwm_tf_count, 0).toLocaleString()}<br>3D model TFs: ${asNumber(node.model_tf_count || node.active_model_tf_count, 0).toLocaleString()}`;
      });
      group.addEventListener("mouseleave", () => {
        if (tooltip) {
          tooltip.hidden = true;
        }
      });

      svg.appendChild(group);

      const labelAnchor = asNumber(node.x, 0) < 410 ? "end" : "start";
      const labelX = asNumber(node.x, 0) + (labelAnchor === "end" ? -(outerRadius + 8) : outerRadius + 8);
      const labelY = asNumber(node.y, 0) + 4;
      addSvgText(svg, node.short_label || node.label, labelX, labelY, "family-tree-node-label", labelAnchor);
      addSvgText(svg, countLabel(node), labelX, labelY + 15, "family-tree-node-count", labelAnchor);
    });

    const first = nodes.find((node) => node.family_id === "2.3") || nodes.find((node) => familyCount(node) > 0) || nodes[0];
    if (first) {
      selectNode(nodes, first.family_id);
    }
  }

  function selectNode(nodes, familyId) {
    const selected = nodes.find((node) => node.family_id === familyId);
    if (!selected) {
      return;
    }
    document.querySelectorAll(".family-tree-node").forEach((node) => {
      node.classList.toggle("is-selected", node.getAttribute("data-family-id") === familyId);
    });
    updateSelectedPanel(selected);
  }

  document.addEventListener("DOMContentLoaded", async () => {
    const nodes = await loadTreeData();
    renderTree(nodes);
  });
})();
