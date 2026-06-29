(() => {
  const pairs = [
    ["Predicted Low: AlphaFold3", "AF3-derived ModCRE"],
    ["AlphaFold3-assisted ModCRE", "AF3-derived ModCRE"],
    ["AlphaFold_ModCRE", "AF3-derived ModCRE"],
    ["Predicted Low: ModCRE", "ModCRE"],
    ["Relatively_Homologous_PWM", "Nearest Neighbor (50–70%)"],
    ["Relatively Homologous PWM", "Nearest Neighbor (50–70%)"],
    ["Distant homologous candidate", "Nearest Neighbor (50–70%)"],
    ["Homologous_PWM", "Nearest Neighbor (>70%)"],
    ["Homologous PWM", "Nearest Neighbor (>70%)"],
    ["Close homologous PWM", "Nearest Neighbor (>70%)"],
    ["Identical_PWM", "Known"],
    ["Identical PWM", "Known"],
    ["Direct PWM", "Known"],
    ["FIMO-ready", "Generated PWM"],
    ["w=0 / no matrix", "Missing MEME"],
    ["Active model files", "3D Models"],
    ["active PDB models", "3D models"],
    ["active PDB model", "3D model"],
    ["active models", "3D models"],
    ["TF links", "TFs using this motif"],
    ["Motif links", "TFs using this motif"]
  ];
  function replaceText(root) {
    const walker = document.createTreeWalker(root, NodeFilter.SHOW_TEXT);
    const nodes = [];
    while (walker.nextNode()) nodes.push(walker.currentNode);
    nodes.forEach(node => {
      if (node.parentElement && ["SCRIPT", "STYLE", "TEXTAREA"].includes(node.parentElement.tagName)) return;
      let text = node.nodeValue;
      pairs.forEach(pair => { text = text.split(pair[0]).join(pair[1]); });
      const trimmed = text.trim();
      if (trimmed === "Families") text = text.replace("Families", "PFAM Families");
      if (trimmed === "Motifs") text = text.replace("Motifs", "Generated PWM");
      node.nodeValue = text;
    });
  }
  replaceText(document.body);
})();
